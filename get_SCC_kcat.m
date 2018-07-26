function [SCC, Status, B_min, B_max, ReactionSet, ODE] = get_SCC_kcat(model,fctable,blk,kcat,method,met_to_be_checked,V)
%
% Preliminaries:
%
%                   1. all reactions should be irreversible
%                   2. no blocked reactions
%                   3. F2C2 toolbox or coupling matrix as input needed
%
%                      F2C2 toolbox (Larhlimi et al., BMC Bioinformatics (2012)) can be downloaded from
%                      https://sourceforge.net/projects/f2c2/files/
%
%
% Input: model: Model in Cobra format (unblocked, irreversible reactions)
%        fctable: coupling matrix obtained from F2C2 (optional)
%        blk: vector indicating blockd reactions obtained from F2C2 (optional)
%        kcat: table of known rate constants (format see Supplementary Table 3)
%        method: for unknown ratios take 'equal','average','median'
%        met_to_be_checked: index of metabolites to be checked for SCC
%        V: flux distribution to get ratio of fully coupled reactions
%        (optional) otherwise fractional LP is used
%
% Output:
%
%           SCC: vector of length m x 1
%                0 if approach is not applicable, 1 if species shows SCC
%           Status: vector of length m x 1
%                vector ndicating numerical problems during calculation
%                3 numerical problems detect, 0 otherwise
%           B_min/B_max: min and max of B_ps (Eq. 1 main text)
%           ReacktionSet: set of reactions around each SCC metabolite (v_p
%           and v_s Eq. 1 main text) 1. colum v_p, 2. colum v_s
%           ODE: ODE from which metabolite was found to be SCC
%

addpath(genpath('F2C2/'))

%% check input data
if ~all(model.rev==0)
    warning('the network contains reversible reactions')
end

% calculate coupling matrix if not given
if isempty(fctable) || isempty(blk)
    T=evalc('[fctable,blk]=F2C2(''glpk'',CobraToF2C2(model));');
end

if isempty(kcat)
    kcat=table(model.rxns, cell(size(model.rxns)), cell(size(model.rxns)), nan(size(model.rxns)));
end

Nplus = model.S;
Nplus(Nplus>0) = 0;

%%
% for method 'average' or 'median' take average/median ratio over kcat
% ratios from reactions sharing a substrate
set=[];
for i=1:size(Nplus,2)
    r = find(sum(abs(Nplus-repmat(Nplus(:,i),1,size(Nplus,2))))==1);
    if ~isempty(r)
        r_n = r(find(sum(Nplus(:,r)-repmat(Nplus(:,i),1,length(r)))==-1));
        if ~isempty(r_n)
            set(end+1:end+length(r_n),:) = [repmat(i,length(r_n),1) r_n'];
        end
    end
end

ratio=[];
for i=1:size(set,1)
    ratio(end+1) = mean(kcat{strcmp(kcat{:,1},model.rxns(set(i,1))),4})./mean(kcat{strcmp(kcat{:,1},model.rxns(set(i,2))),4});
end

AV = mean(ratio(~isnan(ratio)));
ME = median(ratio(~isnan(ratio)));

if ~all(blk==0)
    warning('blocked reactions detected during coupling matrix calculation')
end

%% define variables
Tol_Lower=model.lb;
Tol_Upper=model.ub;
Min_t=min(model.lb);
N = model.S;
num_mets = size(N,1);
num_rxns = size(N,2);
SCC=zeros(size(N,1),1);
B_min=cell(size(N,1),1);B_max=cell(size(N,1),1);
ODE=cell(size(N,1),1);
ReactionSet=cell(size(N,1),1);
Status=zeros(size(N,1),1);
Tol_t=1e+10;

% check minimum flux of all reactions
T=evalc('[vmin,vmax]=fluxVariability(model)');
bs = ones(size(model.rxns));
bs(find(vmin>1)) = -1;

if isempty(V)
    model_V=model;
    model_V.c(:)=1;
    Sol = optimizeCbModel(model_V);
    V = Sol.x;
end

%% start SCC detection
Nplus=N.*(N<0);

%   S_met_idx: index of metabolite S for which we check SCC
for o=1:length(met_to_be_checked)
    S_met_idx = met_to_be_checked(o);
    
    % find reactions where S_met_idx is on the substrate side
    R = find(N(S_met_idx,:)<0);
    
    % find all species occuring in reactions R
    [M,~] = find(N(:,R)~=0);
    M=unique(M);
    
    for P=1:length(M) % for each ode in which S appears ...
        
        % find all substrate complexes including M(P) => Cp
        Cp=unique(Nplus(:,find(Nplus(M(P),:)<0))','rows','stable')'; % each column corresponds to one complex
        
        % remove one molecule of S from complexes in Cp => Cp_S
        vector = zeros(num_mets,1);
        vector(S_met_idx) = 1;
        
        Cp_S=Cp + repmat(vector,1,size(Cp,2));
        
        % check that complex Cp_S is still in the network
        Cp_S_check=cell(1,size(Cp_S,2));
        
        for c = 1:size(Cp_S,2)
            % case 1: zero complex
            C1=Cp_S(:,c);
            if all(C1==0)
                Cp_S_check{c}=[find(all(N<=0)) find(all(N<=0))]; % find import export reactions
                if isempty([find(all(N<=0)) find(all(N<=0))])
                    Cp_S_check{c}=NaN;
                end
                % case 2: complex including -X
            elseif ~all(C1<=0)
                Cp_S_check{c}=NaN;
                % case 3: find complex C1
            elseif all(C1<=0)
                Cp_S_check{c} = find(all(Nplus==repmat(C1,1,num_rxns)));
                if isempty(find(all(Nplus==repmat(C1,1,num_rxns)))) % if C1 is not in the network
                    Cp_S_check{c}=NaN;
                end
            end
        end
        
        % find complexes of reactions where M(P) is on product side
        Cx=unique(Nplus(:,find(N(M(P),:)>0))','rows','stable')';
        
        % find all substrate complex of P - S => Cp_S
        vector = zeros(num_mets,1);
        vector(S_met_idx) = 1;
        
        Cx_S=Cx + repmat(vector,1,size(Cx,2));
        
        % check that complex is still in the network
        Cx_S_check=cell(1,size(Cx_S,2));
        for c = 1:size(Cx_S,2)
            % case 1: zero complex
            C1=Cx_S(:,c);
            if all(C1==0)
                Cx_S_check{c}=[find(all(N<=0)) find(all(N<=0))];
                if isempty([find(all(N<=0)) find(all(N<=0))])
                    Cx_S_check{c}=NaN;
                end
                % case 2: complex including -X
            elseif ~all(C1<=0)
                Cx_S_check{c}=NaN;
                % case 3: find complex C1
            elseif all(C1<=0)
                Cx_S_check{c} = find(all(Nplus==repmat(C1,1,num_rxns)));
                if isempty(find(all(Nplus==repmat(C1,1,num_rxns))))
                    Cx_S_check{c}=NaN;
                end
            end
        end
        
        % check coupling of complexes in Cp & Cx_S if they are all in the network
        if length(find(~isnan(cell2mat(Cx_S_check))))>0 && ~isempty(Cp) && length(find(isnan(cell2mat(Cx_S_check))))==0
            % first Cp - Cx_S
            substrate_check=[];
            product_minus_check=[];
            Result1=[];
            
            for i=1:size(Cp,2)
                for j=1:size(Cp,2)
                    % coupling inside Cp
                    R1 = find(all(Nplus==repmat(Cp(:,i),1,size(Nplus,2))));
                    R2 = find(all(Nplus==repmat(Cp(:,j),1,size(Nplus,2))));
                    
                    substrate_check(end+1)=sum(sum(fctable(R1,R2)==1))>0; % two complexes are fully coupled if at least one pair of reactions is fully coupled
                end
            end
            for i=1:size(Cx_S,2)
                for j=1:size(Cx_S,2)
                    % coupling inside Cx_S
                    R1 = find(all(Nplus==repmat(Cx_S(:,i),1,size(Nplus,2))));
                    R2 = find(all(Nplus==repmat(Cx_S(:,j),1,size(Nplus,2))));
                    
                    product_minus_check(end+1)=sum(sum(fctable(R1,R2)==1))>0;
                end
            end
            if all(substrate_check==1) && all(product_minus_check==1) && ~isempty(substrate_check) && ~isempty(product_minus_check)
                
                SCC(S_met_idx) = 1;
                
                if S_met_idx==M(P)
                    X_i_X_j(S_met_idx) = 1;
                end
                
                Sum_i=[]; Sum_j_min=[]; Sum_j_max=[];
                % calculate SCC
                
                in=[];out=[];B_temp_min=[];B_temp_max=[];ReactionSet_temp=[];
                
                % v_p / v_g(l)
                for i=1:size(Cp,2)
                    for j=1:size(Cx_S,2) % size of Cx - mapped
                        R1 = find(all(Nplus==repmat(Cp(:,i),1,size(Nplus,2))));
                        R2 = find(all(Nplus==repmat(Cx_S(:,j),1,size(Nplus,2))));
                        
                        for gammaj = 1:length(R1) % over p
                            in(end+1,1) = R1(gammaj);
                            out(end+1,1) = R2(1); % v_g(l) - one possible Q
                        end
                        
                    end
                end
                
                % lambda_ii - vk/vp
                for in_i = 1:length(in)
                    s_ii=[];
                    lambda_ii=[];
                    for j=1:size(Cp,2)
                        R2 = find(all(Nplus==repmat(Cp(:,j),1,size(Nplus,2))));
                        
                        if length(R2)>1 % if there are more than two reactions assigned to one complex
                            r=find(abs(N(M(P),R2))~=0);
                            r=r(1);
                        else
                            r=1;
                        end
                        
                        if isempty(V)
                            % fractional LP
                            [x,Tol]=fractional_LP(N,Tol_Lower,Tol_Upper,Min_t,Tol_t,bs,in(in_i),R2(r),-1);
                            
                            if Tol ~=1
                                Status(S_met_idx) = 3;
                            end
                            
                            lambda_ii(end+1) = x(R2(r))/x(in(in_i)); % vk/vp
                            %  lambda_ii(end+1) = model.lb(R2(r))/model.lb(in(in_i)); % vk/vp
                            
                        else
                            lambda_ii(end+1) = V(R2(r))/V(in(in_i));
                        end
                        
                        s_ii(end+1) = abs(N(M(P),R2(r))); % stoichiometry of k
                        
                    end
                    Sum_i(in_i) = sum(s_ii.*lambda_ii);
                end
                
                % lambda_jj
                for out_i = 1:length(out)
                    s_jj=[];
                    lambda_jj=[];k_j_min=[];k_j_max=[];
                    for j=1:size(Cx_S,2)
                        R2 = find(all(Nplus==repmat(Cx_S(:,j),1,size(Nplus,2))));
                        Rs = find(all(Nplus==repmat(Cx(:,j),1,size(Nplus,2)))); % stoichiometry taken from original set of reactions
                        
                        if length(R2)>1 % if there is more than one reaction associated to Cx_S
                            r=find(fctable(R2,out(out_i))); % take one of the coupled reactions
                            r=r(1);
                        else
                            r=1;
                        end
                        
                        for s=1:length(Rs) % over all l
                            
                            if isempty(V)
                                % fractional LP
                                [x,Tol]=fractional_LP(N,Tol_Lower,Tol_Upper,Min_t,Tol_t,bs,out(out_i),R2(r),-1);
                                
                                if Tol ~=1
                                    Status(S_met_idx) = 3;
                                end
                                
                                lambda_jj(end+1) = x(R2(r))/x(out(out_i)); %v_g(l)/v_g(s)
                                % lambda_jj(end+1) = model.lb(R2(r))/model.lb(out(out_i)); %v_g(l)/v_g(s)
                            else
                                lambda_jj(end+1) = V(R2(r))/V(out(out_i));
                            end
                            
                            k_Rs_temp = mean(kcat{strcmp(kcat{:,1},model.rxns(Rs(s))),4});
                            k_j_min(end+1) = min(k_Rs_temp./mean(kcat{strcmp(kcat{:,1},model.rxns(R2(r))),4})); % teta(l)/ teta_g(l)
                            k_j_max(end+1) = max(k_Rs_temp./mean(kcat{strcmp(kcat{:,1},model.rxns(R2(r))),4}));
                            if isnan(k_j_min(end))
                                if strcmp(method,'equal')
                                    k_j_min(end)=1;
                                    k_j_max(end)=1;
                                elseif strcmp(method,'average')
                                    k_j_min(end)=AV;
                                    k_j_max(end)=AV;
                                elseif strcmp(method,'median')
                                    k_j_min(end)=ME;
                                    k_j_max(end)=ME;
                                end
                            end
                            
                            s_jj(end+1) = abs(N(M(P),Rs(s)));
                        end
                        
                    end
                    Sum_j_min(out_i) = sum(s_jj.*lambda_jj.*k_j_min);
                    Sum_j_max(out_i) = sum(s_jj.*lambda_jj.*k_j_max);
                end
                
                B_temp_max = Sum_i./(Sum_j_min);
                B_temp_min = Sum_i./(Sum_j_max);
                ReactionSet_temp(end+1:end+size([in out],1),1:2) = [in out];
                
            end
            
        elseif length(find(~isnan(cell2mat(Cp_S_check))))>0 && ~isempty(Cx) && length(find(isnan(cell2mat(Cp_S_check))))==0
            % Cx - Cp_S
            substrate_minus_check=[];
            product_check=[];
            Result2=[];
            
            for i=1:size(Cx,2)
                for j=1:size(Cx,2)
                    % coupling inside Cx
                    R1 = find(all(Nplus==repmat(Cx(:,i),1,size(Nplus,2))));
                    R2 = find(all(Nplus==repmat(Cx(:,j),1,size(Nplus,2))));
                    
                    substrate_minus_check(end+1)=sum(sum(fctable(R1,R2)==1))>0;
                end
            end
            for i=1:size(Cp_S,2)
                for j=1:size(Cp_S,2)
                    % coupling inside Cp_S
                    R1 = find(all(Nplus==repmat(Cp_S(:,i),1,size(Nplus,2))));
                    R2 = find(all(Nplus==repmat(Cp_S(:,j),1,size(Nplus,2))));
                    
                    product_check(end+1)=sum(sum(fctable(R1,R2)==1))>0;
                end
            end
            
            if all(substrate_minus_check==1) && all(product_check==1) && ~isempty(substrate_minus_check) && ~isempty(product_check)
                
                SCC(S_met_idx) = 1;
                if S_met_idx==M(P)
                    X_i_X_j(S_met_idx) = 1;
                end
                
                Sum_i=[]; Sum_j_max=[]; Sum_j_min=[];
                % calculate SCC
                
                in=[];out=[];B_temp_min=[];B_temp_max=[];ReactionSet_temp=[];
                
                % v_in / v_out
                for i=1:size(Cx,2)
                    for j=1:size(Cp_S,2)
                        R1 = find(all(Nplus==repmat(Cx(:,i),1,size(Nplus,2))));
                        R2 = find(all(Nplus==repmat(Cp_S(:,j),1,size(Nplus,2))));
                        
                        for gammaj = 1:length(R1)
                            
                            in(end+1,1) = R1(gammaj);
                            out(end+1,1) = R2(1);
                            
                        end
                    end
                end
                
                % lambda_ii
                for in_i = 1:length(in)
                    s_ii=[];
                    lambda_ii=[];
                    for j=1:size(Cx,2)
                        R2 = find(all(Nplus==repmat(Cx(:,j),1,size(Nplus,2))));
                        
                        if length(R2)>1
                            r=find(abs(N(M(P),R2))~=0);
                            r=r(1);
                        else
                            r=1;
                        end
                        
                        if isempty(V)
                            % fractional LP
                            [x,Tol]=fractional_LP(N,Tol_Lower,Tol_Upper,Min_t,Tol_t,bs,in(in_i),R2(r),-1);
                            
                            if Tol ~=1
                                Status(S_met_idx) = 3;
                            end
                            
                            lambda_ii(end+1) = x(R2(r))/x(in(in_i)); % vk/vp
                            % lambda_ii(end+1) = model.lb(R2(r))/model.lb(in(in_i)); % vk/vp
                            
                        else
                            lambda_ii(end+1) = V(R2(r))/V(in(in_i));
                        end                      
                        
                        s_ii(end+1) = abs(N(M(P),R2(r)));
                        
                    end
                    Sum_i(in_i) = sum(s_ii.*lambda_ii);
                end
                
                % lambda_jj
                for out_i = 1:length(out)
                    s_jj=[];
                    lambda_jj=[];k_j_min=[];k_j_max=[];
                    for j=1:size(Cp_S,2)
                        R2 = find(all(Nplus==repmat(Cp_S(:,j),1,size(Nplus,2))));
                        Rs = find(all(Nplus==repmat(Cp(:,j),1,size(Nplus,2))));
                        
                        if length(R2)>1
                            r=find(fctable(R2,out(out_i)));
                            r=r(1);
                        else
                            r=1;
                        end
                        
                        for s=1:length(Rs)
                            
                            if isempty(V)
                                % fractional LP
                                [x,Tol]=fractional_LP(N,Tol_Lower,Tol_Upper,Min_t,Tol_t,bs,out(out_i),R2(r),-1);
                                
                                if Tol ~=1
                                    Status(S_met_idx) = 3;
                                end
                                
                                lambda_jj(end+1) = x(R2(r))/x(out(out_i)); %v_g(l)/v_g(s)
                                % lambda_jj(end+1) = model.lb(R2(r))/model.lb(out(out_i)); %v_g(l)/v_g(s)
                            else
                                lambda_jj(end+1) = V(R2(r))/V(out(out_i));
                            end
                            
                            k_Rs_temp = mean(kcat{strcmp(kcat{:,1},model.rxns(Rs(s))),4});
                            k_j_min(end+1) = min(k_Rs_temp./mean(kcat{strcmp(kcat{:,1},model.rxns(R2(r))),4}));
                            k_j_max(end+1) = max(k_Rs_temp./mean(kcat{strcmp(kcat{:,1},model.rxns(R2(r))),4}));
                            if isnan(k_j_min(end))
                                if strcmp(method,'equal')
                                    k_j_min(end)=1;
                                    k_j_max(end)=1;
                                elseif strcmp(method,'average')
                                    k_j_min(end)=AV;
                                    k_j_max(end)=AV;
                                elseif strcmp(method,'median')
                                    k_j_min(end)=ME;
                                    k_j_max(end)=ME;
                                end
                            end
                            
                            s_jj(end+1) = abs(N(M(P),Rs(s)));
                        end
                    end
                    Sum_j_min(out_i) = sum(s_jj.*lambda_jj.*k_j_min);
                    Sum_j_max(out_i) = sum(s_jj.*lambda_jj.*k_j_max);
                end
                
                B_temp_min = Sum_i./Sum_j_max;
                B_temp_max = Sum_i./Sum_j_min;
                ReactionSet_temp(end+1:end+size([in out],1),1:2) = [in out];
            end
        end
        if SCC(S_met_idx)==1
            B_min{S_met_idx}(end+1:end+length(B_temp_min)) = B_temp_min;
            B_max{S_met_idx}(end+1:end+length(B_temp_max)) = B_temp_max;
            ReactionSet{S_met_idx}(end+1:end+size(ReactionSet_temp,1),:) = ReactionSet_temp;
            ODE{S_met_idx}(end+1:end+length(B_temp_min)) = M(P);
            B_temp_min = [];B_temp_max = [];
            ReactionSet_temp = [];
        end
    end
end
end
