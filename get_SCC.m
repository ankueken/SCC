function SCC = get_SCC(model,fctable,blk)
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
% Input: Model in Cobra format (unblocked, irreversible reactions)
%        fctable: coupling matrix obtained from F2C2 (optional)
%        blk: vector indicating blockd reactions obtained from F2C2 (optional)
%
% Output:
%           SCC: vector of length m x 1
%                0 if approach is not applicable, 1 if species shows SCC

w=[];
addpath(genpath('F2C2/'))
network=CobraToF2C2(model);

%% check input data
if ~all(model.rev==0)
    warning('the network contains reversible reactions')
end

% calculate coupling matrix
if isempty(fctable) || isempty(blk)
    T=evalc('[fctable,blk]=F2C2(''glpk'',network);');
end

if ~all(blk==0)
    warning('blocked reactions detected during coupling matrix calculation')
    w=lastwarn;
end

% define variables
N = network.stoichiometricMatrix;
num_mets = size(N,1);
num_rxns = size(N,2);
SCC=zeros(size(N,1),1);

%% start SCC detection

Nplus=N.*(N<0);

%   S_met_idx: index of metabolite S for which we check SCC
for S_met_idx = 1:num_mets
    
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
        Cp_S=unique(Cp_S','rows','stable')';
        
        % check that complex Cp_S is still in the network
        Cp_S_check=cell(1,size(Cp_S,2));
        
        for c = 1:size(Cp_S,2)
            % case 1: zero complex
            C1=Cp_S(:,c);
            if all(C1==0)
                Cp_S_check{c}=[find(all(N<=0)) find(all(N>=0))]; % find import export reactions
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
        Cx_S=unique(Cx_S','rows','stable')';
        
        % check that complex is still in the network
        Cx_S_check=cell(1,size(Cx_S,2));
        for c = 1:size(Cx_S,2)
            % case 1: zero complex
            C1=Cx_S(:,c);
            if all(C1==0)
                Cx_S_check{c}=[find(all(N<=0)) find(all(N>=0))];
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
            end
        end
    end
    
end
end
