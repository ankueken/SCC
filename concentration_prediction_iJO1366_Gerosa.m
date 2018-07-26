changeCobraSolver('tomlab_cplex')

%% read experimental data
Data_c = readtable('Gerosa_data_E_coli_metabolites.csv','Delimiter',',','HeaderLines',1);
SR = readtable('Gerosa_data_E_coli_specific_rates.csv','Delimiter',',','HeaderLines',1);

%% load results of cytosolic SCC metabolites
SCC_naive = load('E_coli_iJO1366_cytosolic_SCC.mat');
SCC_i = Data_c.model_ID;  % index of SCC metabolites comprised in experimental dataset
SCC_if = SCC_naive.SCC_if; % index of all model SCC metabolites
ReactionSet = SCC_naive.ReactionSet; ODE = SCC_naive.ODE;
Model_o = load('Results/Result_E_coli_K12_iJO1366.mat','Model_pro','SCC');
relevant_rxns = unique(cell2mat(ReactionSet));
epsilon=0;epsObj=1e-8;

%% add metFormulas to model
iJO = readCbModel('Models/E_coli_K12_iJO1366.xml');
iJO = removeMetabolites(iJO,setdiff(iJO.mets,Model_o.Model_pro.mets));
Model_o.Model_pro.metFormulas = iJO.metFormulas;

%% specify objective
Model_o.Model_pro.ub(:) = 1000;
Model_o.Model_pro.lb(:) = 1e-10;
Model_o.Model_pro.csense(1:length(Model_o.Model_pro.mets),1) = 'E';
% Model_o.Model_pro.c(:) = 0;
% Model_o.Model_pro.c(1110) = 1; % max ATP synthesis alone
Model_o.Model_pro.c(:) = -1000;
Model_o.Model_pro.c(1110) = 0.001;

%% get sugar import
Sugar=find(cellfun(@isempty,regexp(Model_o.Model_pro.metFormulas,'C[3-7]H[1-9]*O[3-7]*'))==0);
Import_rxns = find(all(Model_o.Model_pro.S>=0));

[~,Import_Sugar] = find(Model_o.Model_pro.S(Sugar,Import_rxns)~=0);
Import_Sugar=Import_rxns(Import_Sugar);

%% learn B from acetate condition

Model_pro = Model_o.Model_pro;
Acetate=find(strcmp(Model_pro.mets,'ac[e]'));
Import_ac = intersect(Import_rxns,find(Model_pro.S(Acetate,:)>0));

% constrain biomass to be the measured growth rate
Model_pro.lb(7) = SR.Acetate(1);Model_pro.ub(7) = SR.Acetate(1);

Export = find(cellfun(@isempty,strfind(Model_pro.rxns,'EX_'))==0);
Model_pro.lb(setdiff(Export,relevant_rxns)) = 0;

% block uptake of other sugars
Model_pro.lb(setdiff(Import_Sugar,Import_ac))=0;
Model_pro.ub(setdiff(Import_Sugar,Import_ac))=0;
Model_pro.ub(Import_ac) = SR.Acetate(2)+SR.Acetate_1(2);
Model_pro.lb(Import_ac) = SR.Acetate(2)-SR.Acetate_1(2);

Sol = optimizeCbModel(Model_pro);
f = Sol.f;
Vc = Sol.x;

T=evalc('[vmin,vmax]=fluxVariability(Model_pro,100,''max'',Model_pro.rxns(unique(cell2mat(ReactionSet))))');
bs = ones(size(Model_pro.rxns));
bs(find(vmin>1)) = -1;
vp_vs_ratio_B=[];

% get min/max ratio of relevant reaction pairs (fractional LP)
for i=1:length(SCC_i)
    for j=1:size(ReactionSet{SCC_i(i)},1)
        Model_pro_FVA = Model_pro;
        
        y_i = zeros(1,size(Model_pro.S,2));
        y_i(ReactionSet{SCC_i(i)}(j,2)) = 1;
        
        Model_pro_FVA.S = [Model_pro.S zeros(size(Model_pro.S,1),1);  % steady state
            Model_o.Model_pro.c' 0; % at optimal solution
            y_i bs(find(y_i==1))  % y_j + bt = 1, b=0
            ];
        
        Model_pro_FVA.csense = [Model_pro.csense; 'G'; 'E'];
        
        Model_pro_FVA.lb = [Model_pro.lb; 0];
        Model_pro_FVA.ub = [Model_pro.ub; 10000];
        Model_pro_FVA.b = [Model_pro.b; f+(f*epsObj); 1];
        Model_pro_FVA.c = zeros(size(Model_pro.S,2)+1,1);
        Model_pro_FVA.c(ReactionSet{SCC_i(i)}(j,1)) = 1;
        
        Sol_r = optimizeCbModel(Model_pro_FVA,'min');
        y=Sol_r.x;
        if Sol_r.stat~=1
            [y, ~, ~,~, ~, ~, ~, Stat] = cplex(Model_pro_FVA.c, Model_pro_FVA.S, Model_pro_FVA.lb, Model_pro_FVA.ub, Model_pro_FVA.b, Model_pro_FVA.b);
            if Stat~=1
                vp_vs_ratio_B{i}(j,1) = round((Vc(ReactionSet{SCC_i(i)}(j,1))/Vc(ReactionSet{SCC_i(i)}(j,2)))*1e9)/1e9;
            else
                x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e10)/1e10;
                vp_vs_ratio_B{i}(j,1) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
            end
        else
            
            x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e10)/1e10;
            vp_vs_ratio_B{i}(j,1) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
        end
        
        Sol_r = optimizeCbModel(Model_pro_FVA,'max');
        y=Sol_r.x;
        if Sol_r.stat~=1
            [y, ~, ~,~, ~, ~, ~, Stat] = cplex(-Model_pro_FVA.c, Model_pro_FVA.S, Model_pro_FVA.lb, Model_pro_FVA.ub, Model_pro_FVA.b, Model_pro_FVA.b);
            
            if Stat~=1
                vp_vs_ratio_B{i}(j,2) = round((Vc(ReactionSet{SCC_i(i)}(j,1))/Vc(ReactionSet{SCC_i(i)}(j,2)))*1e9)/1e9;
            else
                x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e10)/1e10;
                vp_vs_ratio_B{i}(j,2) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
            end
        else
            x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e10)/1e10;
            vp_vs_ratio_B{i}(j,2) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
        end
    end
end

B_c1_max=[];B_c1_min=[];
for i=1:length(SCC_i)
    u = unique(ODE{SCC_i(i)});
    o1=[];u1=[];
    for j=1:length(u)
        temp_min = repmat(Data_c{i,3},size(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),:),1),1)./(vp_vs_ratio_B{i}((ODE{SCC_i(i)}==u(j)),2));
        temp_max = repmat(Data_c{i,3},size(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),:),1),1)./(vp_vs_ratio_B{i}((ODE{SCC_i(i)}==u(j)),1));
        B_c1_min{i}(find(ODE{SCC_i(i)}==u(j)),:) = temp_min;B_c1_max{i}(find(ODE{SCC_i(i)}==u(j)),:) = temp_max;
    end
    B_c1_min{i}(B_c1_min{i}==Inf)=NaN;B_c1_max{i}(B_c1_max{i}==Inf)=NaN;
    B_c1_min{i}(isnan(B_c1_max{i}))=NaN;B_c1_max{i}(isnan(B_c1_min{i}))=NaN;
end

%% concentration prediction for other carbon sources using learned constants
GR = SR{1,3:9};
imported_sugar_idx=[1005 1018 1030 1029 1038 1135 1148];
x_min_B1 = [];x_max_B1 = [];

for sample = 1:7
    Model_pro = Model_o.Model_pro;
    
    Model_pro.lb(7) = GR(sample);Model_pro.ub(7) = GR(sample);
    
    Import_met=imported_sugar_idx(sample);
    Import_met_rxns = intersect(Import_rxns,find(Model_pro.S(Import_met,:)>0));
    
    Export = find(cellfun(@isempty,strfind(Model_pro.rxns,'EX_'))==0);
    Model_pro.lb(setdiff(Export,relevant_rxns)) = 0;
    
    % block uptake of other sugars
    Model_pro.lb(setdiff(Import_Sugar,Import_met_rxns))=0;
    Model_pro.ub(setdiff(Import_Sugar,Import_met_rxns))=0;
    Model_pro.ub(Import_met_rxns) = SR{2,2+sample}+SR{2,11+sample};
    Model_pro.lb(Import_met_rxns) = SR{2,2+sample}-SR{2,11+sample};
    
    % constrain exchange rates by measured values
    if ~isnan(SR{3,2+sample})
        Model_pro.lb(intersect(Export,find(Model_pro.S(Acetate,:)<0))) = SR{3,2+sample}-(SR{3,11+sample});
        Model_pro.ub(intersect(Export,find(Model_pro.S(Acetate,:)<0))) = SR{3,2+sample}+(SR{3,11+sample});
    end
    if ~isnan(SR{4,2+sample})
        Model_pro.lb(intersect(Export,find(Model_pro.S(1070,:)<0))) = SR{4,2+sample}-(SR{4,11+sample});
        Model_pro.ub(intersect(Export,find(Model_pro.S(1070,:)<0))) = SR{4,2+sample}+(SR{4,11+sample});
    end
    if ~isnan(SR{5,2+sample})
        Model_pro.lb(intersect(Export,find(Model_pro.S(1009,:)<0))) = SR{5,2+sample}-(SR{5,11+sample});
        Model_pro.ub(intersect(Export,find(Model_pro.S(1009,:)<0))) = SR{5,2+sample}+(SR{5,11+sample});
    end
    
    Sol = optimizeCbModel(Model_pro);
    
    if ~isempty(Sol.x) % if model has solution upon restriction to carbon source
        Vg = Sol.x;
        f = Sol.f;
        
        % check for fractional LP if min flux of relevant reactions above 1
        T=evalc('[vmin,vmax]=fluxVariability(Model_pro,100,''max'',Model_pro.rxns(unique(cell2mat(ReactionSet))))');
        bs = ones(size(Model_pro.rxns));
        bs(find(vmin>1)) = -1;
        vp_vs_ratio=[];
        
        % min/max ratio of relevant reaction pairs
        for i=1:length(SCC_i)
            for j=1:size(ReactionSet{SCC_i(i)},1)
                Model_pro_FVA = Model_pro;
                
                y_i = zeros(1,size(Model_pro.S,2));
                y_i(ReactionSet{SCC_i(i)}(j,2)) = 1;
                
                Model_pro_FVA.S = [Model_pro.S zeros(size(Model_pro.S,1),1);  % steady state
                    Model_o.Model_pro.c' 0; % at optimal solution
                    y_i bs(find(y_i==1))  % y_j + bt = 1, b=0
                    ];
                
                Model_pro_FVA.csense = [Model_pro.csense; 'G'; 'E'];
                
                Model_pro_FVA.lb = [Model_pro.lb; 1e-10];
                Model_pro_FVA.ub = [Model_pro.ub; 10000];
                Model_pro_FVA.b = [Model_pro.b; f+(f*epsObj); 1];
                Model_pro_FVA.c = zeros(size(Model_pro.S,2)+1,1);
                Model_pro_FVA.c(ReactionSet{SCC_i(i)}(j,1)) = 1;
                
                Sol_r = optimizeCbModel(Model_pro_FVA,'min');
                y=Sol_r.x;
                if Sol_r.stat~=1
                    [y, ~, ~,~, ~, ~, ~, Stat] = cplex(Model_pro_FVA.c, Model_pro_FVA.S, Model_pro_FVA.lb, Model_pro_FVA.ub, Model_pro_FVA.b, Model_pro_FVA.b);
                    if Stat~=1
                        vp_vs_ratio{i}(j,1) = round((Vg(ReactionSet{SCC_i(i)}(j,1))/Vg(ReactionSet{SCC_i(i)}(j,2)))*1e10)/1e10;
                    else
                        x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e10)/1e10;
                        vp_vs_ratio{i}(j,1) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
                    end
                else
                    
                    x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e10)/1e10;
                    vp_vs_ratio{i}(j,1) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
                end
                
                Sol_r = optimizeCbModel(Model_pro_FVA,'max');
                y=Sol_r.x;
                if Sol_r.stat~=1
                    [y, ~, ~,~, ~, ~, ~, Stat] = cplex(-Model_pro_FVA.c, Model_pro_FVA.S, Model_pro_FVA.lb, Model_pro_FVA.ub, Model_pro_FVA.b, Model_pro_FVA.b);
                    
                    if Stat~=1
                        vp_vs_ratio{i}(j,2) = round((Vg(ReactionSet{SCC_i(i)}(j,1))/Vg(ReactionSet{SCC_i(i)}(j,2)))*1e10)/1e10;
                    else
                        x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e10)/1e10;
                        vp_vs_ratio{i}(j,2) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
                    end
                else
                    x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e10)/1e10;
                    vp_vs_ratio{i}(j,2) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
                end
            end
        end
        
        
        % calculate concentration for SCC metabolites using the min/max
        % ratio of the relevant reaction pairs and the constants learned
        % from acetate as carbon source
        B_c_min=B_c1_min;B_c_max=B_c1_max;
        
        for i=1:length(SCC_i)
            vp_vs_ratio{i}(vp_vs_ratio{i}==Inf)=NaN;
            if ~all(isnan(B_c_min{i})) && ~all(isnan(vp_vs_ratio{i}(:)))
                u = unique(ODE{SCC_i(i)});
                o1=[];u1=[];
                for j=1:length(u)
                    if  ~all(isnan(B_c_min{i}(ODE{SCC_i(i)}==u(j)))) && ~all(isnan(vp_vs_ratio{i}(ODE{SCC_i(i)}==u(j),1))) && ~all(isnan(vp_vs_ratio{i}(ODE{SCC_i(i)}==u(j),2)))
                        o1(j) = min(B_c_min{i}(ODE{SCC_i(i)}==u(j)) .* vp_vs_ratio{i}(ODE{SCC_i(i)}==u(j),1))';
                        u1(j) = max(B_c_max{i}(ODE{SCC_i(i)}==u(j)) .* vp_vs_ratio{i}(ODE{SCC_i(i)}==u(j),2))';
                    end
                end
                x_min_B1(i,sample) = max(o1); x_max_B1(i,sample) = min(u1);
            else
                x_min_B1(i,sample) = NaN; x_max_B1(i,sample) = NaN;
            end
        end
    else
        x_min_B1(1:length(SCC_i),sample) = Inf; x_max_B1(1:length(SCC_i),sample) = Inf;
    end
    
    x_min_B1_temp = x_min_B1;
    x_min_B1((round((x_min_B1./x_max_B1)*1e3)/1e3)>1) = NaN;
    x_max_B1((round((x_min_B1_temp./x_max_B1)*1e3)/1e3)>1) = NaN;
end


    
   
