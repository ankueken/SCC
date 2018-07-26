% load model
addpath('Marans_Ecoli_kinetic_core_model/Users/')
load('Marans_Ecoli_kinetic_core_model/Users/Data.mat')

% load simulated data obtained with function perturbationAnalysis using 1, 5, 10 and 20% as input
R1=load('Perturb_all_1pc_1.mat');
R2=load('Perturb_all_5pc_1.mat');
R3=load('Perturb_all_10pc_1.mat');
R4=load('Perturb_all_20pc_1.mat');

Initial = [R1.START R2.START R3.START R4.START]; % initial concentrations used
Ref_flux_x = [R1.Ref_flux R2.Ref_flux R3.Ref_flux R4.Ref_flux]; % simulated steady-state flux
Ref_conc_x = [R1.Ref_conc R2.Ref_conc R3.Ref_conc R4.Ref_conc]; % simulated steady-state concentration

% bring kinetic model to cobra format
Model.S = Network_Data.S_f_b;
Model.lb = ones(size(Network_Data.S_f_b,2),1)*1e-8;
Model.ub = ones(size(Network_Data.S_f_b,2),1)*1000;
Model.rev = zeros(size(Network_Data.S_f_b,2),1);
Model.rxns = Network_Data.rxn_f_b;
Model.mets = Network_Data.metab;
Model.b = zeros(size(Network_Data.S_f_b,1),1);
Model.c = -ones(size(Network_Data.S_f_b,2),1)*0.01; % minimize total flux
Export = find(cellfun(@isempty,strfind(Model.rxns,'EX_'))==0);
Model.lb(Export) = min(Ref_flux_x(Export,:)')'; % constrain import/export flux by minimum/maximum value obtained in simulations
Model.ub(Export) = max(Ref_flux_x(Export,:)')';
Model.c(123) = 1; % maximize ATP synthesis

Sol = optimizeCbModel(Model);
f = Sol.f;
Vc = Sol.x;

kcat=table(Model.rxns,cell(size(Model.rxns)),cell(size(Model.rxns)),Elementary_Kinetic_Param);

% run SCC check
cd ../
    [SCC, Status, B_min, B_max, ReactionSet,ODE] = get_SCC_kcat(Model,[],[],kcat,'equal',1:length(Model.mets));
cd kinetic_model_predictions/

SCC_i = find(SCC==1);

T=evalc('[vmin,vmax]=fluxVariability(Model,100,''max'',Model.rxns(unique(cell2mat(ReactionSet))))');
bs = ones(size(Model.rxns));
bs(find(vmin>1)) = -1;
vp_vs_ratio=[];

for i=1:length(SCC_i) % for each SCC
    for j=1:size(ReactionSet{SCC_i(i)},1) % min/max ratio for each relevant reaction pair v_p and v_s of this SCC
        Model_pro_FVA = Model;
        
        y_i = zeros(1,size(Model.S,2));
        y_i(ReactionSet{SCC_i(i)}(j,2)) = 1;
        
        Model_pro_FVA.S = [Model.S zeros(size(Model.S,1),1);  % steady state
            Model.c' 0; % at optimal solution
            y_i bs(find(y_i==1))  % y_j + bt = 1, b=0
            ];
        
        Model_pro_FVA.lb = [Model.lb; 1e-8];
        Model_pro_FVA.ub = [Model.ub; 10000];
        Model_pro_FVA.b = [Model.b; f; 1];
        Model_pro_FVA.c = zeros(size(Model.S,2)+1,1);
        Model_pro_FVA.c(ReactionSet{SCC_i(i)}(j,1)) = 1;
        
        Sol_r = optimizeCbModel(Model_pro_FVA,'min');
        y=Sol_r.x;
        if Sol_r.stat~=1
            [y, ~, ~,~, ~, ~, ~, Stat] = cplex(Model_pro_FVA.c, Model_pro_FVA.S, Model_pro_FVA.lb, Model_pro_FVA.ub, Model_pro_FVA.b, Model_pro_FVA.b);
            
            if Stat~=1
                vp_vs_ratio{i}(j,1) = round((Vc(ReactionSet{SCC_i(i)}(j,1),1)/Vc(ReactionSet{SCC_i(i)}(j,2),1))*1e9)/1e9;
            else
                x = round(((1/y(end)).*y(1:size(Model.S,2)))*1e9)/1e9;
                vp_vs_ratio{i}(j,1) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
            end
        else
            
            x = round(((1/y(end)).*y(1:size(Model.S,2)))*1e9)/1e9;
            vp_vs_ratio{i}(j,1) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
        end
        Sol_r = optimizeCbModel(Model_pro_FVA,'max');
        y=Sol_r.x;
        if Sol_r.stat~=1
            [y, ~, ~,~, ~, ~, ~, Stat] = cplex(-Model_pro_FVA.c, Model_pro_FVA.S, Model_pro_FVA.lb, Model_pro_FVA.ub, Model_pro_FVA.b, Model_pro_FVA.b);
            
            if Stat~=1
                vp_vs_ratio{i}(j,2) = round((Vc(ReactionSet{SCC_i(i)}(j,1),1)/Vc(ReactionSet{SCC_i(i)}(j,2),1))*1e9)/1e9;
            else
                x = round(((1/y(end)).*y(1:size(Model.S,2)))*1e9)/1e9;
                vp_vs_ratio{i}(j,2) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
            end
        else
            
            x = round(((1/y(end)).*y(1:size(Model.S,2)))*1e9)/1e9;
            vp_vs_ratio{i}(j,2) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
        end
    end
end

% range from kinetic simulations
Simulated_range = max(Ref_conc_x')'-min(Ref_conc_x')';

% range from prediction
V = Ref_flux_x;

SCC_i = find(SCC==1);

% get SCC concentration range using ratios calculated before
for i=1:length(SCC_i)
    if all(~isnan(B_min{i}))
        u = unique(ODE{SCC_i(i)});
        o1=[];u1=[];
        for j=1:length(u)
            o1(j) = min(B_min{SCC_i(i)}(ODE{SCC_i(i)}==u(j)) .* vp_vs_ratio{i}(ODE{SCC_i(i)}==u(j),1))';
            u1(j) = max(B_max{SCC_i(i)}(ODE{SCC_i(i)}==u(j)) .* vp_vs_ratio{i}(ODE{SCC_i(i)}==u(j),2))';
        end
        x_min(i,1) = max(o1); x_max(i,1) = min(u1);
    else
        x_min(i,1) = NaN; x_max(i,1) = NaN;
    end
end

x_min_temp = x_min;
x_min((round((x_min./x_max)*1e3)/1e3)>1) = NaN;
x_max((round((x_min_temp./x_max)*1e3)/1e3)>1) = NaN;

c_lb = min(Ref_conc_x')'; c_ub = max(Ref_conc_x')';

% plot comparison simulated vs. predicted concentrations
for i=1:length(SCC_i)
    hold on
    plot([i-0.2 i-0.2],[c_lb(SCC_i(i)) c_ub(SCC_i(i))],'r')
    hold on
    plot([i-0.3 i-0.1],[c_lb(SCC_i(i)) c_lb(SCC_i(i))],'r')
    hold on
    plot([i-0.3 i-0.1],[c_ub(SCC_i(i)) c_ub(SCC_i(i))],'r')
    hold on
    plot([i+0.2 i+0.2],[x_min(i,1) x_max(i,1)],'k')
    hold on
    plot([i+0.1 i+0.3],[x_min(i,1) x_min(i,1)],'k')
    hold on
    plot([i+0.1 i+0.3],[x_max(i,1) x_max(i,1)],'k')
end

set(gca,'XTick',1:length(SCC_i),'XTickLabel',Model.mets(SCC_i),'XTickLabelRotation',45,'YScale','log')
xlim([0 length(SCC_i)+1])
ylabel('mmol/gDW')

% correlations
[cval_lb,pval_lb]=corr(x_min,c_lb(SCC==1),'type','Pearson','rows','pairwise')
[cval_lb,pval_lb]=corr(x_max,c_ub(SCC==1),'type','Pearson','rows','pairwise')
[cval_lb,pval_lb]=corr(x_min,c_lb(SCC==1),'type','Spearman','rows','pairwise')
[cval_lb,pval_lb]=corr(x_max,c_ub(SCC==1),'type','Spearman','rows','pairwise')
[cval_lb,pval_lb]=corr(x_min,c_lb(SCC==1),'type','Kendall','rows','pairwise')
[cval_ub,pval_ub]=corr(x_max,c_ub(SCC==1),'type','Kendall','rows','pairwise')

[cval_range,pval_range]=corr(x_max-x_min,c_ub(SCC==1)-c_lb(SCC==1),'rows','pairwise')
[cval_rangeS,pval_rangeS]=corr(x_max-x_min,c_ub(SCC==1)-c_lb(SCC==1),'type','Spearman','rows','pairwise')
[cval_rangeS,pval_rangeS]=corr(x_max-x_min,c_ub(SCC==1)-c_lb(SCC==1),'type','Kendall','rows','pairwise')

% average euclidean distance
pdist([x_min(~isnan(x_min))'; c_lb(SCC_i(~isnan(x_min)))'])/length(SCC_i(~isnan(x_min)))
pdist([x_max(~isnan(x_max))'; c_ub(SCC_i(~isnan(x_max)))'])/length(SCC_i(~isnan(x_max)))

% euclidean distance log-transformed data
pdist([log(x_min(~isnan(x_min)))'; log(c_lb(SCC_i(~isnan(x_min))))'])/length(SCC_i(~isnan(x_min)))
pdist([log(x_max(~isnan(x_max)))'; log(c_ub(SCC_i(~isnan(x_max))))'])/length(SCC_i(~isnan(x_max)))

% relative euclidean distance
pdist([x_min(~isnan(x_min))'./max(x_min(~isnan(x_min))); c_lb(SCC_i(~isnan(x_min)))'./max(c_lb(SCC_i(~isnan(x_min))))])/length(SCC_i(~isnan(x_min)))
pdist([x_max(~isnan(x_max))'./max(x_max(~isnan(x_max))); c_ub(SCC_i(~isnan(x_max)))'./max(c_ub(SCC_i(~isnan(x_max))))])/length(SCC_i(~isnan(x_max)))

