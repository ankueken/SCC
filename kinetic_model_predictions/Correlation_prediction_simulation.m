% Code to compare (correlation, euclidean distance) concentration of SCC metabolites in an E. coli kinetic model 
% compare simulated steady-state concentration from 100 different initial
% conditions and concentrations predicted using the respective staedy-state flux
% distributions

%%
addpath('Marans_Ecoli_kinetic_core_model/Users/')
load('Marans_Ecoli_kinetic_core_model/Users/Data.mat')

% put Maranas kinetic model in Cobra format
Model.S = Network_Data.S_f_b;
Model.lb = ones(size(Network_Data.S_f_b,2),1)*1e-6;
Model.ub = ones(size(Network_Data.S_f_b,2),1)*1000;
Model.rev = zeros(size(Network_Data.S_f_b,2),1);
Model.rxns = Network_Data.rxn_f_b;
Model.mets = Network_Data.metab;
Model.c = zeros(size(Network_Data.rxn_f_b));
Model.b = zeros(size(Network_Data.metab));

% create table of kcat values
% 1. column reaction names
% 4. column kcat values
kcat=table(Model.rxns,cell(size(Model.rxns)),cell(size(Model.rxns)),Elementary_Kinetic_Param);

% run SCC check
cd ../
[SCC, Status, B_min, B_max, ReactionSet, ODE] = get_SCC_kcat(Model,[],[],kcat,'equal',1:length(Model.mets));
cd kinetic_model_predictions/

disp('Number of SCC:')
disp(length(find(SCC==1)))

% load simulated data obtained with function perturbationAnalysis using 1, 5, 10 and 20% as input
R1=load('Perturb_all_1pc_1.mat');
R2=load('Perturb_all_5pc_1.mat');
R3=load('Perturb_all_10pc_1.mat');
R4=load('Perturb_all_20pc_1.mat');

Initial = [R1.START R2.START R3.START R4.START]; % initial concentrations used
Ref_flux_x = [R1.Ref_flux R2.Ref_flux R3.Ref_flux R4.Ref_flux]; % simulated steady-state flux
Ref_conc_x = [R1.Ref_conc R2.Ref_conc R3.Ref_conc R4.Ref_conc]; % simulated steady-state concentration

% range from kinetic simulations
Simulated_range = max(Ref_conc_x')'-min(Ref_conc_x')';

% predict concentration range using all simulated flux distributions
V = Ref_flux_x; 
SCC_i = find(SCC==1);  % index of SCC metabolites

for i=1:length(SCC_i)
    if all(~isnan(B_min{i}))
        u = unique(ODE{SCC_i(i)});
        o1=[];u1=[];
        for j=1:length(u)
            o1(j) = min(B_min{SCC_i(i)}(ODE{SCC_i(i)}==u(j)) .* (V(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),1),:)./V(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),2),:)))';
            u1(j) = max(B_max{SCC_i(i)}(ODE{SCC_i(i)}==u(j)) .* (V(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),1),:)./V(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),2),:)))';
        end
        x_min(i,1) = max(o1); x_max(i,1) = min(u1);
    else
        x_min(i,1) = NaN; x_max(i,1) = NaN;
    end
end

x_min_temp = x_min;
x_min((round((x_min./x_max)*1e3)/1e3)>1) = NaN;
x_max((round((x_min_temp./x_max)*1e3)/1e3)>1) = NaN;

% maximum and minimum concentration simulated
c_lb = min(Ref_conc_x')'; c_ub = max(Ref_conc_x')';

% correlation simulated and predicted minimum and maximum concentration
[cval_lb,pval_lb]=corr(x_min,c_lb(SCC==1),'type','Pearson')
[cval_lb,pval_lb]=corr(x_max,c_ub(SCC==1),'type','Pearson')
[cval_lb,pval_lb]=corr(x_min,c_lb(SCC==1),'type','Spearman')
[cval_lb,pval_lb]=corr(x_max,c_ub(SCC==1),'type','Spearman')
[cval_lb,pval_lb]=corr(x_min,c_lb(SCC==1),'type','Kendall')
[cval_ub,pval_ub]=corr(x_max,c_ub(SCC==1),'type','Kendall')

% correlation simulated and predicted concentration range
[cval_range,pval_range]=corr(x_max-x_min,c_ub(SCC==1)-c_lb(SCC==1))
[cval_rangeS,pval_rangeS]=corr(x_max-x_min,c_ub(SCC==1)-c_lb(SCC==1),'type','Spearman')
[cval_rangeS,pval_rangeS]=corr(x_max-x_min,c_ub(SCC==1)-c_lb(SCC==1),'type','Kendall')

% average euclidean distance
pdist([x_min'; c_lb(SCC_i)'])/length(SCC_i)
pdist([x_max'; c_ub(SCC_i)'])/length(SCC_i)

% euclidean distance log-transformed data
pdist([log(x_min)'; log(c_lb(SCC_i))'])/length(SCC_i)
pdist([log(x_max)'; log(c_ub(SCC_i))'])/length(SCC_i)

% relative euclidean distance
pdist([x_min'./max(x_min); c_lb(SCC_i)'./max(c_lb(SCC_i))])/length(SCC_i)
pdist([x_max'./max(x_max); c_ub(SCC_i)'./max(c_ub(SCC_i))])/length(SCC_i)

