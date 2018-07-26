%% compare shadow price and simulated concentration range for model SCC metabolites

% load model
addpath('Marans_Ecoli_kinetic_core_model/Users/')
load('Marans_Ecoli_kinetic_core_model/Users/Data.mat')

% put it to cobra format
Model.S = Network_Data.S_f_b;
Model.lb = ones(size(Network_Data.S_f_b,2),1)*1e-5;
Model.ub = ones(size(Network_Data.S_f_b,2),1)*1000;
Model.rev = zeros(size(Network_Data.S_f_b,2),1);
Model.rxns = Network_Data.rxn_f_b;
Model.mets = Network_Data.metab;
Model.b = zeros(size(Network_Data.S_f_b,1),1);
Model.c = zeros(size(Network_Data.S_f_b,2),1);
Export = find(cellfun(@isempty,strfind(Model.rxns,'EX_'))==0);
Model.c(Export(find(cellfun(@isempty,strfind(Model.rxns(Export),'3_f'))==0))) = 1; % maximize
Model.c(Export(find(cellfun(@isempty,strfind(Model.rxns(Export),'3_b'))==0))) = -1; % minimize

kcat=table(Model.rxns,cell(size(Model.rxns)),cell(size(Model.rxns)),Elementary_Kinetic_Param);

% check SCC 
cd ../
[SCC, concentration_range, Status, B_min, B_max, ReactionSet] = get_SCC_kcat(Model,[],[],kcat,'average',1:length(Model.mets));
cd kinetic_model_predictions/

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
c_ub = max(Ref_conc_x')'; c_lb = min(Ref_conc_x')';

% shadow price
saRequest.obj.index = SCC_i;
[v, ~, shadow_price, ~, f_k, ~, ~, Inform] = ...
    cplex(-Model.c, Model.S, Model.lb, Model.ub, Model.b,Model.b, ...
    [], [], [], [], [], [], [], [], ...
    [], [], [], [], [], [], [], ...
    [], [], saRequest);

% compare results using correlation
[c,p]=corr(shadow_price(SCC_i),Simulated_range(SCC_i)) % range
[c,p]=corr(shadow_price(SCC_i),Simulated_range(SCC_i),'type','Spearman') % range Spearman
[c,p]=corr(shadow_price(SCC_i),log(Simulated_range(SCC_i))) % log range
[c,p]=corr(shadow_price(SCC_i),std(Ref_conc_x(SCC_i,:)')'./mean(Ref_conc_x(SCC_i,:)')') % CV
[c,p]=corr(shadow_price(SCC_i),log(std(Ref_conc_x(SCC_i,:)')'./mean(Ref_conc_x(SCC_i,:)')')) % log CV



