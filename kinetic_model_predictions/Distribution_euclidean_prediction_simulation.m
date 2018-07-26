%% Euclidean distance between predicted metabolite concentration and metabolite concentration simulated with kinetic model
% - We use steady-state metabolite concentrations and flux distributions simulated with an kinetic E. coli model using
%   100 different initial conditions
% - We take each flux distribution simulated and rate constants provided in the kinetic model to predict a concentration
%   for the SCC metabolites
% - The simulated concentration value is compared to the predicted value by Euclidean distance

% -------------------------------------------------------------------------

% load kinetic model
addpath('Marans_Ecoli_kinetic_core_model/Users/')
load('Marans_Ecoli_kinetic_core_model/Users/Data.mat')

% put it to Cobra format
Model.S = Network_Data.S_f_b;
Model.lb = ones(size(Network_Data.S_f_b,2),1)*1e-6;
Model.ub = ones(size(Network_Data.S_f_b,2),1)*1000;
Model.rev = zeros(size(Network_Data.S_f_b,2),1);
Model.rxns = Network_Data.rxn_f_b;
Model.mets = Network_Data.metab;
Model.c = zeros(size(Network_Data.rxn_f_b));
Model.b = zeros(size(Network_Data.metab));

% create table of kcat values used by function get_SCC_kcat
% 1. column reaction names
% 4. column kcat values
kcat=table(Model.rxns,cell(size(Model.rxns)),cell(size(Model.rxns)),Elementary_Kinetic_Param);

% run SCC check
cd ../
    [SCC, Status, B_min, B_max, ReactionSet, ODE] = get_SCC_kcat(Model,[],[],kcat,'equal',1:length(Model.mets));
cd kinetic_model_predictions/

% load simulated data obtained with function perturbationAnalysis using 1, 5, 10 and 20% as input
R1=load('Perturb_all_1pc_1.mat');
R2=load('Perturb_all_5pc_1.mat');
R3=load('Perturb_all_10pc_1.mat');
R4=load('Perturb_all_20pc_1.mat');

Initial = [R1.START R2.START R3.START R4.START]; % initial concentrations used
Ref_flux_x = [R1.Ref_flux R2.Ref_flux R3.Ref_flux R4.Ref_flux]; % simulated steady-state flux
Ref_conc_x = [R1.Ref_conc R2.Ref_conc R3.Ref_conc R4.Ref_conc]; % simulated steady-state concentration
Ref_conc_x = Ref_conc_x(SCC==1,:);

SCC_i = find(SCC==1); % index of SCC metabolites

for ref=1:size(Ref_flux_x,2) % for each flux distribution simulated
    
    V = Ref_flux_x(:,ref);
    
    % calculate concentration of each SCC metabolite
    for i=1:length(SCC_i)
        if all(~isnan(B_min{i}))
            u = unique(ODE{SCC_i(i)});
            o1=[];u1=[];
            for j=1:length(u)
                o1(j) = min(B_min{SCC_i(i)}(ODE{SCC_i(i)}==u(j)) .* (V(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),1))./V(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),2))))';
                u1(j) = max(B_max{SCC_i(i)}(ODE{SCC_i(i)}==u(j)) .* (V(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),1))./V(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),2))))';
            end
                x_min(i,1) = max(o1); x_max(i,1) = min(u1);
        else
            x_min(i,1) = NaN; x_max(i,1) = NaN;
        end
    end
    
    x_min_temp = x_min;
    x_min((round((x_min./x_max)*1e3)/1e3)>1) = NaN;
    x_max((round((x_min_temp./x_max)*1e3)/1e3)>1) = NaN;
        
    % euclidean distance
    euc_av(ref)=pdist([x_min(~isnan(x_min))'; Ref_conc_x(~isnan(x_min),ref)'])/length(x_min(~isnan(x_min)));
  
end

% plot distribution Euclidean distance between predicted and simulated concentration of SCC
% metabolites
[H,b] = hist(euc_av,[0 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5]);
bar(H/sum(H),0.6)
set(gca,'XTickLabel',[0 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5])
colormap([0.5,0.5,0.5])
ylabel('Fraction')
xlabel({'Average Euclidean distance between simulated'; 'and predicted concentration'})