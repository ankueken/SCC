%% Concentration changes upon reaction knock-out
changeCobraSolver('tomlab_cplex','QP')
changeCobraSolver('tomlab_cplex','LP')

% load model
addpath('Marans_Ecoli_kinetic_core_model/Users/')
load('Marans_Ecoli_kinetic_core_model/Users/Data.mat')

%% WT reference
% we use initial conditions provided by the model to simulate a steady-state
% flux distribution
[~,C,V]=solve_ode(Network_Data,0:100,Metabolite_and_Enzyme_Conc,Elementary_Kinetic_Param);

% we calculate SCC concentrations for that flux distribution
Model.S = Network_Data.S_f_b;
Model.lb = V(:,end);
Model.ub = V(:,end);
Model.rev = zeros(size(Network_Data.S_f_b,2),1);
Model.rxns = Network_Data.rxn_f_b;
Model.mets = Network_Data.metab;
Model.c = zeros(size(Network_Data.rxn_f_b));
Model.b = zeros(size(Network_Data.metab));

kcat=table(Network_Data.rxn_f_b,cell(size(Network_Data.rxn_f_b)),cell(size(Network_Data.rxn_f_b)),Elementary_Kinetic_Param);
cd ..
[SCC, Status, B_min, B_max, ReactionSet, ODE] = get_SCC_kcat(Model,[],[],kcat,'equal',1:length(Network_Data.metab),V);
cd kinetic_model_predictions/

SCC_i=find(SCC);

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
x_min((round((x_min./x_max)*1e2)/1e2)>1) = NaN;
x_max((round((x_min_temp./x_max)*1e2)/1e2)>1) = NaN;

table(Model.mets(SCC_i),x_min,x_max)

%% MOMA knock-out mutants
% knock-out reactions one by one
% record change in concentration for the SCC metabolites
% use MOMA to find knock-out flux distribution

for ko_idx=1:length(Model.rxns)
    disp(ko_idx)
    % we keep the input flux constant
    Export = find(all(Model.S<=0));
    Model_ko.S = Network_Data.S_f_b;
    Model_ko.lb = ones(size(Network_Data.S_f_b,2),1)*1e-10;
    Model_ko.ub = ones(size(Network_Data.S_f_b,2),1)*1000;
    Model_ko.rev = zeros(size(Network_Data.S_f_b,2),1);
    Model_ko.rxns = Network_Data.rxn_f_b;
    Model_ko.mets = Network_Data.metab;
    Model_ko.c = zeros(size(Network_Data.rxn_f_b));
    Model_ko.b = zeros(size(Network_Data.metab));
    Model_ko.lb(Export) = V(Export,end);
    Model_ko.ub(Export) = V(Export,end);
    Model_ko.lb(ko_idx) = 0;
    Model_ko.ub(ko_idx) = 0;
    
    toms 1474x1 x
    Vdel=V(:,end);
    constr = {Model_ko.S*x == Model_ko.b
              Model_ko.lb <= x <= Model_ko.ub};
    obj = x-Vdel;    
    options.norm='L2';
    options.prilev=0;
    [solutionDel,R]=ezsolve(obj,constr,[],options);
    
    s_del(:,ko_idx)=solutionDel.x(1:1474);
    
    SCC_i=find(SCC);
    V_del=solutionDel.x;
    Opt_diff(ko_idx)=R.f_k;
    Stat(ko_idx)=R.Inform;

    for i=1:length(SCC_i)
        if all(~isnan(B_min{i}))
            u = unique(ODE{SCC_i(i)});
            o1=[];u1=[];
            for j=1:length(u)
                o1(j) = min(B_min{SCC_i(i)}(ODE{SCC_i(i)}==u(j)) .* (V_del(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),1),:)./V_del(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),2),:)))';
                u1(j) = max(B_max{SCC_i(i)}(ODE{SCC_i(i)}==u(j)) .* (V_del(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),1),:)./V_del(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),2),:)))';
            end
            x_min_del(i,ko_idx) = max(o1); x_max_del(i,ko_idx) = min(u1);
        else
            x_min_del(i,ko_idx) = NaN; x_max_del(i,ko_idx) = NaN;
        end
    end

    x_min_temp = x_min_del;
    x_min_del((round((x_min_del./x_max_del)*1e2)/1e2)>1) = NaN;
    x_max_del((round((x_min_temp./x_max_del)*1e2)/1e2)>1) = NaN;

end

x_min_del(isnan(x_min_del))=0; % if no solution concentration is set to zero
change_SCC=repmat(x_min(~isnan(x_min)),1,size(x_min_del,2))-x_min_del(~isnan(x_min),:);

no_change=sum(abs(change_SCC')<1e-2)'; % no change in concentration 
dec_change=sum(change_SCC'>=1e-2)'; % decrease in concentration upon KO
inc_change=sum(change_SCC'<=-1e-2)'; % increase in concentration upon KO
table(Model.mets(SCC_i(~isnan(x_min))),no_change,dec_change,inc_change)

% fold-changes between WT and MOMA mutant (predictied fold changes)
PredFoldChange = (x_min_del-repmat(x_min,1,size(x_min_del,2)))./repmat(x_min,1,size(x_min_del,2));

save('Results_SCC_conc_upon_rxn_ko_fixed_EX.mat')
