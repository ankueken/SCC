% load SCC results for E. coli model obtained from frunction main()
load('Results/Result_E_coli_K12_iJO1366.mat','Model_pro','SCC')
kcat=readtable('kcat_list_all.csv'); % read list of kcat values

% cytosolic metabolites
cytosolic = find(cellfun(@isempty,strfind(Model_pro.mets,'[c]'))==0);

Model_v = Model_pro;
Model_v.ub(1:length(Model_v.rxns)) = 1e4;
Model_v.lb(1:length(Model_v.rxns)) = 1e-7;
Model_v.c(7) = 1;
Sol = optimizeCbModel(Model_v);

% find SCC from set of cytosolic metabolites
[SCC, Status, B_min, B_max, ReactionSet, ODE] = ...
     get_SCC_kcat(Model_v,[],[],kcat,'equal',cytosolic,Sol.x);

SCC_if = find(SCC==1);

save('E_coli_iJO1366_cytosolic_SCC.mat')