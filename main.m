function main(f)
%% get SCC metabolites for models in folder Models/

% Input: number of the model in folder Models/ to load [1-14]

% - first remove blocked reactions
% - split reversible reactions into irreversible 
% - get SCC metabolites 
% - save result into folder Results/

addpath(genpath('F2C2/'))

% load model
File=dir('Models/*_*.*');
if File(f).name(end)=='t'
    load(strcat('Models/',File(f).name))
    if exist('RECON1')
        model=RECON1;
        model.rules=cell(size(model.rxns));
    end
else
    model=readCbModel(strcat('Models/',File(f).name));
end

% choose one out of the multiple biomass reactions and remove the others
if strcmp(File(f).name,'E_coli_K12_iJO1366.xml')
    model = removeRxns(model,model.rxns(8));
elseif strcmp(File(f).name,'Arabidopsis_CoreModel.xml')
    model = removeRxns(model,model.rxns(546:548));
end

model = checkCobraModelUnique(model);
disp('Reaction(s) with lower and upper bound set to zero:')
disp(model.rxns(intersect(find(model.lb==0),find(model.ub==0))))
model = removeRxns(model,model.rxns(intersect(find(model.lb==0),find(model.ub==0))));

[Exit,Model_pro] = model_preprocessing(model);  % remove blocked rxns, split rxns

SCC = get_SCC(Model_pro,[],[]); % find SCC metabolites

% save result
save(strcat('Results/Result_',File(f).name(1:end-3),'mat'))
end
