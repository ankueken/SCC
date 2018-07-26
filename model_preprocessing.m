function [Exit, model_noblk_splitted] = model_preprocessing(model)
% Function used for model preprocessing
%
% Requires COBRA toolbox, glpk solver, F2C2 toolbox
%
% Input: Model in Cobra format
% 
% Output: 
%   Model without blocked reactions and splitted reversible reactions
%
%   Exit: 0
%   The resulting model is able to carry positive flux in all reactions and 
%   contains no blocked reactions, reversible reactions are splitted 
%
%   Exit: 1
%   The resulting model contains no blocked reactions 
%   and reversible reactions are splitted 
%
%   Exit: 2
%   New blocked reactions after reaction splitting

%% set bounds
model.ub=ones(size(model.ub))*1000;
model.lb=ones(size(model.lb))*0;
model.lb(model.rev==1)=-1000;
 
if ~isfield(model,'grRules') && isfield(model,'rules')
    model.grRules=model.rules;
elseif ~isfield(model,'rules') && isfield(model,'grRules')
    model.rules=model.grRules;
else
    model.rules=cell(size(model.rxns));
    model.grRules=cell(size(model.rxns));
end

net_check.stoichiometricMatrix = full(model.S);
net_check.Metabolites = model.mets;
net_check.Reactions = model.rxns;
net_check.reversibilityVector = model.rev;

% check and remove blocked reactions
[~, blocked] = F2C2('glpk',net_check);

if ~all(blocked==0)
    model_noblk=removeRxns(model,model.rxns(blocked==1),model.rev(blocked==1),'true');
    fprintf('remove %d blocked reactions \n',length(find(blocked==1)))
else
    model_noblk=model;
    disp('*********************************')
    disp('No blocked reactions.')
end
% get biomass reaction index
indBio=find(cellfun(@isempty,strfind(model_noblk.rxns,'Bio'))==0);
if isempty(indBio)
    indBio=find(cellfun(@isempty,strfind(model_noblk.rxns,'BIO'))==0);
end
if isempty(indBio)
    indBio=find(cellfun(@isempty,strfind(model_noblk.rxns,'bio'))==0);
end

% split reversible reactions into two irreversible reactions
disp('split reversible reactions ...')
model_noblk_splitted=split_rxns(model_noblk);

% get reaction coupling matrix using F2C2 toolbox
ACR_check_net.stoichiometricMatrix = full(model_noblk_splitted.S);
ACR_check_net.Metabolites = strcat(model_noblk_splitted.mets,' -- ',model_noblk_splitted.metNames);
ACR_check_net.Reactions = strcat(model_noblk_splitted.rxns,' -- ',model_noblk_splitted.rxnNames);
ACR_check_net.reversibilityVector = model_noblk_splitted.rev;

[fctable,blk]=F2C2('glpk',ACR_check_net);

% check again if there is no blocked reaction otherwise remove them
if ~all(blk==0)
    model_noblk_splitted=removeRxns(model_noblk_splitted,model_noblk_splitted.rxns(blk==1),model_noblk_splitted.rev(blk==1),'true');
    fprintf('remove %d blocked reactions \n',length(find(blocked==1)))
    ACR_check_net.stoichiometricMatrix = full(model_noblk_splitted.S);
    ACR_check_net.Metabolites = strcat(model_noblk_splitted.mets,' -- ',model_noblk_splitted.metNames);
    ACR_check_net.Reactions = strcat(model_noblk_splitted.rxns,' -- ',model_noblk_splitted.rxnNames);
    ACR_check_net.reversibilityVector = model_noblk_splitted.rev;

    [fctable,blk]=F2C2('glpk',ACR_check_net);
end

if sum(blk)==0

    % test for positive steady state flux distribution
    model_noblk_splitted.lb=ones(size(model_noblk_splitted.S,2),1)*1e-5;

    Sol = optimizeCbModel(model_noblk_splitted);

    if ~isempty(Sol.x)
        Exit = 0;
    else
        Exit = 1;
    end
else
    Exit = 2;
end
end


