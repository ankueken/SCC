function Model_MM = MM_split_rxns(Model)
%% function to split reactions into Michealis-Menten like form
% Input: cobra model
% Output: model including elementary reaction steps (Michealis-Menten like mass-acion) of Substrate + Enzyme = Complex
% and Complex -> Product

%% check completeness and set variables

if ~isfield(Model, 'rxns')
    disp('Error: Field "rxns" is missing.')
elseif ~isfield(Model, 'S')
    disp('Error: Field "S" is missing.')
elseif ~isfield(Model, 'rev')
    disp('Error: Field "rev" is missing.')
elseif ~isfield(Model, 'lb')
    disp('Error: Field "lb" is missing.')
elseif ~isfield(Model, 'ub')
    disp('Error: Field "ub" is missing.')
elseif ~isfield(Model, 'c')
    disp('Error: Field "c" is missing.')
end

% splitted Model
Model_MM.S = zeros(size(Model.S,1),0);
Model_MM.rxns = cell(0);
Model_MM.mets = Model.mets;
Model_MM.lb = zeros(0);
Model_MM.ub = zeros(0);
Model_MM.rev = zeros(0);
Model_MM.c = zeros(0);

if isfield(Model, 'rxnNames')
    Model_MM.rxnNames = cell(0);
end
if isfield(Model, 'rxnECNumbers')
    Model_MM.rxnECNumbers = cell(0);
end
if isfield(Model, 'rules')
    Model_MM.rules = cell(0);
end
if isfield(Model, 'rxnGeneMat')
    Model_MM.rxnGeneMat = zeros(0,size(Model.rxnGeneMat,2));
end
if isfield(Model, 'grRules')
    Model_MM.grRules = cell(0);
end
if isfield(Model, 'subSystems')
    Model_MM.subSystems = cell(0);
end
if isfield(Model, 'rxnReferences')
    Model_MM.rxnReferences = cell(0);
end
if isfield(Model, 'rxnNotes')
    Model_MM.rxnNotes = cell(0);
end
if isfield(Model, 'confidenceScore')
    Model_MM.confidenceScore = zeros(0);
end
if isfield(Model, 'proteins')
    Model_MM.proteins = cell(0);
end

if isfield(Model, 'metNames')
    Model_MM.metNames = Model.metNames;
end

S_positive = Model.S;
S_positive(S_positive<0) = 0;
S_negative = Model.S;
S_negative(S_negative>0) = 0;

% set list of enzymes
% fill empty fiels...
if ~isfield(Model, 'grRules')
    Model.grRules = cell(size(Model.rxns)); % every rxn gets one enzyme
    Model_MM.grRules = cell(0);
end    
    
E1=find(cellfun(@isempty,Model.grRules)==1);
E2=find(strcmp(Model.grRules,'-'));
for i=1:length(E1)
    Model.grRules{E1(i)} = strcat('Empty_grRule_',num2str(i));
end
for i=1:length(E2)
    Model.grRules{E2(i)} = strcat('Empty_grRule_',num2str(length(E1)+i));
end

[E,~,IE] = unique(Model.grRules); % Model.grRules = E(IE)

for i=1:length(E)
    Model_MM.S(end+1,:) = 0;
    Model_MM.mets{end+1,1} = strcat('Enzyme_GPR_', num2str(i)); 
    if isfield(Model, 'metNames')
        Model_MM.metNames{end+1,1} = strcat('Enzyme_GPR_', num2str(i));
    end
end


% split into MM-format:
for i=1:length(Model.rxns)
    
    Model_MM.S(1:size(S_negative,1),end+1) = S_negative(:,i); % A + E = ES
    Model_MM.S(size(S_negative,1)+IE(i),end) = -1; % add -E and +ES
    Model_MM.S(end+1,end) = 1; % add row for complex
    
    Model_MM.mets{end+1,1} = strcat(Model_MM.mets{size(S_negative,1)+IE(i)},'_complex');
    
    if isfield(Model, 'metNames')
        Model_MM.metNames{end+1,1} = strcat(Model_MM.metNames{size(S_negative,1)+IE(i)},'_complex');
    end
    
    Model_MM.rev(end+1,1)=  1;
    Model_MM.c(end+1,1)  =  0;
    Model_MM.lb(end+1,1) = -1000;
    Model_MM.ub(end+1,1) =  1000;
    
    Model_MM.S(1:size(S_positive,1),end+1) = S_positive(:,i); % ES -> B + E
    Model_MM.S([size(S_positive,1)+IE(i) size(Model_MM.S,1)],end) = [1; -1]; % add -ES and +E
    
    if Model.rev(i)==1
        Model_MM.rev(end+1,1)=  1;
        Model_MM.lb(end+1,1) = -1000;
        Model_MM.ub(end+1,1) =  1000;
    elseif Model.rev(i)==0
        Model_MM.rev(end+1,1)=  0;
        Model_MM.lb(end+1,1) =  0;
        Model_MM.ub(end+1,1) =  1000;
    else
        disp('Error in setting reversibility.')
    end
    
    if Model.c(i) ~=0
        Model_MM.c(end+1,1) = Model.c(i);
    else
        Model_MM.c(end+1,1)  =  0;
    end
    
    % rxn names
    Model_MM.rxns{end+1,1} = strcat(Model.rxns{i}, '_k1');
    Model_MM.rxns{end+1,1} = strcat(Model.rxns{i}, '_k2');
    
    if isfield(Model, 'rxnNames')
        Model_MM.rxnNames{end+1,1} = strcat(Model.rxnNames{i}, '_k1');
        Model_MM.rxnNames{end+1,1} = strcat(Model.rxnNames{i}, '_k2');
    end
    
    if isfield(Model, 'rxnECNumbers')
        Model_MM.rxnECNumbers{end+1,1} = Model.rxnECNumbers{i};
        Model_MM.rxnECNumbers{end+1,1} = Model.rxnECNumbers{i};
    end
    if isfield(Model, 'rules')
        Model_MM.rules{end+1,1} = Model.rules{i};
        Model_MM.rules{end+1,1} = Model.rules{i};
    end
    if isfield(Model, 'rxnGeneMat')
        Model_MM.rxnGeneMat(end+1,:) = Model.rxnGeneMat(i,:);
        Model_MM.rxnGeneMat(end+1,:) = Model.rxnGeneMat(i,:);
    end
    if isfield(Model, 'grRules')
        Model_MM.grRules{end+1,1} = Model.grRules{i};
        Model_MM.grRules{end+1,1} = Model.grRules{i};
    end
    if isfield(Model, 'subSystems')
        Model_MM.subSystems{end+1,1} = Model.subSystems{i};
        Model_MM.subSystems{end+1,1} = Model.subSystems{i};
    end
    if isfield(Model, 'rxnReferences')
        Model_MM.rxnReferences{end+1,1} = Model.rxnReferences{i};
        Model_MM.rxnReferences{end+1,1} = Model.rxnReferences{i};
    end
    if isfield(Model, 'rxnNotes')
        Model_MM.rxnNotes{end+1,1} = Model.rxnNotes{i};
        Model_MM.rxnNotes{end+1,1} = Model.rxnNotes{i};
    end
    if isfield(Model, 'confidenceScore')
        Model_MM.confidenceScore{end+1,1} = Model.confidenceScore{i};
        Model_MM.confidenceScore{end+1,1} = Model.confidenceScore{i};
    end
    if isfield(Model, 'proteins')
        Model_MM.proteins{end+1,1} = Model.proteins{i};
        Model_MM.proteins{end+1,1} = Model.proteins{i};
    end
end

for i=1:length(E)
    Model_MM.S(1:size(S_negative,1),end+1) = 0; % A + E = ES
    Model_MM.S(size(S_negative,1)+i,end) = -1; % add -E
    
    Model_MM.rev(end+1,1)=  1;
    Model_MM.c(end+1,1)  =  0;
    Model_MM.lb(end+1,1) = -1000;
    Model_MM.ub(end+1,1) =  1000;
    
    % rxn names
    Model_MM.rxns{end+1,1} = strcat(Model_MM.mets{size(S_negative,1)+i}, '_Import_Export');
    
    if isfield(Model, 'rxnNames')
        Model_MM.rxnNames{end+1,1} = strcat(Model_MM.mets{size(S_negative,1)+i}, '_Import_Export');
    end
    
    if isfield(Model, 'rxnECNumbers')
        Model_MM.rxnECNumbers{end+1,1} = '-';

    end
    if isfield(Model, 'rules')
        Model_MM.rules{end+1,1} = '-';
    end
    if isfield(Model, 'rxnGeneMat')
        Model_MM.rxnGeneMat(end+1,:) = 0;
    end
    if isfield(Model, 'grRules')
        Model_MM.grRules{end+1,1} = '-';
    end
    if isfield(Model, 'subSystems')
        Model_MM.subSystems{end+1,1} = '-';
    end
    if isfield(Model, 'rxnReferences')
        Model_MM.rxnReferences{end+1,1} = '-';
    end
    if isfield(Model, 'rxnNotes')
        Model_MM.rxnNotes{end+1,1} = '-';
    end
    if isfield(Model, 'confidenceScore')
        Model_MM.confidenceScore{end+1,1} = '-';
    end
    if isfield(Model, 'proteins')
        Model_MM.proteins{end+1,1} = '-';
    end
end

Model_MM.b = zeros(size(Model_MM.S,1),1);
    
EXITStatus=size(unique(Model_MM.S(1:size(S_positive,1)+length(E),:)','rows'))==size(Model_MM.S,2);

if EXITStatus~=1
    disp('Duplicated Enzyme Metabolite combinations found.')
end
end


