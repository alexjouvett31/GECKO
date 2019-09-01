%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecCorrModel] = getECCorrModel(ecModelsim)
% Main function for generating an enzyme-constrained draw model with
% protein source and sink reactions to use for simulating based on
% maximizing correlation between measured proteomics and model protein use.
%
%   ecModelsim      A simulatable enzyme-constrained-draw model
%	
%   ecCorrModel     An enzyme-constrained draw model that contains source
%                   and sink reactions for each protein and an objective
%                   function to minimize use of source and sink reactions
%   
% Daniel Cook       08/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecCorrModel] = getECCorrModel(ecModelsim)

% Clear objective function and set new objective 
% (obj: minimize deviation from data)
model = ecModelsim;
newModel = ecModelsim;
newModel.c(:) = 0;
for i =1:length(model.rxns)
    temp = strsplit(model.rxns{i},'draw_prot_');
    if numel(temp) == 2
        rxnsToAdd.rxns = {strcat('prot_',temp{2},'_sink'),strcat('prot_',temp{2},'_source')};
        rxnsToAdd.equations = cellstr([strcat(model.mets(find(model.S(:,i) > 0))," => prot_sink");...
            strcat("prot_source => ",model.mets(find(model.S(:,i) > 0)))]);
        rxnsToAdd.rxnNames = rxnsToAdd.rxns;
        rxnsToAdd.lb = [0 0];
        rxnsToAdd.ub = [1000 1000];
        rxnsToAdd.c = [0,0];
        rxnsToAdd.rxnComps = ["s","s"];
        newModel=addRxns(newModel,rxnsToAdd,1,'s',1);
    end
    if numel(temp) == 1
        temp{2}='';
    end
    prot_draw{i} = temp{2};
end

% Add one source and one sink rxn
rxnsToAdd.rxns = {'prot_sink','prot_source'};
rxnsToAdd.equations = cellstr(["prot_sink => ";" => prot_source"]);
rxnsToAdd.rxnNames = rxnsToAdd.rxns;
rxnsToAdd.lb = [0 0];
rxnsToAdd.ub = [1000 1000];
rxnsToAdd.c = [-1,-1];
rxnsToAdd.rxnComps = ["s","s"];
newModel=addRxns(newModel,rxnsToAdd,1,'s',1);

% Set output
ecCorrModel = newModel;
end