function ecModel = add_ecPathway(ecModel,fileName)
% add_ecPathway
%
%
%   ecModel         enzyme-constrained model
%
%   Output:
%   newModel        enzyme-constrained model with the introduced
%                   heterologous pathway
%
% usage: [kcat,rxnIdx,rxnName,MW] = add_ecPathway(ecModel,protein)
%
% Ivan Domenzain  Last edited: 2019-11-23

%Load pathway file as a table
fileName = ['../../databases/' fileName];
if exist(fileName,'file')~= 0
    pathwayTable = readtable(fileName,'delimiter','\t');
end
%Define rxnsToAdd structure
rxnsToAdd.rxns      = pathwayTable.rxns;
rxnsToAdd.rxnNames  = pathwayTable.rxnNames;
rxnsToAdd.grRules   = pathwayTable.grRules;
rxnsToAdd.equations = pathwayTable.formulas;
rxnsToAdd.c  = pathwayTable.c;
rxnsToAdd.lb = pathwayTable.lb;
rxnsToAdd.ub = pathwayTable.ub;
%Get metNames from formulas
[~, mets, ~, ~]=constructS(pathwayTable.formulas);
%Check which metabolites are not in model and incorporate them
for i=1:length(mets)
    metsToAdd = [];
    met  = mets{i};
    comp = met(end-2:end);
    met  = met(1:end-3);
    comp = strrep(comp,'[','');
    comp = strrep(comp,']','');
    comp = find(strcmp(ecModel.comps,comp));
    metIndxs = find(strcmpi(ecModel.metNames,met));
    if ~isempty(metIndxs)
        %If metName is already in model but not found in the same
        %compartment then it should be added
        if ~any(ecModel.metComps(metIndxs)==comp)
            metId = generateNewIds(ecModel,'mets','s_',1);
            metsToAdd.mets         = metId;
            metsToAdd.metNames     = {met};
            metsToAdd.compartments = ecModel.comps(comp);
        end
    else
        metId = generateNewIds(ecModel,'mets','s_',1);
        metsToAdd.mets         = metId;
        metsToAdd.metNames     = {met};
        metsToAdd.compartments = ecModel.comps(comp);
    end
    if ~isempty(metsToAdd)
        ecModel = addMets(ecModel,metsToAdd);
    end
end
%Same with genes
for i=1:height(pathwayTable)
    grRule = pathwayTable.grRules{i};
    if ~isempty(grRule)
        grRule = strrep(grRule,'(','');
        grRule = strrep(grRule,')','');
        grRule = strrep(grRule,' or ',' ');
        grRule = strrep(grRule,' and ',' ');
        grRule = strsplit(grRule,' ');
        for gene=grRule
            if ~ismember(gene,ecModel.genes)
                genesToAdd.genes = gene;
                genesToAdd.geneShortNames = gene;
                ecModel = addGenesRaven(ecModel,genesToAdd);
            end
        end
    end
end
%Add rxns
newModel = addRxns(ecModel,rxnsToAdd,3);
%Standardize gene related fields
[grRules, rxnGeneMat] = standardizeGrRules(newModel,true);
newModel.grRules      = grRules;
newModel.rxnGeneMat   = rxnGeneMat;
end
