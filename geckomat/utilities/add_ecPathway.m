function newModel = add_ecPathway(ecModel,fileName)
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
% Ivan Domenzain  Last edited: 2020-06-20

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
    brackets = strfind(met,'[');
    brackets = brackets(end);
    comp = met(brackets+1:end-1);
    met  = met(1:brackets-1);
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
%Add enzymes
enzymes = pathwayTable.proteins;
kcats   = pathwayTable.kcats;
MWs     = pathwayTable.MWs;
for i=1:height(pathwayTable)
    rxnIdx = strcmpi(newModel.rxns,pathwayTable.rxns(1));
    if ~isempty(enzymes{i})
        if ~ismember(enzymes{i},newModel.enzymes)
            newModel.enzymes  = [newModel.enzymes; enzymes(i)];
            newModel.MWs      = [newModel.MWs; MWs(i)];
            newModel.enzGenes = [newModel.enzGenes; pathwayTable.grRules(i)];
            newModel.enzNames = [newModel.enzNames; pathwayTable.grRules(i)];
            newModel.pathways = [newModel.pathways; {''}];
            newModel.concs = [newModel.concs; NaN];
            metsToAdd = [];
            metId     = ['prot_' enzymes{i}];
            metsToAdd.mets         = {metId};
            metsToAdd.metNames     = {metId};
            metsToAdd.compartments = {'c'};
            newModel = addMets(newModel,metsToAdd);
        end
        metIdx = strcmpi(newModel.mets,metId);
        %assign Kcat as pseudo-stoichiometric coeff.
        newModel.S(metIdx,rxnIdx) = -1/(3600*kcats(i));
        %Add protein draw reaction
        rxnsToAdd = [];
        rxnsToAdd.rxns      = {['draw_' metId]};
        rxnsToAdd.rxnNames  = {['draw_' metId]};
        rxnsToAdd.grRules   = pathwayTable.grRules(i);
        rxnsToAdd.equations = {[num2str(MWs(i)) ' prot_pool[c] => ' metId '[c]']};
        rxnsToAdd.c  = 0;
        rxnsToAdd.lb = 0;
        rxnsToAdd.ub = 1000;
        newModel = addRxns(newModel,rxnsToAdd,3);
    end
end    
end
