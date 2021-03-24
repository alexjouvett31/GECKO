function kcats_table = find_kcat_discrepancies(ecModel)
data = readtable('../../Databases/manual_data.txt','delimiter','\t');
kcatStrs = [];
rxnStrs  = [];
rxnNames = [];
grRules  = [];
cd ../utilities
for i=1:height(data)
    protein = data.Var2{i};
    if ~contains(protein,' + ')
        [kcat,rxnIdx,rxnName] = getKcat(ecModel,protein);
        kcat = num2cell(kcat);
        kcat = cellfun(@num2str,kcat,'UniformOutput',false);
        kcat = strjoin(kcat,'//');
        
        rxns = ecModel.rxns(rxnIdx);
        rxns = strjoin(rxns,'//');
        
        rxnName = strjoin(rxnName,'//');
        
        rule = ecModel.grRules(rxnIdx);
        rule = strjoin(rule,'//');
        
        kcatStrs = [kcatStrs;{kcat}];
        rxnStrs  = [rxnStrs;{rxns}];
        rxnNames = [rxnNames;{rxnName}];
        grRules  = [grRules;{rule}];
    else
        kcatStrs = [kcatStrs;{''}];
        rxnStrs  = [rxnStrs;{''}];
        rxnNames = [rxnNames;{''}];
        grRules  = [grRules;{''}];
    end
end
cd ../test
kcats_table = data;
kcats_table(:,6) = [];
kcats_table.Properties.VariableNames = {'name' 'proteins' 'ecNumber' 'genes' 'kcat'};
kcats_table.kcat_model = kcatStrs;
kcats_table.grRules = grRules;
kcats_table.rxnNames = rxnNames;
kcats_table.rxns = rxnStrs;


end