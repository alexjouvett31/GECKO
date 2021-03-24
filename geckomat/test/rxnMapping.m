function mappedIndxs = rxnMapping(rxnID,rxnsList,revFlag)
% rxnMapping
% 
% Function that maps a metabolic rxn from a GEM on its enzyme constrained 
% version (irreversible model).
%
%    rxnID         rxnsList(i) that is going to be mapped
%    model         ecModel in which the rxn is going to be searched
%    revFlag       True if the searched rxn is reversible
%
%    mappedIndxs   Cell array that contains the index of the corresponding
%                  metabolic rxn in the ecModel (the arm rxn in the case of
%                  isoenzymes). If the original rxn is reversible then it
%                  contains the indexes of the forward and backward rxns
%                  respectively.
%
% Usage: mappedIndxs = rxnMapping(rxnID,model,revFlag)
%
% Ivan Domenzain.      Last edited: 2019-03-04


indexes = find(contains(rxnsList,rxnID));
if revFlag
    backwardIndxs = indexes(find(contains(rxnsList(indexes),'_REV')));
    forwardIndxs  = setdiff(indexes,backwardIndxs);
    backwardIndxs = findArmRxns(backwardIndxs,rxnsList);
    forwardIndxs  = findArmRxns(forwardIndxs,rxnsList);
    mappedIndxs   = vertcat(forwardIndxs,backwardIndxs);
else
    mappedIndxs  = findArmRxns(indexes,rxnsList);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ArmIndex = findArmRxns(rxnIndexes,rxnsList)
if length(rxnIndexes)>1
    idx = find(contains(rxnsList(rxnIndexes),'arm_'));
    if ~isempty(idx)
        ArmIndex = rxnIndexes(contains(rxnsList(rxnIndexes),'arm_'));
    else
        ArmIndex = rxnIndexes;
    end
else
    ArmIndex = rxnIndexes;
end
if isempty(ArmIndex)
    rxnsList(rxnIndexes)
    pause
end
end