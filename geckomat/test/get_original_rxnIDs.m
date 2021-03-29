function original_IDs = get_original_rxnIDs(rxnIDs)
original_IDs = unique(rxnIDs);
original_IDs = strrep(original_IDs,'_REV','');
original_IDs = strrep(original_IDs,'arm_','');
original_IDs = strrep(original_IDs,'_arm_','');
original_IDs = regexprep(original_IDs,'No(\d)','');
original_IDs = unique(original_IDs(~contains(original_IDs,'prot_')));
original_IDs = strtrim(original_IDs);

original_IDs = unique(original_IDs);

end