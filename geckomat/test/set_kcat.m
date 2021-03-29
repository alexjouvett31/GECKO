function model = set_kcat(model,protein,kcat)
cd  ../utilities
[~,rxnIdx] = getKcat(model,protein);
metIdx = find(strcmpi(model.metNames,['prot_' protein]));
kcat = -1/(kcat*3600);
model.S(metIdx,rxnIdx) = kcat;
cd ../test
end