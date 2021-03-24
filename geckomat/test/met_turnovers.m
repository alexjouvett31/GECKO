function metTurnovers = met_turnovers(fluxTable,model)
%Get model S matrix
Smat = model.S;
%Just keep those reactions that are also part of the fluxTable
[iA,iB] = ismember(fluxTable.rxns,model.rxns);
if sum(iA) == length(iA)
    Smat = Smat(:,iB);
end
%compute met turnover numbers 
[m,~] = size(Smat);
turnoverNumbers = zeros(m,1);
for i=1:m
    turnoverNumbers(i) = 0.5*(sum(abs(Smat(i,:)'.*fluxTable.netFluxes)));     
end
compartments = model.compNames(model.metComps);
metTurnovers    = table(model.mets,model.metNames,compartments,turnoverNumbers);
end
