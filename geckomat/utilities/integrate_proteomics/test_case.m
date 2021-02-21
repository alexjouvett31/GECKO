%
regime = {'low' 'high'};
for i = 1:length(regime)
    %substitute regime-specific dataset files in GECKO
    fileName = ['../../../Databases/abs_proteomics_' regime{i} '.txt'];
    copyfile(fileName,'../../../Databases/abs_proteomics.txt');
    fileName = ['../../../Databases/fermentationData_' regime{i} '.txt'];
    copyfile(fileName,'../../../Databases/fermentationData.txt');
    %load case-specific parameters
    [pIDs,protData,fermParameters,byProds] = load_Prot_Ferm_Data(3);
    %substitute parameters in getModelParameters
    
    %Check if model can achieve condition specific growth rate with the
    %specified GUR
    temp = setParam(ecModel_batch,'obj','r_2111',1);
    temp = setParam(temp,'ub','r_1714_REV',1.05*fermParameters.GUR);
    sol = solveLP(temp,1);
    if -sol.f<fermParameters.Drate
        cd ../../limit_proteins
        f = measureAbundance(temp.enzymes);
        cd ../kcat_sensitivity_analysis
        OptSigma = sigmaFitter(temp,fermParameters.Ptot,fermParameters.Drate,f);
        pbound = fermParameters.Ptot*OptSigma*f;
        temp = setParam(ecModel_batch,'ub','prot_pool_exchange',pbound);
        sol = solveLP(temp,1);
        printFluxes(temp,sol.x)
        cd ../utilities/integrate_proteomics
    end
    name = ['ecYeastGEM_prot_' regime{i}];
    generate_protModels(ecModel,3,name,temp)
end
 