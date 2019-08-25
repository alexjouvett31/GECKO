%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,ecModelSim] = getECDrawModel(model,org_name,model_name,version,Ptot)
% Main function for generating an enzyme-constrained draw model.
%
%   model           Model following the RAVEN structural convention
%   sigma           Average saturation factor. Default = 0.5.
%   Ptot            Total protein content [g/gDW]. Default = 0.67.
% 	org_name        Name of your organism (e.g. 'homo sapiens')
% 	name            Name to assign the finished model (e.g. 'ec_human')
%   version         Numeric identifyer of model version (e.g. '1.0.1')
%	
%   ecModel         Enzyme-constrained draw model
%   ecModelSim      Simulatable EC-draw model 
%
% Daniel Cook       08/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecModel,ecModelSim] = getECDrawModel(model,org_name,name,version,Ptot)

%Convert model to RAVEN for easier visualization later on:
format short e
if isfield(model,'rules')
    initCobraToolbox
    model = ravenCobraWrapper(model);
end

%Remove blocked rxns + correct model.rev:
cd change_model
[model,name,version] = preprocessModel(model,name,version);

%Retrieve kcats & MWs for each rxn in model:
cd ../get_enzyme_data
model_data = getEnzymeCodes(model,'add');
kcats      = matchKcats(model_data,org_name);

%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);
% [ecModel,modifications] = manualModifications(ecModel); % No manual modifications

% Constrain protein pool
cd ../limit_proteins
total_protein_mass = Ptot;
non_measured = ones(length(ecModel.enzymes),1);
ecModel = constrainPool(ecModel,non_measured,total_protein_mass);

% Make simulatable
ecModelSim = simplifyModel(ecModel,true,false,true,true);
end

