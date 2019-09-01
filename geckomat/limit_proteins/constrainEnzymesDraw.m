%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model,sol,enzUsages,modifications,MW,counter] = constrainEnzymesDraw(model,Ptot,sigma,f,pIDs,data,MW,counter)
% Main function for overlaying proteomics data on an enzyme-constrained
% model. 
%
%   model           ecModel.
%   sigma           Average saturation factor. Default = 0.5.
%   Ptot            Total protein content [g/gDW]. Default = 0.67.
% 	f				(Opt) Estimated mass fraction of enzymes in model.
% 	pIDs			(Opt) Protein IDs from proteomics data.
%	data			(Opt) Raw protein abundances from proteomics data [counts or abundance].
%   MW              (Opt) MW for all proteins in pIDs
%   counter         (Opt) Match between model.proteins & pIDs
%
%   model           ecModel with calibrated enzyme usage upper bounds
%   enzUsages       Calculated enzyme usages after final calibration 
%                   (enzyme_i demand/enzyme_i upper bound)
%   modifications   Table with all the modified values 
%                   (Protein ID/old value/Flexibilized value)
%   MW              MW for all proteins in pIDs (calculated if not provided)
%   counter         Match between model.proteins & pIDs (calculated if not provided)
%
% Daniel Cook       08/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,sol,MW_out,counter_out] = constrainEnzymesDraw(model,Ptot,sigma,f,pIDs,data,MW,counter)

%Compute f if not provided:
if nargin < 4
    [f,~] = [];
end

%No UB will be changed if no data is available -> pool = all enzymes(FBAwMC)
if nargin < 5
    pIDs = cell(0,1);
    data = zeros(0,1);
end

%Remove negative values
for i=1:length(data)
    if data(i)< 0
        data(i) = NaN;
    end
end

if sum(isnan(data)) ~= 0
    fprintf('Warning: data contains negative numbers.')
end

if ~exist('MW') || ~exist('counter') % Check this statement
    databases = load('../../databases/ProtDatabase.mat');
    swissprot = databases.swissprot;
    for i = 1:length(swissprot)
        swissprot{i,1} = strsplit(swissprot{i,1},' ');
    end
    
    %Main loop: grab MW for all proteins in dataset
    MW_ave  = mean(cell2mat(swissprot(:,5)));
    concs   = zeros(size(pIDs));
    MW   = zeros(size(pIDs));
    counter = false(size(pIDs));
    
    for i = 1:size(swissprot,1)
        swissprot_vector(i) = swissprot{i,1};
    end
    
    % Match pIDs to swissprot database
    [matched_proteins,swissprot_indx] = ismember(pIDs,swissprot_vector);
    
    % Get MW
    MW(matched_proteins ~= 0) = swissprot_indx(swissprot_indx ~= 0);
    MW(matched_proteins == 0) = MW_ave;
    
    % Identify proteins in the model (N.B. This variable should be renamed.)
    counter = ismember(pIDs,model.proteins);
    
    % Assign output variables
    MW_out = MW;
    counter_out = counter;
    
%     %Main loop: grab MW for all proteins in dataset
%     MW_ave  = mean(cell2mat(swissprot(:,5)));
%     concs   = zeros(size(pIDs));
%     MW   = zeros(size(pIDs));
%     counter = false(size(pIDs));
%     for i = 1:length(pIDs)
%         MW(i) = MW_ave;
%         %Find gene in swissprot database:
%         for j = 1:length(swissprot)
%             if sum(strcmp(swissprot{j,1},pIDs{i})) > 0
%                 MW(i) = swissprot{j,5};	%g/mol
%                 %Check if uniprot is in model:
%                 if sum(strcmp(model.proteins,swissprot{j,1})) > 0
%                     counter(i) = true;
%                 end
%             end
%         end
%         if rem(i,100) == 0
%             disp(['Calculating total abundance: Ready with ' num2str(i) '/' ...
%                 num2str(length(pIDs)) ' genes '])
%         end
%     end
% end

% Set up data
total_protein_mass = Ptot; % User input
concs = MW.*data;     %g/mol(tot prot)
concs_sum = sum(concs(~isnan(concs)));
mass_fracs = concs/concs_sum;

% Constrain model fluxes
for i =1:length(model.rxns)
    temp = strsplit(model.rxns{i},'draw_prot_');
    if numel(temp) == 1
        temp{2}='';
    end
    prot_draw{i} = temp{2};
end
[~, loc]=ismember(pIDs,prot_draw); % match indexes

j = find(loc~=0);
[~,loc2] = ismember(pIDs(j),model.enzymes);
model.ub(loc(j)) = mass_fracs(j)*total_protein_mass*sigma./model.MWs(loc2);

% Unconstrain protein fluxes set to zero
model.ub(intersect(loc(j),find(model.ub==0))) = min(model.ub(intersect(loc(j),find(model.ub~=0))));

% Calculate metabolic protein fraction
if ~exist(f)
    f_resid = sum(mass_fracs(counter==1));
else
    f_resid = f; % User-specified value
end
model.ub(strcmp('prot_pool_exchange',model.rxns)) = total_protein_mass*f_resid*sigma;

% Esure the constrained model is solvable
sol = solveLP(model);
if isempty(sol.f)
    [sol,gR,relaxed_index] = ...
        relax_constraints(model,pIDs,data,MW,loc,counter,total_protein_mass,sigma);
end
end

