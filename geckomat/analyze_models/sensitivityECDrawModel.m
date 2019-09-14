%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sensRes] = sensitivityECDrawModel(ecModelsim,data,Ptot,sigma,MW,counter,analysis)
% First run constrainEnzymesDraw.m to generate a constrained, simulatable
% model
% (EDIT BELOW)
% Main function for generating an enzyme-constrained draw model with
% protein source and sink reactions to use for simulating based on
% maximizing correlation between measured proteomics and model protein use.
%
%   ecModelsim      A simulatable enzyme-constrained-draw model
%	
%   ecCorrModel     An enzyme-constrained draw model that contains source
%                   and sink reactions for each protein and an objective
%                   function to minimize use of source and sink reactions
%   
% Daniel Cook       08/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sensRes] = sensitivityECDrawModel(ecModelsim,data,Ptot,sigma,MW,counter,loc,analysis)
%% Set up function inputs
%If no single analysis is selected, run all
if nargin < 5
    MW = [];
    counter = [];
    loc = [];
end

%If no single analysis is selected, run all
if nargin < 8
    analysis = [1,2,3,4];
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

% Set up model
model = ecModelsim;

% Get MW and counter if not supplied
if isempty('MW') || isempty('counter') % Check this statement
    databases = load('../databases/ProtDatabase.mat');
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
end

%% Sensitivity Analysis 1: Enzyme knockdown
if ismember(1,analysis)
    % Save UB from original model
    ub_orig = model.ub;
    % Get location of proteins in model (if not provided)
    if isempty(loc)
        for i =1:length(model.rxns)
            temp = strsplit(model.rxns{i},'draw_prot_');
            if numel(temp) == 1
                temp{2}='';
            end
            prot_draw{i} = temp{2};
        end
        [~, loc]=ismember(pIDs,prot_draw); % match indexes
    end
    
    % Perform sensitivity analysis
    j = find(loc~=0);
    for i = 1:length(j)
        model.ub = ub_orig;
        model.ub(loc(j(i))) = 0; % knock-out each enzyme
        sol = solveLP(model);
        if ~isempty(sol.f)
            gR(i) = -sol.f;
        else
            gR(i) = NaN;
        end
    end
    sensRes.EKO.rxns = model.rxns(loc(j));
    sensRes.EKO.gR = gR;
end
%% Sensitivity Analysis 2: Metabolite inhibition (anti-metabolite)
if ismember(2,analysis)
    for i = 1:length(model.mets)
        rxns = find(model.S(i,:)); % find all rxns for the metabolite
        rxns_sub = find(model.S(i,rxns)<0); % find all rxns where metabolite is a substrate
        model.ub(rxns_sub) = 0; % Set UB to 0
        sol = solveLP(model);
        if ~isempty(sol.f)
            gR(i) = -sol.f;
        else
            gR(i) = NaN;
        end
    end
    sensRes.MKO.met = model.mets;
    sensRes.MKO.gR = gR;
end

%% Sensitivity Analysis 3: Enzyme treatment
if ismember(3,analysis)
    % Save UB from original model
    ub_orig = model.ub;
    lb_orig = model.lb;
    % Get location of proteins in model (if not provided)
    if isempty(loc)
        for i =1:length(model.rxns)
            temp = strsplit(model.rxns{i},'draw_prot_');
            if numel(temp) == 1
                temp{2}='';
            end
            prot_draw{i} = temp{2};
        end
        [~, loc]=ismember(pIDs,prot_draw); % match indexes
    end
    
    % Perform sensitivity analysis
    j = find(loc~=0);
    for i = 1:length(j)
        model.ub = ub_orig;
        model.lb = lb_orig;
        model.ub(loc(j(i))) = 3*ub_orig(loc(j(i))); % Push each enzyme
        model.lb(loc(j(i))) = 2*ub_orig(loc(j(i))); % Push each enzyme
        sol = solveLP(model);
        if ~isempty(sol.f)
            gR(i) = -sol.f;
        else
            gR(i) = NaN;
        end
    end
    sensRes.EP.rxns = model.rxns(loc(j));
    sensRes.EP.gR = gR;
end

%% Sensitivity Analysis 4: Metabolite treatment
if ismember(4,analysis)
    % Should this be something with the b vector?
    b_orig = model.b;
    met_end = (min(find(startsWith(model.mets,'pmet'),1),find(startsWith(model.mets,'prot'),1))-1);
    for i = 1:met_end % End of metabolites, start of proteins
        model.b = b_orig;
        model.b(1,i) = .01;
        sol = solveLP(model);
        if ~isempty(sol.f)
            gR(i) = -sol.f;
        else
            gR(i) = NaN;
        end
    end
    sensRes.MP.mets = model.mets(1:met_end);
    sensRes.MP.gR = gR;
end





end