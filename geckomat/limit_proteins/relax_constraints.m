function [sol,gR,relaxed_index] = relax_constraints(model,pIDs,data,MW,loc,counter,total_protein_mass,sigma)
% Run this function when the growth solver returns an infeasible solution.
% This function checks to see if there is a single index keeping the model
% from solving the objective function. 

% % Solve once and check to see if the model runs
% concs = MW.*data;     %g/mol(tot prot)
% concs_sum = sum(concs(~isnan(concs)));
% mass_fracs = concs/concs_sum;
% 
% % Constrain fluxes
% model_batch = model;
% j = find(loc~=0);
% [~,loc2] = ismember(pIDs(j),model_batch.enzymes);
% model_batch.ub(loc(j)) = mass_fracs(j)*total_protein_mass*sigma./model_batch.MWs(loc2);
% 
% % Unconstrain protein fluxes set to zero
% model_batch.ub(intersect(loc(j),find(model_batch.ub==0))) = min(model_batch.ub(intersect(loc(j),find(model_batch.ub~=0))));
% 
% % Calculate metabolic protein fraction
% f_resid = sum(mass_fracs(counter==1));
% model_batch.ub(strcmp('prot_pool_exchange',model_batch.rxns)) = total_protein_mass*f_resid*sigma;
% 
% % Solve for maximizing growth
% sol = solveLP(model_batch);
% 
% % Search for a single reaction causing an issue with the solver
% if isempty(sol.f)
%     
%     indxs_ok = [];
%     solved = 0;
%     
%     while solved == 0
%         %Split into two
%         indxs = find(loc~=0);
%         
%         indxs_split = indxs(~ismember(indxs,indxs_ok));
%         
%         indxs1 = indxs_split(1:(round(length(indxs_split)/2)));
%         indxs2 = indxs_split((round(length(indxs_split)/2)+1):end);
%         
%         model_batch = model;
%         
%         % Constrain fluxes
%         j = [indxs_ok;indxs1];
%         [~,loc2] = ismember(pIDs(j),model_batch.enzymes);
%         model_batch.ub(loc(j)) = mass_fracs(j)*total_protein_mass*sigma./model_batch.MWs(loc2);
%         
%         % Unconstrain protein fluxes set to zero
%         model_batch.ub(intersect(loc(j),find(model_batch.ub==0))) = min(model_batch.ub(intersect(loc(j),find(model_batch.ub~=0))));
%         
%         % Calculate metabolic protein fraction
%         f_resid = sum(mass_fracs(counter==1));
%         model_batch.ub(strcmp('prot_pool_exchange',model_batch.rxns)) = total_protein_mass*f_resid*sigma;
%         
%         % Solve for maximizing growth
%         sol = solveLP(model_batch);
%         if ~isempty(sol.f)
%             indxs_ok = [indxs_ok;indxs1];
%         else
%             indxs_ok = [indxs_ok;indxs2];
%         end
%         
%         if length(indxs(~ismember(indxs,indxs_ok))) == 1
%             solved = 1;
%             relaxed_index = indxs(~ismember(indxs,indxs_ok));
%             relaxed_flux = sol.x(relaxed_index);
%             gR = -sol.f;
%         end
%     end
% else
%     disp("Solution found. Exiting.");
% end
% end

% Solve once and check to see if the model runs
concs = MW.*data;     %g/mol(tot prot)
concs_sum = sum(concs(~isnan(concs)));
mass_fracs = concs/concs_sum;

% Constrain fluxes
model_batch = model;
j = find(loc~=0);
[~,loc2] = ismember(pIDs(j),model_batch.enzymes);
model_batch.ub(loc(j)) = mass_fracs(j)*total_protein_mass*sigma./model_batch.MWs(loc2);

% Unconstrain protein fluxes set to zero
model_batch.ub(intersect(loc(j),find(model_batch.ub==0))) = min(model_batch.ub(intersect(loc(j),find(model_batch.ub~=0))));

% Calculate metabolic protein fraction
f_resid = sum(mass_fracs(counter==1));
model_batch.ub(strcmp('prot_pool_exchange',model_batch.rxns)) = total_protein_mass*f_resid*sigma;

% Solve for maximizing growth
sol = solveLP(model_batch);

% Search for a single reaction causing an issue with the solver
if isempty(sol.f)
    
    indxs_ok = [];
    solved = 0;
    
    %Split constrained indexes into batches of 100
    indxs = find(loc~=0);
    batches = 1:100:length(indxs);
    batches(end+1) = length(indxs);
    
    for i = 1:(length(batches)-1)
        model_batch = model;
        
        % Constrain fluxes
        j = [indxs_ok;indxs(batches(i):batches(i+1))];
        [~,loc2] = ismember(pIDs(j),model_batch.enzymes);
        model_batch.ub(loc(j)) = mass_fracs(j)*total_protein_mass*sigma./model_batch.MWs(loc2);
        
        % Unconstrain protein fluxes set to zero
        model_batch.ub(intersect(loc(j),find(model_batch.ub==0))) = min(model_batch.ub(intersect(loc(j),find(model_batch.ub~=0))));
        
        % Calculate metabolic protein fraction
        f_resid = sum(mass_fracs(counter==1));
        model_batch.ub(strcmp('prot_pool_exchange',model_batch.rxns)) = total_protein_mass*f_resid*sigma;
        
        % Solve for maximizing growth
        sol = solveLP(model_batch);
        if ~isempty(sol.f)
            indxs_ok = j;
        end
    end
    
    %Split constrained indexes into batches of 10
    indxs_10 = indxs(~ismember(indxs,indxs_ok));
    batches = 1:10:length(indxs_10);
    batches(end+1) = length(indxs_10);
    
    for i = 1:(length(batches)-1)
        model_batch = model;
        
        % Constrain fluxes
        j = [indxs_ok;indxs_10(batches(i):batches(i+1))];
        [~,loc2] = ismember(pIDs(j),model_batch.enzymes);
        model_batch.ub(loc(j)) = mass_fracs(j)*total_protein_mass*sigma./model_batch.MWs(loc2);
        
        % Unconstrain protein fluxes set to zero
        model_batch.ub(intersect(loc(j),find(model_batch.ub==0))) = min(model_batch.ub(intersect(loc(j),find(model_batch.ub~=0))));
        
        % Calculate metabolic protein fraction
        f_resid = sum(mass_fracs(counter==1));
        model_batch.ub(strcmp('prot_pool_exchange',model_batch.rxns)) = total_protein_mass*f_resid*sigma;
        
        % Solve for maximizing growth
        sol = solveLP(model_batch);
        if ~isempty(sol.f)
            indxs_ok = j;
        end
    end
    
    %Split constrained indexes into batches of 1
    indxs_1 = indxs(~ismember(indxs,indxs_ok));
    batches = 1:1:length(indxs_1);
    
    for i = 1:(length(batches))
        model_batch = model;
        
        % Constrain fluxes
        j = [indxs_ok;indxs_1(batches(i))];
        [~,loc2] = ismember(pIDs(j),model_batch.enzymes);
        model_batch.ub(loc(j)) = mass_fracs(j)*total_protein_mass*sigma./model_batch.MWs(loc2);
        
        % Unconstrain protein fluxes set to zero
        model_batch.ub(intersect(loc(j),find(model_batch.ub==0))) = min(model_batch.ub(intersect(loc(j),find(model_batch.ub~=0))));
        
        % Calculate metabolic protein fraction
        f_resid = sum(mass_fracs(counter==1));
        model_batch.ub(strcmp('prot_pool_exchange',model_batch.rxns)) = total_protein_mass*f_resid*sigma;
        
        % Solve for maximizing growth
        sol = solveLP(model_batch);
        if ~isempty(sol.f)
            indxs_ok = j;
        end
    end
    
    % Solve with trouble rxn bounds removed
    model_batch = model;
    j = indxs_ok;
    [~,loc2] = ismember(pIDs(j),model_batch.enzymes);
    model_batch.ub(loc(j)) = mass_fracs(j)*total_protein_mass*sigma./model_batch.MWs(loc2);
    
    % Unconstrain protein fluxes set to zero
    model_batch.ub(intersect(loc(j),find(model_batch.ub==0))) = min(model_batch.ub(intersect(loc(j),find(model_batch.ub~=0))));
    
    % Calculate metabolic protein fraction
    f_resid = sum(mass_fracs(counter==1));
    model_batch.ub(strcmp('prot_pool_exchange',model_batch.rxns)) = total_protein_mass*f_resid*sigma;
    
    % Solve for maximizing growth
    sol = solveLP(model_batch);
    
    % Create outputs
    relaxed_index = indxs(~ismember(indxs,indxs_ok));
%     relaxed_flux = sol.x(relaxed_index);
    gR = -sol.f;
else
    disp("Solution found. Exiting.");
end
end





