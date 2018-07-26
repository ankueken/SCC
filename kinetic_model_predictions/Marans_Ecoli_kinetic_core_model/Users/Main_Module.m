%------------------------------------------------------------
% This code integrates the ODEs of mass conservation for the metabolites and
% enzymes (free and complexes), following decomposition of reactions into
% their elementary stpes (PMID:24928774).
%
%
% INPUTS:
% -------
%          x1: The id of the rection(s) that gets perturbed
%          x2: The level of the perturbation as follows:
%                       0: knock out
%             less than 1: down regulation
%          greater than 1: up regulation (e.g., 2 means 2-fold overexpression)
%               
% 
% OUTPUTS:
% ------
%           Vnet_perturb: Rate of the reactions at different time intervals
%  concentration_perturb: Metabolite concentrations at different time intervals
%
%
% Ali Khodayari, Ali R. Zomorrodi, Costas Maranas Lab @ Penn State
%------------------------------------------------------------

function [Vnet_perturb,concentration_perturb]=Main_Module(x1,x2)

% Loading model structure and parameters
load Data.mat
    
% Time interval for ODEs integration
t_interval = [0:1e4];


if nargin<2;
    
   fprintf('\n-------------------------------\n');
   fprintf('\nTwo inputs are required to run this module, in form of:\n');
   fprintf('Main_Module(reaction ID, Perturbation level)\n');
   fprintf('Default: reaction ID=1, Perturbation level=1\n');
   x1=1;
   x2=1;
   
end

fprintf('\n-------------------------------\n');

       % Initial metabolite and enzyme ocncentrations
       y0 = perturbModel(Network_Data,Metabolite_and_Enzyme_Conc,x1,x2);
       
       % Integration of ODEs
       [Vnet_perturb,concentration_perturb] = solve_ode(Network_Data,t_interval,y0,Elementary_Kinetic_Param);

end
