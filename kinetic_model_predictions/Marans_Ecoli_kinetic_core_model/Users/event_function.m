
function [value,isterminal,direction] = event_function(t_interval,y,CPU_t,kinetic_param,S_f_b)

% This function evaluates conservation of mass (for metabolites, enzymes and enzyme complexes) at each time interval and terminates ODE's
% integration once the model reaches to a steady-state solution (S.V~tol). 
% The function allows the integration to be continineud on CPU for 600 seconds by which it will be stopped (PMID:24928774).
%
% INPUTS:
% -------
%     t_interval: The time interval for simulations
%              y: Vector of metabolite and enzyme concentrations
%          CPU_t: CPU time
%  kinetic_param: Vector of elementary kinetic parameters
%          S_f_b: Stoichiometric matrix for the elementary rxns
%
%
% OUTPUTS:
% ------
%          value: Value of the event function
%     isterminal: 1 if the integration is to terminate at a zero of this event function, otherwise, 0
%      direction: 0 if all zeros are to be located (the default), +1 if only zeros where the event 
%                 function is increasing, and -1 if only zeros where the event function is decreasing
%
%
% Ali Khodayari, Ali R. Zomorrodi, Costas Maranas Lab @ Penn State

    % Elapsed CPU time
    c=t_interval-CPU_t;
    
    % Mass Conservation tolerance
    Error=max(abs(mass_balance_ode(0,y,kinetic_param,S_f_b)));
    direction = 0;
    tol=1;
    
    if (c<600 && Error>tol)
        value=Error;
        isterminal = 0;
    else
        value=0;
        isterminal = 1;
    end
end