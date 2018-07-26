function [T,conc] = solve_ode_example(Model,t_interval,y0,kinetic_param)

S_f_b = Model.S;

[T,Y]=ode45(@(t,y)mass_balance_ode_example(t,y,kinetic_param,S_f_b),t_interval,y0);
conc = Y'; 
end
