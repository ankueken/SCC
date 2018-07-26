function [Ref_flux,Ref_conc,is_sst] = perturbationAnalysis(percent,savename)
%% function to obtain different sst flux distributions from perturbed initial conditions

% load model
addpath('Marans_Ecoli_kinetic_core_model/Users/')
load('Marans_Ecoli_kinetic_core_model/Users/Data.mat')
is_sst=[];Ref_conc=[];Ref_flux=[];

global met_idx
met_idx=[];
[V,C]=solve_ode(Network_Data,0:100,Metabolite_and_Enzyme_Conc,Elementary_Kinetic_Param);
steady_state_concentration=C(:,end);
steady_state_diff_abs=max(abs(C(:,end-1)-C(:,end)));
steady_state_diff_rel=max(abs(1-(C(:,end-1)./C(:,end))));

% network metabolites 1:93 unexpanded, 94:231 enzymes, 232:830 enzyme complex

for run=1:50
    for met_idx=1:size(Network_Data.S_f_b,1)
        % pertubate metabolite concentration
        % a) plus x%
        disp([run met_idx])
        start=steady_state_concentration;
        start(met_idx) = steady_state_concentration(met_idx)+(steady_state_concentration(met_idx)*percent); % perturbe initial concentrations
        [Vplus,Cplus]=solve_ode(Network_Data,0:1000,start,Elementary_Kinetic_Param);
        r=0;
        while max(abs(1-(Cplus(:,end-1)./Cplus(:,end)))) > 1e-5 && r < 5
            [~,Cplus]=solve_ode(Network_Data,0:1000,Cplus(:,end-1),Elementary_Kinetic_Param);
            r=r+1;
        end
        is_sst(end+1) = max(abs(1-(Cplus(:,end-1)./Cplus(:,end)))); % check if  solution is steady state
        Ref_conc(:,end+1) = Cplus(:,end);
        Ref_flux(:,end+1) = Vplus(:,end);
        
        % b) minus x%
        start=steady_state_concentration;
        start(met_idx) = steady_state_concentration(met_idx)-(steady_state_concentration(met_idx)*percent);
        [Vminus(:,run),Cminus(:,run)]=solve_ode(Network_Data,0:1000,start,Elementary_Kinetic_Param);
        r=0;
        while max(abs(1-(Cminus(:,end-1)./Cminus(:,end)))) > 1e-5 && r < 5
            disp(r)
            [~,Cminus]=solve_ode(Network_Data,0:100,Cminus(:,end-1),Elementary_Kinetic_Param);
            r=r+1;
        end
        is_sst(end+1) = max(abs(1-(Cminus(:,end-1)./Cminus(:,end))));
        Ref_conc(:,end+1) = Cminus(:,end);
        Ref_flux(:,end+1) = Vminus(:,end);
    end
    save(savename)
end
end
