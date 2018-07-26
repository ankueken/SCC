addpath('Marans_Ecoli_kinetic_core_model/Users/')
load('Marans_Ecoli_kinetic_core_model/Users/Data.mat')
[~,C,~]=solve_ode(Network_Data,0:100,Metabolite_and_Enzyme_Conc,Elementary_Kinetic_Param);

%% WT reference
% we use initial conditions provided by the model to simulate a steady-state
% flux distribution

SST = ones(size(Network_Data_mutant.S_f_b),2,1);

for i=1:size(Network_Data_mutant.S_f_b,2)
	% KO reaction
	Network_Data_mutant = Network_Data;
	Network_Data_mutant.S_f_b(:,i)=0;

	% same initial conditions as used to obtain WT
	[~,Cmt,V]=solve_ode(Network_Data_mutant,0:100,C(:,end),Elementary_Kinetic_Param);

	if max(abs(Cmt(:,end-1)-Cmt(:,end))) > 1e-6 && size(Cmt,2)==101
	    SST(i,1)=0;
	    [~,Cmt,V]=solve_ode(Network_Data_mutant,0:1000,C(:,end),Elementary_Kinetic_Param);
	    if max(abs(Cmt(:,end-1)-Cmt(:,end))) > 1e-6 && size(Cmt,2)==1001
		SST(i,1)=0;
	    else
		SST(i,1)=1;
	    end
	end

	% fold-changes between WT and simulated mutant from kinetic model (observed fold changes)
	CmA(i,:) = Cmt(:,end)';
	V_ko(:,i)=V(:,end);
end
save('SimFoldChanges_All.mat')

