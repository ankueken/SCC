% Input: percent: numeric value between 0 and 1
%        RS: indices of relevant rate constants (can be obtained from
%        structural_robustness_kcat_lp output "ReactionSet")
% Output: correlation of predicted and simulated ranges over 100 random
% removals of SCC relevant kcat values
function removal_kcat_knowledge(percent,RS,savename)
addpath('Marans_Ecoli_kinetic_core_model/Users/')
load('Marans_Ecoli_kinetic_core_model/Users/Data.mat')

Model.S = Network_Data.S_f_b;
Model.lb = ones(size(Network_Data.S_f_b,2),1)*1e-6;
Model.ub = ones(size(Network_Data.S_f_b,2),1)*1000;
Model.rev = zeros(size(Network_Data.S_f_b,2),1);
Model.rxns = Network_Data.rxn_f_b;
Model.mets = Network_Data.metab;
Model.c = zeros(size(Network_Data.rxn_f_b));
Model.b = zeros(size(Network_Data.metab));

%% simulate missing knowledge of rate constants

for run=1:100
    % full table
    kcat=table(Model.rxns,cell(size(Model.rxns)),cell(size(Model.rxns)),Elementary_Kinetic_Param);
    % mi
    kcat{RS(randperm(length(RS),ceil(length(RS)*percent))),4} = NaN;
    cd ..
    [SCC, Status, B_min, B_max, ReactionSet, ODE] = get_SCC_kcat(Model,[],[],kcat,'equal',1:length(Model.mets));
    cd kinetic_model_predictions
    %% concentration from flux sampling
    
    % load simulation results
    R1=load('Perturb_all_1pc_1.mat');
    R2=load('Perturb_all_5pc_1.mat');
    R3=load('Perturb_all_10pc_1.mat');
    
    t=1e-2;
    Ref_flux_x = [R1.Ref_flux(:,find(R1.is_sst<t)) R2.Ref_flux(:,find(R2.is_sst<t)) R3.Ref_flux(:,find(R3.is_sst<t))];
    Ref_conc_x = [R1.Ref_conc(:,find(R1.is_sst<t)) R2.Ref_conc(:,find(R2.is_sst<t)) R3.Ref_conc(:,find(R3.is_sst<t))];
    
    V = Ref_flux_x;
    clear x_min x_max
    for i=1:length(SCC_i)
        if all(~isnan(B_min{i}))
            u = unique(ODE{SCC_i(i)});
            o1=[];u1=[];
            for j=1:length(u)
                o1(j) = min(B_min{SCC_i(i)}(ODE{SCC_i(i)}==u(j)) .* (V(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),1),:)./V(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),2),:)))';
                u1(j) = max(B_max{SCC_i(i)}(ODE{SCC_i(i)}==u(j)) .* (V(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),1),:)./V(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),2),:)))';
            end
            x_min(i,1) = max(o1); x_max(i,1) = min(u1);
        else
            x_min(i,1) = NaN; x_max(i,1) = NaN;
        end
    end
    
    x_min_temp = x_min;
    x_min((round((x_min./x_max)*1e3)/1e3)>1) = NaN;
    x_max((round((x_min_temp./x_max)*1e3)/1e3)>1) = NaN;
    
    c_lb = min(Ref_conc_x')'; c_ub = max(Ref_conc_x')';
    
    x_range=abs(x_max-x_min);
    
    % compare ranges using correlation and euclidean distance
    corval(run,1)=corr(x_range,Simulated_range(SCC_i),'rows','pairwise');
    corval_lb(run,1)=corr(x_min,c_lb(SCC_i),'rows','pairwise');
    corval_lb_s(run,1)=corr(x_min,c_lb(SCC_i),'type','Spearman','rows','pairwise');
    corval_ub(run,1)=corr(x_max,c_ub(SCC_i),'rows','pairwise');
    corval_ub_s(run,1)=corr(x_max,c_ub(SCC_i),'type','Spearman','rows','pairwise');
    
    distvall(run,1)=pdist([x_min(~isnan(x_min))';c_lb(SCC_i(~isnan(x_min)))'])/length(SCC_i(~isnan(x_min)));
    distvall_m(run,1)=pdist([x_min(~isnan(x_min))'./(max(x_min(~isnan(x_min))));c_lb(SCC_i(~isnan(x_min)))'./max(c_lb(SCC_i(~isnan(x_min))))])/length(SCC_i(~isnan(x_min)));
    distvall_l(run,1)=pdist([log(x_min(~isnan(x_min))');log(c_lb(SCC_i(~isnan(x_min))))'])/length(SCC_i(~isnan(x_min)));
    
    distvalu(run,1)=pdist([x_max(~isnan(x_max))';c_ub(SCC_i(~isnan(x_max)))'])/length(SCC_i(~isnan(x_max)));
    distvalu_m(run,1)=pdist([x_max(~isnan(x_max))'./(max(x_max(~isnan(x_max))));c_ub(SCC_i(~isnan(x_max)))'./max(c_ub(SCC_i(~isnan(x_max))))])/length(SCC_i(~isnan(x_max)));
    distvalu_l(run,1)=pdist([log(x_max(~isnan(x_max))');log(c_ub(SCC_i(~isnan(x_max))))'])/length(SCC_i(~isnan(x_max)));
    
end
save(savename)
end
