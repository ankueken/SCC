addpath('Marans_Ecoli_kinetic_core_model/Users/')
load('Marans_Ecoli_kinetic_core_model/Users/Data.mat')

%% WT reference
% we use initial conditions provided by the model to simulate a steady-state
% flux distribution
[~,C,V]=solve_ode(Network_Data,0:100,Metabolite_and_Enzyme_Conc,Elementary_Kinetic_Param);

%% predicted fold change
load('Results_SCC_conc_upon_rxn_ko_fixed_EX.mat','x_min','x_min_del','SCC_i','s_del')
x_min_del(isnan(x_min_del))=0;x_min_del(x_min_del==Inf)=1e10;
PredFoldChange = (x_min_del-repmat(C(SCC_i,end),1,1474))./repmat(C(SCC_i,end),1,1474);
Flux_diff_pred = abs(s_del - repmat(V(:,end),1,1474));

%% simulated fold change
load('SimFoldChanges_All.mat','CmA','SST','V_ko')
SST(min(V_ko)<-1e-8)=0;
V_ko(:,SST==0)=0;

SimFoldChange = (CmA(:,SCC_i)-repmat(C(SCC_i,end)',1474,1))./repmat(C(SCC_i,end)',1474,1);
SimFoldChange = SimFoldChange';
Flux_diff_sim = abs(V_ko - repmat(V(:,end),1,1474));

%%
EDGES=[-Inf -1e3 -1e1 -1 -1e-1 -1e-3 0 1e-3 1e-1 1 1e1 1e3 Inf];
H=histc(SimFoldChange(:),EDGES);
H2=histc(PredFoldChange(:),EDGES);
bar([H(1:end-1) H2(1:end-1)]/length(SimFoldChange(:)),1)
set(gca,'XTickLabel',{'(Inf, -10^{3})' '[-10^3, -10^1)' '[-10^1, -1)' '[-1, -10^{-1})' '[-10^{-1}, -10^{-3})' '[-10^{-3}, 0)' '[0, 10^{-3})' '[10^{-3}, 10^{-1})' '[10^{-1}, 1)' '[1, 10^1)' '[10^1, 10^3)' '[10^3, Inf)'})
set(gca,'XTickLabelRotation',45)
xlabel('Fold change')
ylabel('Fraction')
legend('Simulated fold change','Predicted fold change','Orientation','horizontal','Location','North')
ylim([0 1])
colormap('hot')

for i=1:12
    EDGES=[-Inf -1e3 -1e1 -1 -1e-1 -1e-3 0 1e-3 1e-1 1 1e1 1e3 Inf];
    subplot(3,4,i)
    H=histc(SimFoldChange(i,:),EDGES);
    H2=histc(PredFoldChange(i,:),EDGES);
    bar([H(1:end-1);H2(1:end-1)]'./length(SimFoldChange(i,:)),1,'grouped')
    set(gca,'XTickLabel',[])
    ylim([0 1])
    title(Network_Data.metab(SCC_i(i)))
end
for i=13:23
    EDGES=[-Inf -1e3 -1e1 -1 -1e-1 -1e-3 0 1e-3 1e-1 1 1e1 1e3 Inf];
    subplot(3,4,i-12)
    H=histc(SimFoldChange(i,:),EDGES);
    H2=histc(PredFoldChange(i,:),EDGES);
    bar([H(1:end-1);H2(1:end-1)]'./length(SimFoldChange(i,:)),1,'grouped')
    set(gca,'XTickLabel',[])
    ylim([0 1])
    title(Network_Data.metab(SCC_i(i)))
end

plot(Flux_diff_sim(:),Flux_diff_pred(:),'.','MarkerSize',12)
hold on 
plot([-500000 .35e8],[-500000 .35e8],'k')
xlim([-500000 .35e8])
ylim([-500000 .35e8])
xlabel({'Difference flux distribution';'simulated knock-out to wild type'})
ylabel({'Difference flux distribution'; 'MOMA knock-out to wild type'})
colormap([0.5 0.5 0.5])

sum(SST)/length(SST)
P=PredFoldChange;
PA=PredFoldChange_adjusted;
S=SimFoldChange;

CORR_FC_ALL(1,1) = corr(P(:),S(:),'rows','pairwise');
CORR_FC_ALL(1,2) = corr(P(:),S(:),'type','Spearman','rows','pairwise');
CORR_FC_ALL(1,3) = corr(PA(:),S(:),'rows','pairwise');
CORR_FC_ALL(1,4) = corr(PA(:),S(:),'type','Spearman','rows','pairwise')
    
    