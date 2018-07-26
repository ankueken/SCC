changeCobraSolver('tomlab_cplex')

%% read experimental data
Data_c = readtable('Ishii_Quantitative_data_E_coli_metabolites.csv','Delimiter',',','HeaderLines',2);
Data_c{:,3:end} = Data_c{:,3:end}.*0.0023; % from mmol/L to mmol/gDW
SR = readtable('Ishii_Quantitative_data_E_coli_specific_rates.csv','Delimiter',',','HeaderLines',1);

%% load results of cytosolic SCC metabolites
SCC_naive = load('E_coli_iJO1366_cytosolic_SCC.mat');
SCC_i = SCC_naive.SCC_if(Data_c.ModelID); % index of SCC metabolites comprised in experimental dataset
SCC_if = SCC_naive.SCC_if; % index of all SCC metabolites
ReactionSet = SCC_naive.ReactionSet; ODE = SCC_naive.ODE;
Model_o = load('Results/Result_E_coli_K12_iJO1366.mat','Model_pro','SCC');
relevant_rxns = unique(cell2mat(ReactionSet));
epsilon=0;epsObj=1e-8;

%% specify objective
Model_o.Model_pro.ub(:) = 1000;
Model_o.Model_pro.lb(:) = 1e-8;
Model_o.Model_pro.csense(1:length(Model_o.Model_pro.mets),1) = 'E';
Model_o.Model_pro.c(:) = -1; % minimize total flux
Model_o.Model_pro.c(1110) = 0.1; % maximize ATP synthesis

%% learn B from RF 1:3
for sampleB=1:3
    Model_pro = Model_o.Model_pro;
    
    Model_pro.lb(7) = 0.2;Model_pro.ub(7) = 0.2; % fix biomass according to measurement
    
    % constrain exchange rates according to measurements
    Export = find(cellfun(@isempty,strfind(Model_pro.rxns,'EX_'))==0);
    Model_pro.lb(setdiff(Export,relevant_rxns)) = 0;
    Model_pro.lb(find(strcmp(Model_pro.rxns,'EX_glc(e)_rev'))) = SR{1,6+sampleB}-(SR{1,6+sampleB}*epsilon);
    Model_pro.lb(find(strcmp(Model_pro.rxns,'EX_o2(e)_rev'))) = SR{2,6+sampleB}-(SR{2,6+sampleB}*epsilon);
    Model_pro.lb(find(strcmp(Model_pro.rxns,'EX_co2(e)'))) = SR{3,6+sampleB}-(SR{3,6+sampleB}*epsilon);
    Model_pro.ub(find(strcmp(Model_pro.rxns,'EX_glc(e)_rev'))) = SR{1,6+sampleB}+(SR{1,6+sampleB}*epsilon);
    Model_pro.ub(find(strcmp(Model_pro.rxns,'EX_o2(e)_rev'))) = SR{2,6+sampleB}+(SR{2,6+sampleB}*epsilon);
    Model_pro.ub(find(strcmp(Model_pro.rxns,'EX_co2(e)'))) = SR{3,6+sampleB}+(SR{3,6+sampleB}*epsilon);
    
    Sol = optimizeCbModel(Model_pro);
    f = Sol.f;
    Vc(:,sampleB) = Sol.x;
    
    % flux range for relevant reactions
    T=evalc('[vmin,vmax]=fluxVariability(Model_pro,100,''max'',Model_pro.rxns(unique(cell2mat(ReactionSet))))');
    bs = ones(size(Model_pro.rxns)); % one if min < 1
    bs(find(vmin>1)) = -1; % otherwise -1 (needed in fractional LP)
    vp_vs_ratio_B=[];
    
    % get min/max ratio of relevant reaction pairs (fractional LP)
    for i=1:length(SCC_i)
        for j=1:size(ReactionSet{SCC_i(i)},1)
            Model_pro_FVA = Model_pro;
            
            y_i = zeros(1,size(Model_pro.S,2));
            y_i(ReactionSet{SCC_i(i)}(j,2)) = 1;
            
            Model_pro_FVA.S = [Model_pro.S zeros(size(Model_pro.S,1),1);  % steady state
                Model_o.Model_pro.c' 0; % at optimal solution
                y_i bs(find(y_i==1))  % y_j + bt = 1, b=0
                ];
            
            Model_pro_FVA.csense = [Model_pro.csense; 'G'; 'E'];
            
            Model_pro_FVA.lb = [Model_pro.lb; 1e-8];
            Model_pro_FVA.ub = [Model_pro.ub; 10000];
            Model_pro_FVA.b = [Model_pro.b; f+(f*epsObj); 1];
            Model_pro_FVA.c = zeros(size(Model_pro.S,2)+1,1);
            Model_pro_FVA.c(ReactionSet{SCC_i(i)}(j,1)) = 1;
            
            Sol_r = optimizeCbModel(Model_pro_FVA,'min');
            y=Sol_r.x;
            if Sol_r.stat~=1 % check with different solver if this is due to numerical issues
                [y, ~, ~,~, ~, ~, ~, Stat] = cplex(Model_pro_FVA.c, Model_pro_FVA.S, Model_pro_FVA.lb, Model_pro_FVA.ub, Model_pro_FVA.b, Model_pro_FVA.b);
                if Stat~=1
                    vp_vs_ratio_B{i}(j,1) = round((Vc(ReactionSet{SCC_i(i)}(j,1),sampleB)/Vc(ReactionSet{SCC_i(i)}(j,2),sampleB))*1e9)/1e9;
                else
                    x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e9)/1e9;
                    vp_vs_ratio_B{i}(j,1) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
                end
            else
                x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e9)/1e9;
                vp_vs_ratio_B{i}(j,1) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
            end
            
            Sol_r = optimizeCbModel(Model_pro_FVA,'max');
            y=Sol_r.x;
            if Sol_r.stat~=1
                [y, ~, ~,~, ~, ~, ~, Stat] = cplex(-Model_pro_FVA.c, Model_pro_FVA.S, Model_pro_FVA.lb, Model_pro_FVA.ub, Model_pro_FVA.b, Model_pro_FVA.b);
                
                if Stat~=1
                    vp_vs_ratio_B{i}(j,2) = round((Vc(ReactionSet{SCC_i(i)}(j,1),sampleB)/Vc(ReactionSet{SCC_i(i)}(j,2),sampleB))*1e9)/1e9;
                else
                    x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e9)/1e9;
                    vp_vs_ratio_B{i}(j,2) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
                end
            else
                x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e9)/1e9;
                vp_vs_ratio_B{i}(j,2) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
            end
        end
    end
    
    % get constants from each sample using measured concentration and
    % min/max ratio of relevant reaction pairs
    if sampleB==1
        V=Vc(:,1);
        B_c1_max=[];B_c1_min=[];
        for i=1:length(SCC_i)
            u = unique(ODE{SCC_i(i)});
            o1=[];u1=[];
            for j=1:length(u)
                temp_min = repmat(Data_c{i,8},size(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),:),1),1)./(vp_vs_ratio_B{i}((ODE{SCC_i(i)}==u(j)),2));
                temp_max = repmat(Data_c{i,8},size(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),:),1),1)./(vp_vs_ratio_B{i}((ODE{SCC_i(i)}==u(j)),1));
                B_c1_min{i}(find(ODE{SCC_i(i)}==u(j)),:) = temp_min;B_c1_max{i}(find(ODE{SCC_i(i)}==u(j)),:) = temp_max;
            end
        end
    elseif sampleB==2
        V=Vc(:,2);
        B_c2_max=[];B_c2_min=[];
        for i=1:length(SCC_i)
            u = unique(ODE{SCC_i(i)});
            o1=[];u1=[];
            for j=1:length(u)
                temp_min = repmat(Data_c{i,9},size(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),:),1),1)./(vp_vs_ratio_B{i}((ODE{SCC_i(i)}==u(j)),2));
                temp_max = repmat(Data_c{i,9},size(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),:),1),1)./(vp_vs_ratio_B{i}((ODE{SCC_i(i)}==u(j)),1));
                B_c2_min{i}(find(ODE{SCC_i(i)}==u(j)),:) = temp_min;B_c2_max{i}(find(ODE{SCC_i(i)}==u(j)),:) = temp_max;
            end
        end
    elseif sampleB==3
        V=Vc(:,3);
        B_c3_max=[];B_c3_min=[];
        for i=1:length(SCC_i)
            u = unique(ODE{SCC_i(i)});
            o1=[];u1=[];
            for j=1:length(u)
                temp_min = repmat(Data_c{i,10},size(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),:),1),1)./(vp_vs_ratio_B{i}((ODE{SCC_i(i)}==u(j)),2));
                temp_max = repmat(Data_c{i,10},size(ReactionSet{SCC_i(i)}((ODE{SCC_i(i)}==u(j)),:),1),1)./(vp_vs_ratio_B{i}((ODE{SCC_i(i)}==u(j)),1));
                B_c3_min{i}(find(ODE{SCC_i(i)}==u(j)),:) = temp_min;B_c3_max{i}(find(ODE{SCC_i(i)}==u(j)),:) = temp_max;
            end
        end
    end
end

%% concentration prediction for other GR using learned constants
GR = [0.1 0.4 0.5 0.7];

for sample = 2:4
    x_min_B1 = [];x_max_B1 = [];
    Model_pro = Model_o.Model_pro;
    
    Model_pro.lb(7) = GR(sample);Model_pro.ub(7) = GR(sample); % fix biomass to measurement
    
    % constrain exchange flux according to experimental data
    Export = find(cellfun(@isempty,strfind(Model_pro.rxns,'EX_'))==0);
    Model_pro.lb(setdiff(Export,relevant_rxns)) = 0;
    Model_pro.lb(find(strcmp(Model_pro.rxns,'EX_glc(e)_rev'))) = SR{1,1+sample}-(SR{1,1+sample}*epsilon);
    Model_pro.lb(find(strcmp(Model_pro.rxns,'EX_o2(e)_rev'))) = SR{2,1+sample}-(SR{2,1+sample}*epsilon);
    Model_pro.lb(find(strcmp(Model_pro.rxns,'EX_co2(e)'))) = SR{3,1+sample}-(SR{3,1+sample}*epsilon);
    Model_pro.ub(find(strcmp(Model_pro.rxns,'EX_glc(e)_rev'))) = SR{1,1+sample}+(SR{1,1+sample}*epsilon);
    Model_pro.ub(find(strcmp(Model_pro.rxns,'EX_o2(e)_rev'))) = SR{2,1+sample}+(SR{2,1+sample}*epsilon);
    Model_pro.ub(find(strcmp(Model_pro.rxns,'EX_co2(e)'))) = SR{3,1+sample}+(SR{3,1+sample}*epsilon);
    
    Sol = optimizeCbModel(Model_pro);
    Vg = Sol.x;
    f = Sol.f;
    
    % check for fractional LP if min flux of relevant reactions above 1
    T=evalc('[vmin,vmax]=fluxVariability(Model_pro,100,''max'',Model_pro.rxns(unique(cell2mat(ReactionSet))))');
    bs = ones(size(Model_pro.rxns));
    bs(find(vmin>1)) = -1;
    vp_vs_ratio=[];
    
    % min/max ratio of relevant reaction pairs
    for i=1:length(SCC_i)
        for j=1:size(ReactionSet{SCC_i(i)},1)
            Model_pro_FVA = Model_pro;
            
            y_i = zeros(1,size(Model_pro.S,2));
            y_i(ReactionSet{SCC_i(i)}(j,2)) = 1;
            
            Model_pro_FVA.S = [Model_pro.S zeros(size(Model_pro.S,1),1);  % steady state
                Model_o.Model_pro.c' 0; % at optimal solution
                y_i bs(find(y_i==1))  % y_j + bt = 1, b=0
                ];
            
            Model_pro_FVA.csense = [Model_pro.csense; 'G'; 'E'];
            
            Model_pro_FVA.lb = [Model_pro.lb; 1e-8];
            Model_pro_FVA.ub = [Model_pro.ub; 10000];
            Model_pro_FVA.b = [Model_pro.b; f+(f*epsObj); 1];
            Model_pro_FVA.c = zeros(size(Model_pro.S,2)+1,1);
            Model_pro_FVA.c(ReactionSet{SCC_i(i)}(j,1)) = 1;
            
            Sol_r = optimizeCbModel(Model_pro_FVA,'min');
            y=Sol_r.x;
            if Sol_r.stat~=1
                [y, ~, ~,~, ~, ~, ~, Stat] = cplex(Model_pro_FVA.c, Model_pro_FVA.S, Model_pro_FVA.lb, Model_pro_FVA.ub, Model_pro_FVA.b, Model_pro_FVA.b);
                if Stat~=1
                    vp_vs_ratio{i}(j,1) = round((Vg(ReactionSet{SCC_i(i)}(j,1))/Vg(ReactionSet{SCC_i(i)}(j,2)))*1e9)/1e9;
                else
                    x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e9)/1e9;
                    vp_vs_ratio{i}(j,1) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
                end
            else
                
                x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e9)/1e9;
                vp_vs_ratio{i}(j,1) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
            end
            
            Sol_r = optimizeCbModel(Model_pro_FVA,'max');
            y=Sol_r.x;
            if Sol_r.stat~=1
                [y, ~, ~,~, ~, ~, ~, Stat] = cplex(-Model_pro_FVA.c, Model_pro_FVA.S, Model_pro_FVA.lb, Model_pro_FVA.ub, Model_pro_FVA.b, Model_pro_FVA.b);
                if Stat~=1
                    vp_vs_ratio{i}(j,2) = round((Vg(ReactionSet{SCC_i(i)}(j,1))/Vg(ReactionSet{SCC_i(i)}(j,2)))*1e9)/1e9;
                else
                    x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e9)/1e9;
                    vp_vs_ratio{i}(j,2) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
                end
            else
                
                x = round(((1/y(end)).*y(1:size(Model_pro.S,2)))*1e9)/1e9;
                vp_vs_ratio{i}(j,2) = x(ReactionSet{SCC_i(i)}(j,1))/x(ReactionSet{SCC_i(i)}(j,2));
            end
        end
    end
    
    for sampleB=1:3 % using the learned constant from the three different samples
        if sampleB==1
            B_c_min=B_c1_min;B_c_max=B_c1_max;
        elseif sampleB==2
            B_c_min=B_c2_min;B_c_max=B_c2_max;
        elseif sampleB==3
            B_c_min=B_c3_min;B_c_max=B_c3_max;
        end
        
        % calculate concentration for SCC metabolites using the min/max
        % ratio of the relevant reaction pairs and the constants learned
        % from one out of the three experimental samples
        for i=1:length(SCC_i)
            if all(~isnan(B_c_min{i}))
                u = unique(ODE{SCC_i(i)});
                o1=[];u1=[];
                for j=1:length(u)
                    o1(j) = min(B_c_min{i}(ODE{SCC_i(i)}==u(j)) .* vp_vs_ratio{i}(ODE{SCC_i(i)}==u(j),1))';
                    u1(j) = max(B_c_max{i}(ODE{SCC_i(i)}==u(j)) .* vp_vs_ratio{i}(ODE{SCC_i(i)}==u(j),2))';
                end
                x_min_B1(i,sampleB) = max(o1); x_max_B1(i,sampleB) = min(u1);
            else
                x_min_B1(i,sampleB) = NaN; x_max_B1(i,sampleB) = NaN;
            end
        end
    end
    
    x_min_B1_temp = x_min_B1;
    x_min_B1((round((x_min_B1./x_max_B1)*1e3)/1e3)>1) = NaN;
    x_max_B1((round((x_min_B1_temp./x_max_B1)*1e3)/1e3)>1) = NaN;
    
    %% plot result
    if sample==2
        x_min_B1_new_s2 = x_min_B1; x_max_B1_new_s2 = x_max_B1;
        
        subplot(1,3,1)
        plot(1:length(SCC_i),Data_c{:,2+sample},'rx')
        
        x_min2 = min(x_min_B1')';x_max2 = max(x_max_B1')';
        
        for i=1:length(SCC_i)
            hold on
            plot([i i],[x_min2(i,1) x_max2(i,1)],'k')
            hold on
            plot([i-0.2 i+0.2],[x_min2(i,1) x_min2(i,1)],'k')
            hold on
            plot([i-0.2 i+0.2],[x_max2(i,1) x_max2(i,1)],'k')
        end
        for i=1:length(SCC_i)
            XL(i,1) = strcat(num2str(i),'-',Model_pro.mets(SCC_i(i)));
        end
        set(gca,'XTick',1:length(SCC_i),'XTickLabel',Model_pro.mets(SCC_i),'XTickLabelRotation',45,'YScale','log')
        xlim([0 length(SCC_i)+1])
        title('GR 0.4')
        
        colorstr = ['b','g','m'];
        for i=1:length(SCC_i)
            for j=1:3
                hold on
                plot([i+(j*0.1+0.2);i+(j*0.1+0.2)],[x_min_B1(i,j); x_max_B1(i,j)],'Color',colorstr(j))
                hold on
                plot([i+(j*0.1+0.15);i+(j*0.1+0.25)],[x_min_B1(i,j); x_min_B1(i,j)],'Color',colorstr(j))
                hold on
                plot([i+(j*0.1+0.15);i+(j*0.1+0.25)],[x_max_B1(i,j); x_max_B1(i,j)],'Color',colorstr(j))
            end
        end
        
    elseif sample==3
        x_min_B1_new_s3 = x_min_B1; x_max_B1_new_s3 = x_max_B1;
        
        subplot(1,3,2)
        plot(1:length(SCC_i),Data_c{:,2+sample},'rx')
        
        x_min2 = min(x_min_B1')';x_max2 = max(x_max_B1')';
        
        for i=1:length(SCC_i)
            hold on
            plot([i i],[x_min2(i,1) x_max2(i,1)],'k')
            hold on
            plot([i-0.2 i+0.2],[x_min2(i,1) x_min2(i,1)],'k')
            hold on
            plot([i-0.2 i+0.2],[x_max2(i,1) x_max2(i,1)],'k')
        end
        for i=1:length(SCC_i)
            XL(i,1) = strcat(num2str(i),'-',Model_pro.mets(SCC_i(i)));
        end
        set(gca,'XTick',1:length(SCC_i),'XTickLabel',Model_pro.mets(SCC_i),'XTickLabelRotation',45,'YScale','log')
        xlim([0 length(SCC_i)+1])
        ylabel('mmol/gDW')
        title('GR 0.5')
        
        colorstr = ['b','g','m'];
        for i=1:length(SCC_i)
            for j=1:3
                hold on
                plot([i+(j*0.1+0.2);i+(j*0.1+0.2)],[x_min_B1(i,j); x_max_B1(i,j)],'Color',colorstr(j))
                hold on
                plot([i+(j*0.1+0.15);i+(j*0.1+0.25)],[x_min_B1(i,j); x_min_B1(i,j)],'Color',colorstr(j))
                hold on
                plot([i+(j*0.1+0.15);i+(j*0.1+0.25)],[x_max_B1(i,j); x_max_B1(i,j)],'Color',colorstr(j))
            end
        end
        
    elseif sample==4
        x_min_B1_new_s4 = x_min_B1; x_max_B1_new_s4 = x_max_B1;
        
        subplot(1,3,3)
        plot(1:length(SCC_i),Data_c{:,2+sample},'rx')
        
        x_min2 = min(x_min_B1')';x_max2 = max(x_max_B1')';
        
        for i=1:length(SCC_i)
            hold on
            plot([i i],[x_min2(i,1) x_max2(i,1)],'k')
            hold on
            plot([i-0.2 i+0.2],[x_min2(i,1) x_min2(i,1)],'k')
            hold on
            plot([i-0.2 i+0.2],[x_max2(i,1) x_max2(i,1)],'k')
        end
        for i=1:length(SCC_i)
            XL(i,1) = strcat(num2str(i),'-',Model_pro.mets(SCC_i(i)));
        end
        set(gca,'XTick',1:length(SCC_i),'XTickLabel',Model_pro.mets(SCC_i),'XTickLabelRotation',45,'YScale','log')
        xlim([0 length(SCC_i)+1])
        title('GR 0.7')
        
        colorstr = ['b','g','m'];
        for i=1:length(SCC_i)
            for j=1:3
                hold on
                plot([i+(j*0.1+0.2);i+(j*0.1+0.2)],[x_min_B1(i,j); x_max_B1(i,j)],'Color',colorstr(j))
                hold on
                plot([i+(j*0.1+0.15);i+(j*0.1+0.25)],[x_min_B1(i,j); x_min_B1(i,j)],'Color',colorstr(j))
                hold on
                plot([i+(j*0.1+0.15);i+(j*0.1+0.25)],[x_max_B1(i,j); x_max_B1(i,j)],'Color',colorstr(j))
            end
        end
    end
    
end


save('Results_Ishii_max_ATP.mat')