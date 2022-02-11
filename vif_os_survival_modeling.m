clc
close all
clear all
home


%% Load the radiomics features
data=readtable('/home/data/radiomics_features.xlsx','ReadVariableNames', true);
data_features = data(:,7:end);
data_cli = data(:,1:6);
features_name=data_features.Properties.VariableNames; 

features = data_features{:,:};
features_org = features;
features = zscore(features);

% VIF calculation and use less than values 10
vif_features = vif_edit(features);
vif_high_idx = vif_features>=10;
vif_low_idx = ~vif_high_idx;
org_features_selected_name = features_name(:,vif_low_idx);
org_features_selected = features_org(:,vif_low_idx);
vif_features = array2table(org_features_selected,'VariableNames',org_features_selected_name);
vif_total = [data_cli, vif_features];

%% Load the data partitioning for 5-fold cross validation 
c=cvpartition(size(data,1),'k',5);
save('/home/cv.mat','cv')

% cv_298 = readtable('/home/cv.mat');
% for i = 1:5
%     cv(i).training = logical(table2array(cv_298(:,i)));
%     cv(i).test = ~cv(i).training;
% end

%% Load the VIF features
data =vif_total;
start_num=7;
features = data(:, start_num:end);
features_name= features.Properties.VariableNames; % 269 : 232 || 260: 223
features = features{:,:};

data =readtable('/home/vif_data/'); % radiomics features <= VIF 10
columns=data.Properties.VariableNames;

% features starts at index 7start_num=7; 
features = data(:, start_num:end);
features_name= features.Properties.VariableNames; 
features = features{:,:};

%% Clinical variables (age, sex, tumor size)
clinic_features = data(:, 4:6); 
clinic_features = clinic_features{:,:};

% SZM only
save_path = '/home/model/srm';
features = features(:,253:end);
features_name_vif = features_name(:,253:end);

% Radiomics only
save_path = '/home/model/rm';
features = features(:,1:228);
features_name_vif = features_name(:,1:228);

% Radiomics + SZM
save_path = '/home/model/erm';
features_name_vif = features_name;


%% Setting
% Lasso selection loading
load_mode=0;
load_mode=1;

% with Clinical parameters
clinical_mode = 0;
clinical_mode = 1; 

%% 5-fold cross validation modeling
result = struct();
cnt=0;
i = 1;
while i<=5
    disp(i)
    disp('fold start!')
    
    TrainingFeature = features(cv(i).training,:); 
    TestFeature =features(cv(i).test,:);

    % z-normalization
    TrainingFeature = (TrainingFeature-mean(TrainingFeature))./std(TrainingFeature);
    TestFeature = (TestFeature-mean(TrainingFeature))./std(TrainingFeature);
    
    % set for death=0, 1=censored
    TrainingLabel =[data.Osperiod(cv(i).training,1), ~data.OS(cv(i).training,1)];  
    TestLabel = [data.Osperiod(cv(i).test,1),~data.OS(cv(i).test,1)]; 
    
    % feature selection by LASSO for Cox
    cox_lasso = cvglmnet(TrainingFeature,TrainingLabel,'cox');
    selectedFeature_idx = find(cox_lasso.glmnet_fit.beta(:,(cox_lasso.lambda == cox_lasso.lambda_min)));
    result(i).selectedFeatureName = features_name_vif(:,selectedFeature_idx);
    result(i).selectedFeature_coef =cox_lasso.glmnet_fit.beta(selectedFeature_idx,(cox_lasso.lambda == cox_lasso.lambda_min));

    selectedFeature = TrainingFeature(:,selectedFeature_idx);
    selectedFeature_te = TestFeature(:,selectedFeature_idx);
    
    % With clinical variables
    if clinical_mode==1
        TrainingFeature_cli = clinic_features(cv(i).training,:);
        TestFeature_cli = clinic_features(cv(i).test,:);    
        
        TrainingFeature_cli = (TrainingFeature_cli-mean(TrainingFeature_cli))./std(TrainingFeature_cli);
        TestFeature_cli = (TestFeature_cli-mean(TrainingFeature_cli))./std(TrainingFeature_cli);

        selectedFeature = [TrainingFeature_cli, selectedFeature];
        selectedFeature_te = [TestFeature_cli, selectedFeature_te];
    end
        
    
    %% COX modeling
    [b,logl,H,stats] = coxphfit(selectedFeature, TrainingLabel(:,1),'Censoring',TrainingLabel(:,2));
    result(i).b_tr=b;
    result(i).rad_score = selectedFeature * b;
    result(i).tr_m_rad_pfs = median(result(i).rad_score);

    [b_te,logl,H,stats_te] = coxphfit(selectedFeature_te, TestLabel(:, 1),'Censoring',TestLabel(:,2)); 
    result(i).b_te = b_te;
    result(i).rad_score_te = selectedFeature_te *b_te;
    result(i).te_m_rad_pfs = median(result(i).rad_score_te);
    
    
    %% Evaluate and save the cox model
    [lambda1, lambda2, result(i).trHR, result(i).trHRci, UL, SUL, z, result(i).p_pfs_tr, alpha] = logrank([TrainingLabel(result(i).rad_score>=result(i).tr_m_rad_pfs,1),TrainingLabel(result(i).rad_score>=result(i).tr_m_rad_pfs,2)],[TrainingLabel(result(i).rad_score<result(i).tr_m_rad_pfs,1),TrainingLabel(result(i).rad_score<result(i).tr_m_rad_pfs,2)]);
    result(i).trHRci = round(result(i).trHRci,2);
    
    [lambda1, lambda2, result(i).teHR, result(i).teHRci, UL, SUL, z, result(i).p_pfs_te, alpha] = logrank([TestLabel(result(i).rad_score_te>=result(i).te_m_rad_pfs,1),TestLabel(result(i).rad_score_te>=result(i).te_m_rad_pfs,2)], [TestLabel(result(i).rad_score_te<result(i).te_m_rad_pfs,1), TestLabel(result(i).rad_score_te<result(i).te_m_rad_pfs,2)]);
    result(i).teHRci = round(result(i).teHRci,2);

    save(strcat(save_path,'\variable_coxlasso_',num2str(i),'-fold.mat'),'-struct','cox_lasso')    
    saveas(gcf, strcat(save_path,'\logrank_',num2str(i),'-fold_p_,',num2str(round(result(i).p_pfs_te,4)),'.png'));
    result_1=result(i);
    save(strcat(save_path,'\variable_result',num2str(i),'-fold.mat'),'-struct','result_1')

    
    i = i+1; 
end


