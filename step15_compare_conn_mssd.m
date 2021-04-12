% Compare whole connectivity and/or MSSD to classiciation ASD/CON
clc;
clear all;
close all;

%% Settings

% add path
addpath(pwd);

% set path
cd('../..');
Dirworking = pwd;
Dirlog = [pwd '/log'];
Dirdocu = [pwd '/documents'];
Dirdata = [pwd '/subs/data'];
Dirmask = [pwd '/mask'];
Diranalysis = [pwd '/analysis/UTHSC_data/data'];

% Set parameters
option.method = 'svm'; %'kernel_fda', 'svm'
% option.metric = {'auc', 'accuracy', 'dval', 'tval'};
option.metric = {'auc'};
option.preprocess = {'oversample'};
option.statmetric = 'accuracy';

option.windowsize = 40;         % current TR parameter
option.slidingwindow = 'robust';  % current sliding window type 'sw'. 'tsw'...
option.signaltype = 'Sig';  % Raw Dfc signals
option.metricindex = 4;     % 1 is Raw Dfc, 2 is mean, 3 is SD, 4 is MSSD1,...


option.allROI = 1:64; % bilateral AAL2 atals number
option.outROI = [48:64];

% option.allROI = 1:120;
% option.outROI = 95:120;
option.incROI = setdiff(option.allROI, option.outROI);

% Downsampling rate
TR_down = 3;

%% Get subjects Information
Script_BasicInfo;

% Divide the age group
Sub_agegroup1 = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 11.56;  % child
Sub_agegroup2 = BasicInfo.AgeAtScan > 11.56 & BasicInfo.AgeAtScan <= 14.44; % young adolescent
Sub_agegroup3 = BasicInfo.AgeAtScan > 14.44 & BasicInfo.AgeAtScan <= 19.71; % older adolescent
Sub_agegroup4 = BasicInfo.AgeAtScan > 19.71;                                % adult
Sub_agegroup1to2 = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 14.44; % child to young adolescent
Sub_agegroup1to3 = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 19.71; % child to older adolescent
Sub_agegroup2to3 = BasicInfo.AgeAtScan > 11.56 & BasicInfo.AgeAtScan <= 19.71; % child to older adolescent


% Sub_agegroup1 = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 12;  % child
% % Sub_agegroup1 = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 11.56 | BasicInfo.AgeAtScan >= 12 & BasicInfo.AgeAtScan <= 14.44;
% Sub_agegroup2 = BasicInfo.AgeAtScan > 12 & BasicInfo.AgeAtScan <= 20; % adolescent
% Sub_agegroup3 = BasicInfo.AgeAtScan > 20 & BasicInfo.AgeAtScan <= 45; % adult
% Sub_agegroup4 = BasicInfo.AgeAtScan > 45;                                % older
% Sub_agegroup1to2 = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 20; % child to adolescent
% Sub_agegroup1to3 = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 45; % child to older
% Sub_agegroup2to3 = BasicInfo.AgeAtScan > 12 & BasicInfo.AgeAtScan <= 45; % child to older

Sub_agegroup = Sub_agegroup1 + Sub_agegroup2.*2 + Sub_agegroup3.*3 + Sub_agegroup4.*4;
fprintf('Child: %d / Younger Adolescent: %d / Older Adolescent: %d / Adult: %d \n',sum(Sub_agegroup1),sum(Sub_agegroup2),sum(Sub_agegroup3),sum(Sub_agegroup4));

Sub_Con = strcmp(BasicInfo.SubjectType,'CONTROL');
Sub_Pat = strcmp(BasicInfo.SubjectType,'PATIENT');
Sub_group = Sub_Con + Sub_Pat.*2;

BasicInfo = addvars(BasicInfo, Sub_agegroup, 'After','AgeAtScan');
BasicInfo = addvars(BasicInfo, Sub_group, 'After', 'SubjectType');
BasicInfo = addvars(BasicInfo, Sub_agegroup1to2, 'After', 'Sub_agegroup');
BasicInfo = addvars(BasicInfo, Sub_agegroup1to3, 'After', 'Sub_agegroup1to2');
BasicInfo = addvars(BasicInfo, Sub_agegroup2to3, 'After', 'Sub_agegroup1to3');

%% load motion parameter information
List_subj = BasicInfo.subID;

cd(Dirdata);
for nSubj = 1:length(List_subj)
   nSubj;
   cd(List_subj{nSubj});
   cd('raw/func_resting2/motionparam');
   TMP_motionFD = load('tsMotion_FD.txt');
   Idx_over025 = find(TMP_motionFD >= 0.2);
   Idx_portion = (length(Idx_over025)/length(TMP_motionFD))*100;
   FDportion(nSubj,1) = Idx_portion;
   
   cd(Dirdata);
   clear Idx_over025 Idx_portion
end

BasicInfo = addvars(BasicInfo, FDportion, 'After','FD');
Idx_con = find(BasicInfo.Sub_group==1 & BasicInfo.FDportion <= 20);
Idx_pat = find(BasicInfo.Sub_group==2 & BasicInfo.FDportion <= 20);
% Idx_con = find(BasicInfo.Sub_group==1 & BasicInfo.Sub_agegroup==1);
% Idx_pat = find(BasicInfo.Sub_group==2 & BasicInfo.Sub_agegroup==1);
% Idx_con = find(BasicInfo.Sub_group==1 & BasicInfo.FDportion <= 20 & BasicInfo.Sub_agegroup==4);
% Idx_pat = find(BasicInfo.Sub_group==2 & BasicInfo.FDportion <= 20 & BasicInfo.Sub_agegroup==4);
% Idx_con = find(BasicInfo.Sub_group==1 & BasicInfo.FDportion <= 25 & BasicInfo.Sub_agegroup1to2==1);
% Idx_pat = find(BasicInfo.Sub_group==2 & BasicInfo.FDportion <= 25 & BasicInfo.Sub_agegroup1to2==1);

List_subj_specific = List_subj(sort([Idx_con; Idx_pat]));

%% Load data
cd(Diranalysis);


for nSubj = 1:length(List_subj)
    nSubj;
    % load LAIRD time series
    cd(List_subj{nSubj});

    TMP_AAL2_signal = load('ts_aal2_bilateral_func_resting2_nosmoothed_denoised_mni_epinorm.txt');
%     TMP_AAL2_signal = load('ts_aal2_func_resting2_nosmoothed_denoised_mni_epinorm.txt');
    TMP_AAL2_signal = TMP_AAL2_signal(:,option.incROI);
    
    
    Flip_TMP_signal = [flip(TMP_AAL2_signal); TMP_AAL2_signal; flip(TMP_AAL2_signal)];
    
    % Downsampling to TR = 3sec
    TR_orig = BasicInfo.TRSec(nSubj);
    L = ceil(size(TMP_AAL2_signal,1)*(floor(TR_orig*1000))/(floor(TR_down*1000)));
    TMP_AAL2_signal = resample(Flip_TMP_signal,floor(TR_orig*1000),floor(TR_down*1000));
    AAL2_signal = TMP_AAL2_signal([L:2*L-1],:);
    
    Sig_allsubj{nSubj,1} = AAL2_signal;
        
    clear AAL2_signal TMP_AAL2_signal
    cd(Diranalysis);
end

%% Calculate Dynamic connectivty using sliding window (varying window size)
windowsize = option.windowsize;

for nSubj = 1:length(List_subj)  % When you use 'DCC' change to for instaed of parfor

nSubj;

if strcmp(option.slidingwindow,'sw')
    [TMP_Dfc] = Script_Dfc(Sig_allsubj{nSubj,1},windowsize,1);
elseif strcmp(option.slidingwindow,'tsw')
    [TMP_Dfc] = Script_Dfc(Sig_allsubj{nSubj,1},windowsize,2);
elseif strcmp(option.slidingwindow,'DCC')
    [TMP_Dfc] = Script_Dfc(Sig_allsubj{nSubj,1},windowsize,3);
elseif strcmp(option.slidingwindow,'robust')
    [TMP_Dfc] = Script_Dfc(Sig_allsubj{nSubj,1},windowsize,4);
end
    
[Mat_MEAN Mat_SD Mat_MSSD1 Mat_MSSD2 Mat_VSD] = Script_calMSSD(TMP_Dfc);

Mat_MEAN = tril(Mat_MEAN,-1);
Mat_SD = tril(Mat_SD,-1);
Mat_MSSD1 = tril(Mat_MSSD1,-1);
Mat_MSSD2 = tril(Mat_MSSD2,-1);
Mat_VSD = tril(Mat_VSD,-1);

Dfc_metrics1{nSubj,1} = TMP_Dfc;
Dfc_metrics2{nSubj,1} = Mat_MEAN(Mat_MEAN~=0);
Dfc_metrics3{nSubj,1} = Mat_SD(Mat_SD~=0);
Dfc_metrics4{nSubj,1} = Mat_MSSD1(Mat_MSSD1~=0);
Dfc_metrics5{nSubj,1} = Mat_MSSD2(Mat_MSSD2~=0);
Dfc_metrics6{nSubj,1} = Mat_VSD(Mat_VSD~=0);
           
Mat_MEAN=[]; Mat_SD=[]; Mat_MSSD1=[]; Mat_MSSD2=[]; Mat_VSD=[];
TMP_Dfc = [];TMP_Dfc_tsw=[];
     
end

% Summarize to one variable
for nIdx = 1:6
    eval(['Dfc_metrics(:,' num2str(nIdx) ') = Dfc_metrics' num2str(nIdx) '(:,1);']);
end

clear Dfc_metrics1 Dfc_metrics2 Dfc_metrics3 Dfc_metrics4 Dfc_metrics5 Dfc_metrics6

%% Calculating static connectivity strength
if strcmp(option.slidingwindow,'sw') || strcmp(option.slidingwindow,'tsw')
    parfor nSubj = 1:length(List_subj)
        TMP_sig = Sig_allsubj{nSubj,1};
        TMP_corr = tril( atanh(corr(TMP_sig)),-1 );
        fc_static{nSubj,1} = TMP_corr(TMP_corr~=0);
        TMP_sig = [];
    end
elseif strcmp(option.slidingwindow,'DCC')
    for nSubj = 1:length(List_subj)
        TMP_corr = mean(Dfc_metrics{nSubj,1},3);
        TMP_corr = tril(TMP_corr,-1);
        fc_static{nSubj,1} = TMP_corr(TMP_corr~=0);
    end    
elseif strcmp(option.slidingwindow,'robust')
    parfor nSubj = 1:length(List_subj)
        TMP_sig = Sig_allsubj{nSubj,1};
        [T,p] = size(TMP_sig);
        TMP_corr = NaN(p,p);

        for j = 1:size(TMP_sig,2)
            for k = 1:size(TMP_sig,2)
                data1 = TMP_sig(:,j); data2 = TMP_sig(:,k);
                [~,~,rob_corrw] = andlab_robustfit( data1, data2);
                TMP_corr(j,k) = atanh(rob_corrw);
            end
        end
        
        TMP_corr = tril(TMP_corr,-1);
        fc_static{nSubj,1} = TMP_corr(TMP_corr~=0);
        TMP_sig = [];        
    end        
    save fc_static_bilateral_robust fc_static;
end

%% Reshape data matrix

parfor nSubj = 1:size(fc_static,1)
%     Data_subj = Data{nSubj,option.metricindex}(:,:,1:min_size);
    Data_fc = fc_static{nSubj};
    Mat_fc_concat(nSubj,:) = zscore(Data_fc);
    
    Data_metric = Dfc_metrics{nSubj,option.metricindex};
    Mat_Dfc_concat(nSubj,:) = zscore(Data_metric);
        
    Mat_combine(nSubj,:) = cat(1,zscore(Data_fc),zscore(Data_metric));
    Data_fc = []; Data_metric = [];
end
fprintf('Reshape is finished \n');

Mat_con = Mat_fc_concat(Idx_con,:);
Mat_pat = Mat_fc_concat(Idx_pat,:);
Mat_fc_whole = cat(1,Mat_con,Mat_pat);

Mat_con = Mat_Dfc_concat(Idx_con,:);
Mat_pat = Mat_Dfc_concat(Idx_pat,:);
Mat_Dfc_whole = cat(1,Mat_con,Mat_pat);

Mat_con = Mat_combine(Idx_con,:);
Mat_pat = Mat_combine(Idx_pat,:);
Mat_combine_whole = cat(1,Mat_con,Mat_pat);

%% Dimension Reduction
% pparam = [];
% pparam.n = nPCA;
% pparam.is_train_set = 1;
% pparam.target_dimension = 1;
% pparam.feature_dimension = 2;
% pparam.normalize = 1;
% pparam.design = [ones(size(Mat_con,1),1); 2*ones(size(Mat_pat,1),1)];
% [pparam, PCA_fc_whole, clabel] = mv_preprocess_pca(pparam, Mat_fc_whole, pparam.design);
% [pparam, PCA_Dfc_whole, clabel] = mv_preprocess_pca(pparam, Mat_Dfc_whole, pparam.design);
% [pparam, PCA_combine_whole, clabel] = mv_preprocess_pca(pparam, Mat_combine_whole, pparam.design);

% -----------------------------------------
% [PCA_coeff, PCA_score, PCA_latent, ~, PCA_explain] = pca(Mat_fc_whole(nSubj,:)');
% Recon_fc_concat = PCA_score*PCA_coeff';
% PCA_fc_whole(nSubj,:) = Recon_fc_concat';
% 
% [PCA_coeff, PCA_score, PCA_latent, ~, PCA_explain] = pca(Mat_Dfc_whole(nSubj,:)');
% Recon_Dfc_concat = PCA_score*PCA_coeff';
% PCA_Dfc_whole(nSubj,:) = Recon_Dfc_concat';
% 
% [PCA_coeff, PCA_score, PCA_latent, ~, PCA_explain] = pca(Mat_combine_whole(nSub,:)');
% Recon_combine = PCA_score*PCA_coeff';
% Mat_combine_whole(nSubj,:) = Recon_combine;

% ----------------------------------------- LASSO
% Idx_predictor = cell(1,size(Mat_fc_whole,2));
% Idx_subject = [zeros(size(Mat_con,1),1); ones(size(Mat_pat,1),1)];
% Idx_tmp = ones(size(Idx_subject,1),1); %1:size(Idx_subject,1);
% for nIdx = 1:length(Idx_predictor), Idx_predictor{nIdx} = num2str(nIdx); end
% [B,FitInfo] = lassoglm(Mat_fc_whole,Idx_subject,'binomial','CV',3); %,'PredictorNames',Idx_predictor);
% idxLambdaMinMSE = FitInfo.IndexMinMSE;
% minMSEModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinMSE)~=0);
% idxLambda1SE = FitInfo.Index1SE;
% sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0);
% result = cellfun(@str2num, sparseModelPredictors(:,1:end));

if exist('B_fc.mat') > 1 && exist('B_Dfc.mat') > 1 && exist('B_combine.mat') > 1
    load B_fc; load B_Dfc; load B_combine;
    load FitInfo_fc; load FitInfo_Dfc; load FitInfo_combine;
else
%     Idx_subject = BasicInfo.Sub_group;
%     Idx_subject(Idx_subject==2) = 0;
    Idx_subject = cat(1,zeros(length(Idx_con),1), ones(length(Idx_pat),1));
    
    % Lasso regression on whole data (all subject not limited to age group)
    [B_fc,FitInfo_fc] = lassoglm(Mat_fc_whole,Idx_subject,'binomial','CV',5);
    [B_Dfc,FitInfo_Dfc] = lassoglm(Mat_Dfc_whole,Idx_subject,'binomial','CV',5);
    [B_combine,FitInfo_combine] = lassoglm(Mat_combine_whole,Idx_subject,'binomial','CV',5);

    save B_fc B_fc; save FitInfo_fc FitInfo_fc;
    save B_Dfc B_Dfc; save FitInfo_Dfc FitInfo_Dfc;
    save B_combine B_combine; save FitInfo_combine FitInfo_combine;
end

% Predictor_fc = find(B_fc(:,FitInfo_fc.IndexMinDeviance));
% Predictor_Dfc = find(B_Dfc(:,FitInfo_Dfc.IndexMinDeviance));
% Predictor_combine = find(B_combine(:,FitInfo_combine.IndexMinDeviance));
n_feature = 100;
for ifeature = 1:n_feature

    [~,Predictor_fc] = maxk(abs(B_fc(:,1)), ifeature);
    [~,Predictor_Dfc] = maxk(abs(B_Dfc(:,1)), ifeature);
    [~,Predictor_combine] = maxk(abs(B_combine(:,1)), ifeature);

    % Predictor_fc = find(B_fc(:,1)~=0);
    % Predictor_Dfc = find(B_Dfc(:,1)~=0);
    % Predictor_combine = find(B_combine(:,1)~=0);

    % Predictor_combine = [Predictor_fc; Predictor_Dfc+size(B_fc,1)];

    % Make sparse matrices
    i_Mat_fc_whole = Mat_fc_whole(:,Predictor_fc);
    i_Mat_Dfc_whole = Mat_Dfc_whole(:,Predictor_Dfc);
    i_Mat_combine_whole = Mat_combine_whole(:,Predictor_combine);

    %% MVPA (SVM)
    cfg = [];
    % cfg.hyperparameter = param;
    cfg.classifier = option.method;
    cfg.preprocess = option.preprocess;
    %cfg.preprocess_param = {};
    %cfg.preprocess_param{2} = {'n', nPCA};
    cfg.metric     = option.metric;
    cfg.k          = 5;
    cfg.repeat     = 15;
    cfg.sample_dimension = 1; % dimension(subject)
    cfg.feature_dimension = 2; % dimension(feature)
    cfg.dimension_names = {'Subjects','metric'};
    cfg.design     = [ones(size(Mat_con,1),1); 2*ones(size(Mat_pat,1),1)];

    % --------------------------
    % % classification by functional connectivity (whole length)
    % [perf_fc, mvpa_fc] = mv_classify(cfg, PCA_fc_whole, cfg.design);
    % Cell_perf_fc{nPCA} = perf_fc; Cell_mvpa_fc{nPCA} = mvpa_fc;
    % % classification by MSSD
    % [perf_Dfc, mvpa_Dfc] = mv_classify(cfg, PCA_Dfc_whole, cfg.design);
    % Cell_perf_Dfc{nPCA} = perf_Dfc; Cell_mvpa_Dfc{nPCA} = mvpa_Dfc;
    % % classification by combine (functional connectivity & MSSD)
    % [perf_combine, mvpa_combine] = mv_classify(cfg, PCA_combine_whole, cfg.design);
    % Cell_perf_combine{nPCA} = perf_combine; Cell_mvpa_combine{nPCA} = mvpa_combine;
    % --------------------------
    % 
    [perf_fc, mvpa_fc] = mv_classify(cfg, i_Mat_fc_whole, cfg.design);
    Cell_perf_fc{ifeature} = perf_fc; Cell_mvpa_fc{ifeature} = mvpa_fc;
    [perf_Dfc, mvpa_Dfc] = mv_classify(cfg, i_Mat_Dfc_whole, cfg.design);
    Cell_perf_Dfc{ifeature} = perf_Dfc; Cell_mvpa_Dfc{ifeature} = mvpa_Dfc;
    [perf_combine, mvpa_combine] = mv_classify(cfg, i_Mat_combine_whole, cfg.design);
    Cell_perf_combine{ifeature} = perf_combine; Cell_mvpa_combine{ifeature} = mvpa_combine;

%     %Select AUC
%     if length(option.metric) >= 2
%         mvpa_fc_stat{nPCA} = mv_select_result(Cell_mvpa_fc{nPCA}, option.statmetric);
%         mvpa_Dfc_stat{nPCA} = mv_select_result(Cell_mvpa_Dfc{nPCA}, option.statmetric);
%         mvpa_combine_stat{nPCA} = mv_select_result(Cell_mvpa_combine{nPCA}, option.statmetric);
%     end

    clear i_Mat_fc_whole i_Mat_Dfc_whole i_Mat_combine_whole
end

Mat_perf_fc = cell2mat(Cell_perf_fc);
Mat_perf_Dfc = cell2mat(Cell_perf_Dfc);
Mat_perf_combine = cell2mat(Cell_perf_combine);
figure;
plot([1:ifeature],Mat_perf_fc); hold on
plot([1:ifeature],Mat_perf_Dfc); hold on
plot([1:ifeature],Mat_perf_combine);
legend('Static','MSSD','Combined');


% % figure; bar([1 2 3], [mvpa_fc_stat.perf mvpa_Dfc_stat.perf mvpa_combine_stat.perf]);
% figure; h = bar([perf_fc, perf_Dfc, perf_combine]);
% h.FaceColor = 'flat'; 
% h.CData = jet(numel([perf_fc, perf_Dfc, perf_combine]));
% set(gca,'XTick', 1:3, 'XTickLabels', {'Static','MSSD','Combine'})
% set(gca,'Ylim',[0.4 1]);
