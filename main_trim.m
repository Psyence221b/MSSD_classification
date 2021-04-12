% Variability of connectivity between mPFC and BLA (basolateral Amyg)
clear all
close all
clc

% add path
addpath(pwd);

% set path
cd('../..');
Dirworking = pwd;
Dirlog = [pwd '/log'];
Dirdocu = [pwd '/documents'];
Dirdata = [pwd '/analysis/UTHSC_data/data'];

%% Get subjects Information
Script_BasicInfo;

%% Get Data from 'analysis' folder
cd(Dirdata);

for nSubj = 1:size(BasicInfo,1)
    nSubj
    Fol_subj = BasicInfo.subID(nSubj);
    cd(Fol_subj{1});
%     File_BLA = dir('*_BLA_*'); File_BLA = File_BLA.name;
    File_mPFC = dir('*_find_mpfc_epi_norm_*'); File_mPFC = File_mPFC.name;
    File_BLA = dir('*_BLA_epi_norm_*'); File_BLA = File_BLA.name;
    File_VTA = dir('*_vta_epi_norm_*'); File_VTA = File_VTA.name;
    %File_whole = dir('*aal120_epi_*'); File_whole = File_whole.name;
    File_whole = dir('ts_aal2_bilateral_func_resting2_nosmoothed_denoised_mni_epinorm*'); File_whole = File_whole.name;
    
    % Load data
    formatSpec = '%f';
    fileID = fopen(File_mPFC,'r');
    Data_mPFC = fscanf(fileID, formatSpec);
    fileID = fopen(File_BLA,'r');
    Data_BLA = fscanf(fileID, formatSpec);
    fileID = fopen(File_VTA,'r');
    Data_VTA = fscanf(fileID, formatSpec);    
    fileID = fopen(File_whole,'r'); 
    fclose(fileID);
    Data_whole = load(File_whole);
    %Data_whole = reshape(Data_whole,120,[])';
    %Data_Amyg = mean(Data_whole(:,[45 46]),2); % 45 & 46 is Amyg (left & right)
    Data_Amyg = Data_whole(:,23); % 23 is Amyg
    %Data_Insular = mean(Data_whole(:,[33,34]),2); % 33 & 34 is Insular (left & right)
    Data_Insula = Data_whole(:,17); % 17 is Insula
    
    Data_areas{nSubj} = [Data_mPFC Data_BLA Data_Amyg Data_VTA Data_Insula];
    Idx_size(nSubj) = size(Data_areas{nSubj},1);
    clear File_BLA File_mPFC File_whole formatSpec fileID
    
    cd(Dirdata);
end

%% Trim signal length
min_size = min(Idx_size);
for nSubj = 1:size(BasicInfo,1)
   Data_areas{nSubj} = Data_areas{nSubj}(end-min_size+1:end,:); 
    
end


%% Calculate variables

for nSubj = 1:size(BasicInfo,1)
    nSubj
    %% Calculate Measures with raw signals
    [Raw_mPFC(nSubj,1:5)] = calculate_5variables(Data_areas{nSubj}(:,1),1);
    [Raw_BLA(nSubj,1:5)] = calculate_5variables(Data_areas{nSubj}(:,2),1);
    [Raw_Amyg(nSubj,1:5)] = calculate_5variables(Data_areas{nSubj}(:,3),1);
    [Raw_VTA(nSubj,1:5)] = calculate_5variables(Data_areas{nSubj}(:,4),1);
    [Raw_Insular(nSubj,1:5)] = calculate_5variables(Data_areas{nSubj}(:,5),1);
    
    %% Calculate Dynamic connectivty using sliding window (varying window size)
    Var_window = [50,40,30];
    Var_methods = {'sw','tsw','dcc'};
    
    for nWindow = 1:3
    windowsize = Var_window(nWindow);
    [TMP_Dfc_sw] = sliding_window(Data_areas{nSubj}, windowsize); % sliding window
    Idx_nan = isnan(TMP_Dfc_sw);
    TMP_Dfc_sw(Idx_nan) = [];
    TMP_Dfc_sw = reshape(TMP_Dfc_sw,5,5,[]);
    TMP_Dfc_sw = atanh(TMP_Dfc_sw); % Fisher normalization r-to-z
    clear Idx_nan
    
    [TMP_Dfc_tsw] = tapered_sliding_window(Data_areas{nSubj},windowsize,3); % tapered sliding window
    Idx_nan = isnan(TMP_Dfc_tsw);
    TMP_Dfc_tsw(Idx_nan) = [];
    TMP_Dfc_tsw = reshape(TMP_Dfc_tsw,5,5,[]);
    TMP_Dfc_tsw = atanh(TMP_Dfc_tsw); % Fisher normalization r-to-z

    [H,TMP_Dfc_dcc] = DCC_X_mex(Data_areas{nSubj},1,0); % DCC
    TMP_Dfc_dcc = atanh(TMP_Dfc_dcc); % Fisher normalization r-to-z
    
    for nMethod = 1:3
    % mPFC-BLA
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_mPFC_BLA(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(1,2,:)),2);']); 
    % mPFC-Amyg
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_mPFC_Amyg(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(1,3,:)),2);']); 
    % mPFC-VTA
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_mPFC_VTA(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(1,4,:)),2);']); 
    % mPFC-Insular
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_mPFC_Insular(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(1,5,:)),2);']); 
    % BLA-Amyg
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_BLA_Amyg(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(2,3,:)),2);']); 
    % BLA-VTA
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_BLA_VTA(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(2,4,:)),2);']); 
    % BLA-Insular
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_BLA_Insular(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(2,5,:)),2);']); 
    % Amyg-VTA
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_Amyg_VTA(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(3,4,:)),2);']); 
    % Amyg-Insular
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_Amyg_Insular(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(3,5,:)),2);']);     
    % VTA-Insular
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_VTA_Insular(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(4,5,:)),2);']); 
    
    end
    
    
    clear TMP_Dfc_sw TMP_Dfc_tsw TMP_Dfc_dcc
    
    end
    
    
    cd(Dirdata);
end

%% Save temporaly
save TMP_MSSD

%% Add results to the Table
% Divide the age group
Sub_agegroup1 = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 11.56;  % child
Sub_agegroup2 = BasicInfo.AgeAtScan > 11.56 & BasicInfo.AgeAtScan <= 14.44; % young adolescent
Sub_agegroup3 = BasicInfo.AgeAtScan > 14.44 & BasicInfo.AgeAtScan <= 19.71; % older adolescent
Sub_agegroup4 = BasicInfo.AgeAtScan > 19.71;                                % adult
Sub_agegroup1to2 = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 14.44; % child to young adolescent
Sub_agegroup1to3 = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 19.71; % child to older adolescent
Sub_agegroup = Sub_agegroup1 + Sub_agegroup2.*2 + Sub_agegroup3.*3 + Sub_agegroup4.*4;
fprintf('Child: %d / Younger Adolescent: %d / Older Adolescent: %d / Adult: %d \n',sum(Sub_agegroup1),sum(Sub_agegroup2),sum(Sub_agegroup3),sum(Sub_agegroup4));

Sub_Con = strcmp(BasicInfo.SubjectType,'CONTROL');
Sub_Pat = strcmp(BasicInfo.SubjectType,'PATIENT');
Sub_group = Sub_Con + Sub_Pat.*2;

BasicInfo = addvars(BasicInfo, Sub_agegroup, 'After','AgeAtScan');
BasicInfo = addvars(BasicInfo, Sub_group, 'After', 'SubjectType');
BasicInfo = addvars(BasicInfo, Sub_agegroup1to2, 'After', 'Sub_agegroup');
BasicInfo = addvars(BasicInfo, Sub_agegroup1to3, 'After', 'Sub_agegroup1to2');

%% Make Output Table
ResultTable_Raw_mPFC    = addvars(BasicInfo, Raw_mPFC);
ResultTable_Raw_Amyg    = addvars(BasicInfo, Raw_Amyg);
ResultTable_Raw_VTA     = addvars(BasicInfo, Raw_VTA);
ResultTable_Raw_Insular = addvars(BasicInfo, Raw_Insular);
ResultTable_TR50sw      = addvars(BasicInfo, TR50sw_mPFC_BLA, TR50sw_mPFC_Amyg, TR50sw_mPFC_VTA, TR50sw_mPFC_Insular, TR50sw_BLA_Amyg, TR50sw_BLA_VTA, TR50sw_BLA_Insular, TR50sw_Amyg_VTA, TR50sw_Amyg_Insular, TR50sw_VTA_Insular);
ResultTable_TR40sw      = addvars(BasicInfo, TR40sw_mPFC_BLA, TR40sw_mPFC_Amyg, TR40sw_mPFC_VTA, TR40sw_mPFC_Insular, TR40sw_BLA_Amyg, TR40sw_BLA_VTA, TR40sw_BLA_Insular, TR40sw_Amyg_VTA, TR40sw_Amyg_Insular, TR40sw_VTA_Insular);
ResultTable_TR30sw      = addvars(BasicInfo, TR30sw_mPFC_BLA, TR30sw_mPFC_Amyg, TR30sw_mPFC_VTA, TR30sw_mPFC_Insular, TR30sw_BLA_Amyg, TR30sw_BLA_VTA, TR30sw_BLA_Insular, TR30sw_Amyg_VTA, TR30sw_Amyg_Insular, TR30sw_VTA_Insular);
ResultTable_TR50tsw     = addvars(BasicInfo, TR50tsw_mPFC_BLA, TR50tsw_mPFC_Amyg, TR50tsw_mPFC_VTA, TR50tsw_mPFC_Insular, TR50tsw_BLA_Amyg, TR50tsw_BLA_VTA, TR50tsw_BLA_Insular, TR50tsw_Amyg_VTA, TR50tsw_Amyg_Insular, TR50tsw_VTA_Insular);
ResultTable_TR40tsw     = addvars(BasicInfo, TR40tsw_mPFC_BLA, TR40tsw_mPFC_Amyg, TR40tsw_mPFC_VTA, TR40tsw_mPFC_Insular, TR40tsw_BLA_Amyg, TR40tsw_BLA_VTA, TR40tsw_BLA_Insular, TR40tsw_Amyg_VTA, TR40tsw_Amyg_Insular, TR40tsw_VTA_Insular);
ResultTable_TR30tsw     = addvars(BasicInfo, TR30tsw_mPFC_BLA, TR30tsw_mPFC_Amyg, TR30tsw_mPFC_VTA, TR30tsw_mPFC_Insular, TR30tsw_BLA_Amyg, TR30tsw_BLA_VTA, TR30tsw_BLA_Insular, TR30tsw_Amyg_VTA, TR30tsw_Amyg_Insular, TR30tsw_VTA_Insular);
ResultTable_TR50dcc     = addvars(BasicInfo, TR50dcc_mPFC_BLA, TR50dcc_mPFC_Amyg, TR50dcc_mPFC_VTA, TR50dcc_mPFC_Insular, TR50dcc_BLA_Amyg, TR50dcc_BLA_VTA, TR50dcc_BLA_Insular, TR50dcc_Amyg_VTA, TR50dcc_Amyg_Insular, TR50dcc_VTA_Insular);
ResultTable_TR40dcc     = addvars(BasicInfo, TR40dcc_mPFC_BLA, TR40dcc_mPFC_Amyg, TR40dcc_mPFC_VTA, TR40dcc_mPFC_Insular, TR40dcc_BLA_Amyg, TR40dcc_BLA_VTA, TR40dcc_BLA_Insular, TR40dcc_Amyg_VTA, TR40dcc_Amyg_Insular, TR40dcc_VTA_Insular);
ResultTable_TR30dcc     = addvars(BasicInfo, TR30dcc_mPFC_BLA, TR30dcc_mPFC_Amyg, TR30dcc_mPFC_VTA, TR30dcc_mPFC_Insular, TR30dcc_BLA_Amyg, TR30dcc_BLA_VTA, TR30dcc_BLA_Insular, TR30dcc_Amyg_VTA, TR30dcc_Amyg_Insular, TR30dcc_VTA_Insular);

% Split Variable
ResultTable_Raw_mPFC    = splitvars(ResultTable_Raw_mPFC,19);
ResultTable_Raw_Amyg    = splitvars(ResultTable_Raw_Amyg,19);
ResultTable_Raw_VTA     = splitvars(ResultTable_Raw_VTA,19);
ResultTable_Raw_Insular = splitvars(ResultTable_Raw_Insular,19);
ResultTable_TR50sw      = splitvars(ResultTable_TR50sw,[19:28]);
ResultTable_TR40sw      = splitvars(ResultTable_TR40sw,[19:28]);
ResultTable_TR30sw      = splitvars(ResultTable_TR30sw,[19:28]);
ResultTable_TR50tsw     = splitvars(ResultTable_TR50tsw,[19:28]);
ResultTable_TR40tsw     = splitvars(ResultTable_TR40tsw,[19:28]);
ResultTable_TR30tsw     = splitvars(ResultTable_TR30tsw,[19:28]);
ResultTable_TR50dcc     = splitvars(ResultTable_TR50dcc,[19:28]);
ResultTable_TR40dcc     = splitvars(ResultTable_TR40dcc,[19:28]);
ResultTable_TR30dcc     = splitvars(ResultTable_TR30dcc,[19:28]);


%% Exclude outliers (Criteria : mean Correlation, 2.5SD)
Conditions = {'_Raw_mPFC','_Raw_Amyg','_Raw_VTA','_Raw_Insular','_TR50sw','_TR50tsw','_TR50dcc','_TR40sw','_TR40tsw','_TR40dcc','_TR30sw','_TR30tsw','_TR30dcc'};
Var_mean = {'Raw_mPFC_1','Raw_Amyg_1','Raw_VTA_1','Raw_Insular_1','TR50sw_mPFC_BLA_1','TR50tsw_mPFC_BLA_1','TR50dcc_mPFC_BLA_1','TR40sw_mPFC_BLA_1','TR40tsw_mPFC_BLA_1','TR40dcc_mPFC_BLA_1','TR30sw_mPFC_BLA_1','TR30tsw_mPFC_BLA_1','TR30dcc_mPFC_BLA_1'};
Var_mssd = {'Raw_mPFC_2','Raw_Amyg_2','Raw_VTA_2','Raw_Insular_2','TR50sw_mPFC_BLA_2','TR50tsw_mPFC_BLA_2','TR50dcc_mPFC_BLA_2','TR40sw_mPFC_BLA_2','TR40tsw_mPFC_BLA_2','TR40dcc_mPFC_BLA_2','TR30sw_mPFC_BLA_2','TR30tsw_mPFC_BLA_2','TR30dcc_mPFC_BLA_2'};

% for nCondition = 1:length(Conditions)
%     filename = ['ResultTable' Conditions{nCondition}];
%     eval(['TMP_std1 = std(' filename '.' Var_mean{nCondition} ');']);
%     eval(['TMP_mean1 = mean(' filename '.' Var_mean{nCondition} ');']);
%     eval(['Outlier1 = ' filename '.' Var_mean{nCondition} '> TMP_mean1+2.5*TMP_std1 | ' filename '.' Var_mean{nCondition} '< TMP_mean1-2.5*TMP_std1;']);
%     
% %     eval(['TMP_std2 = std(' filename '.' Var_mssd{nCondition} ');']);
% %     eval(['TMP_mean2 = mean(' filename '.' Var_mssd{nCondition} ');']);
% %     eval(['Outlier2 = ' filename '.' Var_mssd{nCondition} '> TMP_mean2+3*TMP_std2 | ' filename '.' Var_mssd{nCondition} '< TMP_mean2-3*TMP_std2;']);
%     
%       Outlier = Outlier1;
% %     Outlier = Outlier1 + Outlier2;
% %     Idx_over = find(Outlier>=2);
% %     if ~isempty(Idx_over), Outlier(Idx_over) = 1; end
% %     Outlier = logical(Outlier);
%     
%     eval([filename '(Outlier,:) = [];']);
%     
%     % save the information about excluded subjects
%     List_exclude{nCondition,1} = filename;
%     List_exclude{nCondition,2} = Outlier;
%     
%     clear Outlier1 Outlier2 TMP_mean1 TMP_std1 TMP_mean2 TMP_std2
% end
% 
% save List_exclude List_exclude

%% Save Results
cd(Dirdata); cd('..');
mkdir('Result_MSSD_fd025_trim');
cd('Result_MSSD_fd025_trim');

for nCondition = 1:length(Conditions)
    filename = ['ResultTable' Conditions{nCondition}];
    % Find index of age group again because we exclude outliers
    eval(['Sub_agegroup1 = ' filename '.AgeAtScan > 6 & ' filename '.AgeAtScan <= 11.56;']);
    eval(['Sub_agegroup2 = ' filename '.AgeAtScan > 11.56 & ' filename '.AgeAtScan <= 14.44;']);
    eval(['Sub_agegroup3 = ' filename '.AgeAtScan > 14.44 & ' filename '.AgeAtScan <= 19.71;']);
    eval(['Sub_agegroup4 = ' filename '.AgeAtScan > 19.71;']);
    eval(['Sub_agegroup1to2 = ' filename '.AgeAtScan > 6 & ' filename '.AgeAtScan <= 14.44;']);
    eval(['Sub_agegroup1to3 = ' filename '.AgeAtScan > 6 & ' filename '.AgeAtScan <= 19.71;']);
    
    eval(['mkdir(''' filename ''');']);
    eval(['cd(''' filename ''');']);
    
    % Divide age groups
    eval([filename '_agegroup1 = ' filename '(Sub_agegroup1,:);']);
    eval([filename '_agegroup2 = ' filename '(Sub_agegroup2,:);']);
    eval([filename '_agegroup3 = ' filename '(Sub_agegroup3,:);']);
    eval([filename '_agegroup4 = ' filename '(Sub_agegroup4,:);']);
    eval([filename '_agegroup1to2 = ' filename '(Sub_agegroup1to2,:);']);
    eval([filename '_agegroup1to3 = ' filename '(Sub_agegroup1to3,:);']);
    
    eval(['writetable(' filename ');']); % save as text file
    eval(['writetable(' filename '_agegroup1);']);
    eval(['writetable(' filename '_agegroup2);']);
    eval(['writetable(' filename '_agegroup3);']);
    eval(['writetable(' filename '_agegroup4);']);
    eval(['writetable(' filename '_agegroup1to2);']);
    eval(['writetable(' filename '_agegroup1to3);']);
    
%     % Plotting and save
%     if nCondition < 4, option1 = 'Raw'; else option1 = 'Dfc'; end
%     
%     mkdir('whole'); cd('whole');
%     eval(['plot_5variables(' filename ', ''' option1 ''', ''' Idx_variable{nCondition} ''');']);
%     cd('..');
%     mkdir('child'); cd('child');
%     eval(['plot_5variables(' filename '_child' ', ''' option1 ''', ''' Idx_variable{nCondition} ''');']);
%     cd('..');
%     mkdir('adolescent'); cd('adolescent');
%     eval(['plot_5variables(' filename '_adolescent' ', ''' option1 ''', ''' Idx_variable{nCondition} ''');']);
%     cd('..');        
%     mkdir('younger'); cd('younger');
%     eval(['plot_5variables(' filename '_younger' ', ''' option1 ''', ''' Idx_variable{nCondition} ''');']);
%     cd('..');        
%     mkdir('middleadult'); cd('middleadult');
%     eval(['plot_5variables(' filename '_middleadult' ', ''' option1 ''', ''' Idx_variable{nCondition} ''');']);
%     cd('..');  
%     mkdir('olderadult'); cd('olderadult');
%     eval(['plot_5variables(' filename '_olderadult' ', ''' option1 ''', ''' Idx_variable{nCondition} ''');']);
%     cd('..');  
%     mkdir('childtoadolescent'); cd('childtoadolescent');
%     eval(['plot_5variables(' filename '_childtoadolescent' ', ''' option1 ''', ''' Idx_variable{nCondition} ''');']);
%     cd('..');       
%     mkdir('childtoyounger'); cd('childtoyounger');
%     eval(['plot_5variables(' filename '_childtoyounger' ', ''' option1 ''', ''' Idx_variable{nCondition} ''');']);
%     cd('..');      
    
    
    
    cd('..');
end

disp('==========Finish=========');


