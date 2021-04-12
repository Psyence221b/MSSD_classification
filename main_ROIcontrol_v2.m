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
    File_mPFC = dir('*_mpfc_epi_norm_*'); File_mPFC = File_mPFC.name;
    File_BLA = dir('*_BLA_epi_norm_*'); File_BLA = File_BLA.name;
    %File_whole = dir('*aal120_epi_*'); File_whole = File_whole.name;
    File_whole = dir('ts_aal2_bilateral_func_resting2_nosmoothed_denoised_mni_epinorm*'); File_whole = File_whole.name;
    
    % Load data
    formatSpec = '%f';
    fileID = fopen(File_mPFC,'r');
    Data_mPFC = fscanf(fileID, formatSpec);
    fileID = fopen(File_BLA,'r');
    Data_BLA = fscanf(fileID, formatSpec);    
    fclose(fileID);
    
    Data_whole = load(File_whole);
    Data_Amyg = Data_whole(:,23); % 23 is Amyg
    Data_calcarine = Data_whole(:,24); % calcarine is V1
    Data_heschl = Data_whole(:,42); % auditory 
    Data_postcentral = Data_whole(:,31);
    Data_olfactory = Data_whole(:,9);
    
    Data_areas = [Data_mPFC Data_BLA Data_Amyg Data_calcarine Data_heschl Data_postcentral Data_olfactory];
    clear File_BLA File_mPFC File_whole formatSpec fileID
    
    %% Calculate Measures with raw signals
    [Raw_mPFC(nSubj,1:5)] = calculate_5variables(Data_mPFC,1);
    [Raw_BLA(nSubj,1:5)] = calculate_5variables(Data_BLA,1);
    [Raw_Amyg(nSubj,1:5)] = calculate_5variables(Data_Amyg,1);
    [Raw_calcarine(nSubj,1:5)] = calculate_5variables(Data_calcarine,1);
    [Raw_heschl(nSubj,1:5)] = calculate_5variables(Data_heschl,1);
    [Raw_postcentral(nSubj,1:5)] = calculate_5variables(Data_postcentral,1);
    [Raw_olfactory(nSubj,1:5)] = calculate_5variables(Data_olfactory,1);
    
    %% Calculate Dynamic connectivty using sliding window (varying window size)
    Var_window = [50,40,30];
    Var_methods = {'sw','tsw','dcc'};
    
    for nWindow = 1:3
    windowsize = Var_window(nWindow);
    [TMP_Dfc_sw] = sliding_window(Data_areas, windowsize); % sliding window
    Idx_nan = isnan(TMP_Dfc_sw);
    TMP_Dfc_sw(Idx_nan) = [];
    TMP_Dfc_sw = reshape(TMP_Dfc_sw,size(Data_areas,2),size(Data_areas,2),[]);
    TMP_Dfc_sw = atanh(TMP_Dfc_sw); % Fisher normalization r-to-z
    clear Idx_nan
    
    [TMP_Dfc_tsw] = tapered_sliding_window(Data_areas,windowsize,3); % tapered sliding window
    Idx_nan = isnan(TMP_Dfc_tsw);
    TMP_Dfc_tsw(Idx_nan) = [];
    TMP_Dfc_tsw = reshape(TMP_Dfc_tsw,size(Data_areas,2),size(Data_areas,2),[]);
    TMP_Dfc_tsw = atanh(TMP_Dfc_tsw); % Fisher normalization r-to-z

    [H,TMP_Dfc_dcc] = DCC_X_mex(Data_areas,1,0); % DCC
    TMP_Dfc_dcc = atanh(TMP_Dfc_dcc); % Fisher normalization r-to-z
    
    for nMethod = 1:3
    % mPFC-BLA
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_mPFC_BLA(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(1,2,:)),2);']); 
    % mPFC-Amyg
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_mPFC_Amyg(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(1,3,:)),2);']); 
    % mPFC-calcarine
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_mPFC_calcarine(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(1,4,:)),2);']); 
    % mPFC-heschl
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_mPFC_heschl(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(1,5,:)),2);']);     
    % mPFC-postcentral
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_mPFC_postcentral(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(1,6,:)),2);']); 
    % mPFC-olfactory
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_mPFC_olfactory(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(1,7,:)),2);']); 
    
    % BLA-Amyg
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_BLA_Amyg(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(2,3,:)),2);']); 
    % BLA-calcarine
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_BLA_calcarine(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(2,4,:)),2);']); 
    % BLA-heschl
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_BLA_heschl(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(2,5,:)),2);']);     
    % BLA-postcentral
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_BLA_postcentral(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(2,6,:)),2);']); 
    % BLA-olfactory
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_BLA_olfactory(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(2,7,:)),2);']); 

    % Amyg-calcarine
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_Amyg_calcarine(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(3,4,:)),2);']); 
    % Amyg-heschl
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_Amyg_heschl(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(3,5,:)),2);']);     
    % Amyg-postcentral
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_Amyg_postcentral(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(3,6,:)),2);']); 
    % Amyg-olfactory
    eval(['[TR' num2str(windowsize) Var_methods{nMethod} '_Amyg_olfactory(nSubj,1:5)] = calculate_5variables(squeeze(TMP_Dfc_' Var_methods{nMethod} '(3,7,:)),2);']); 

    end
    
    
    clear TMP_Dfc_sw TMP_Dfc_tsw TMP_Dfc_dcc
    
    end
    
    
    cd(Dirdata);
end

%% Save temporaly
save TMP_MSSD_ROIcontrol2

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
ResultTable_TR50sw      = addvars(BasicInfo, TR50sw_mPFC_BLA, TR50sw_mPFC_Amyg, TR50sw_mPFC_calcarine, TR50sw_mPFC_heschl, TR50sw_mPFC_postcentral, TR50sw_mPFC_olfactory, TR50sw_BLA_Amyg, TR50sw_BLA_calcarine, TR50sw_BLA_heschl, TR50sw_BLA_postcentral, TR50sw_BLA_olfactory, TR50sw_Amyg_calcarine, TR50sw_Amyg_heschl, TR50sw_Amyg_postcentral, TR50sw_Amyg_olfactory);
ResultTable_TR40sw      = addvars(BasicInfo, TR40sw_mPFC_BLA, TR40sw_mPFC_Amyg, TR40sw_mPFC_calcarine, TR40sw_mPFC_heschl, TR40sw_mPFC_postcentral, TR40sw_mPFC_olfactory, TR40sw_BLA_Amyg, TR40sw_BLA_calcarine, TR40sw_BLA_heschl, TR40sw_BLA_postcentral, TR40sw_BLA_olfactory, TR40sw_Amyg_calcarine, TR40sw_Amyg_heschl, TR40sw_Amyg_postcentral, TR40sw_Amyg_olfactory); 
ResultTable_TR30sw      = addvars(BasicInfo, TR30sw_mPFC_BLA, TR30sw_mPFC_Amyg, TR30sw_mPFC_calcarine, TR30sw_mPFC_heschl, TR30sw_mPFC_postcentral, TR30sw_mPFC_olfactory, TR30sw_BLA_Amyg, TR30sw_BLA_calcarine, TR30sw_BLA_heschl, TR30sw_BLA_postcentral, TR30sw_BLA_olfactory, TR30sw_Amyg_calcarine, TR30sw_Amyg_heschl, TR30sw_Amyg_postcentral, TR30sw_Amyg_olfactory);

ResultTable_TR50tsw      = addvars(BasicInfo, TR50tsw_mPFC_BLA, TR50tsw_mPFC_Amyg, TR50tsw_mPFC_calcarine, TR50tsw_mPFC_heschl, TR50tsw_mPFC_postcentral, TR50tsw_mPFC_olfactory, TR50tsw_BLA_Amyg, TR50tsw_BLA_calcarine, TR50tsw_BLA_heschl, TR50tsw_BLA_postcentral, TR50tsw_BLA_olfactory, TR50tsw_Amyg_calcarine, TR50tsw_Amyg_heschl, TR50tsw_Amyg_postcentral, TR50tsw_Amyg_olfactory);
ResultTable_TR40tsw      = addvars(BasicInfo, TR40tsw_mPFC_BLA, TR40tsw_mPFC_Amyg, TR40tsw_mPFC_calcarine, TR40tsw_mPFC_heschl, TR40tsw_mPFC_postcentral, TR40tsw_mPFC_olfactory, TR40tsw_BLA_Amyg, TR40tsw_BLA_calcarine, TR40tsw_BLA_heschl, TR40tsw_BLA_postcentral, TR40tsw_BLA_olfactory, TR40tsw_Amyg_calcarine, TR40tsw_Amyg_heschl, TR40tsw_Amyg_postcentral, TR40tsw_Amyg_olfactory); 
ResultTable_TR30tsw      = addvars(BasicInfo, TR30tsw_mPFC_BLA, TR30tsw_mPFC_Amyg, TR30tsw_mPFC_calcarine, TR30tsw_mPFC_heschl, TR30tsw_mPFC_postcentral, TR30tsw_mPFC_olfactory, TR30tsw_BLA_Amyg, TR30tsw_BLA_calcarine, TR30tsw_BLA_heschl, TR30tsw_BLA_postcentral, TR30tsw_BLA_olfactory, TR30tsw_Amyg_calcarine, TR30tsw_Amyg_heschl, TR30tsw_Amyg_postcentral, TR30tsw_Amyg_olfactory);


% Split Variable
ResultTable_TR50sw      = splitvars(ResultTable_TR50sw,[21:35]);
ResultTable_TR40sw      = splitvars(ResultTable_TR40sw,[21:35]);
ResultTable_TR30sw      = splitvars(ResultTable_TR30sw,[21:35]);
ResultTable_TR50tsw     = splitvars(ResultTable_TR50tsw,[21:35]);
ResultTable_TR40tsw     = splitvars(ResultTable_TR40tsw,[21:35]);
ResultTable_TR30tsw     = splitvars(ResultTable_TR30tsw,[21:35]);



%% Exclude outliers (Criteria : mean Correlation, 2.5SD)
Conditions = {'_TR50sw','_TR50tsw','_TR40sw','_TR40tsw','_TR30sw','_TR30tsw'};
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
mkdir('Result_MSSD_ROIcontrol_v2');
cd('Result_MSSD_ROIcontrol_v2');

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


