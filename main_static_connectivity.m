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

%% Parameters
p.allROI = [1 8 12 13 16]; % LAIRD ICs
p.outROI = [19,20]; % 19,20 are noise ICs
p.incROI = setdiff(p.allROI, p.outROI);

min_size = 110;

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
    Data_postcentral = Data_whole(:,31);
    
    % load LAIRD time series
    cd('LAIRD2011');
    List_txt = dir('ts_IC*');
    for nIC = 1:length(p.incROI)
        TMP_signal = load(List_txt(p.incROI(nIC)).name);
        LAIRD_whole(:,nIC) = TMP_signal';
        TMP_signal = [];
    end
    cd('..');
    
    Data_limbic = LAIRD_whole(:,1); 
    Data_motor = LAIRD_whole(:,2);
    Data_visual = LAIRD_whole(:,3);
    Data_dmn = LAIRD_whole(:,4);
    Data_auditory = LAIRD_whole(:,5);
    
    Data_areas = [Data_mPFC Data_BLA Data_Amyg Data_calcarine Data_postcentral Data_limbic Data_motor Data_visual Data_dmn Data_auditory];
%     Data_areas = Data_areas(end-min_size+1:end,:); 
    clear File_BLA File_mPFC File_whole formatSpec fileID LAIRD_whole

    
    %% Calculate Dynamic connectivty using sliding window (varying window size)

    TMP_sFC = corr(Data_areas);
    TMP_sFC = atanh(TMP_sFC); % Fisher normalization r-to-z

    % mPFC-BLA
    eval(['mPFC_BLA(nSubj,1) = TMP_sFC(1,2);']); 
    % mPFC-Amyg
    eval(['mPFC_Amyg(nSubj,1) = TMP_sFC(1,3);']); 
    % mPFC-calcarine
    eval(['mPFC_calcarine(nSubj,1) = TMP_sFC(1,4);']); 
    % mPFC-postcentral
    eval(['mPFC_postcentral(nSubj,1) = TMP_sFC(1,5);']); 
    
    % Amyg-BLA
    eval(['Amyg_BLA(nSubj,1) = TMP_sFC(3,2);']); 
    % Amyg-calcarine
    eval(['Amyg_calcarine(nSubj,1) = TMP_sFC(3,4);']); 
    % Amyg-postcentral
    eval(['Amyg_postcentral(nSubj,1) = TMP_sFC(3,5);']); 

    % limbic-DMN
    eval(['limbic_DMN(nSubj,1) = TMP_sFC(6,9);']); 
    % limbic-visual
    eval(['limbic_visual(nSubj,1) = TMP_sFC(6,8);']); 
    % limbic-motor
    eval(['limbic_motor(nSubj,1) = TMP_sFC(6,7);']); 

    % DMN-visual
    eval(['DMN_visual(nSubj,1) = TMP_sFC(9,8);']); 
    % DMN-motor
    eval(['DMN_motor(nSubj,1) = TMP_sFC(9,7);']); 
    
    clear TMP_sFC 
    
    cd(Dirdata);
end


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
ResultTable_sFC      = addvars(BasicInfo, mPFC_BLA, mPFC_Amyg, mPFC_calcarine, mPFC_postcentral, Amyg_BLA, Amyg_calcarine, Amyg_postcentral, limbic_DMN, limbic_visual, limbic_motor, DMN_visual, DMN_motor);

%% Exclude outliers (Criteria : mean Correlation, 2.5SD)
% Conditions = {'_TR50sw','_TR50tsw','_TR40sw','_TR40tsw','_TR30sw','_TR30tsw'};
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
mkdir('Result_sFC_20210329');
cd('Result_sFC_20210329');

    filename = ['ResultTable_sFC'];
    % Find index of age group again because we exclude outliers
    eval(['Sub_agegroup1 = ' filename '.AgeAtScan > 6 & ' filename '.AgeAtScan <= 11.56;']);
    eval(['Sub_agegroup2 = ' filename '.AgeAtScan > 11.56 & ' filename '.AgeAtScan <= 14.44;']);
    eval(['Sub_agegroup3 = ' filename '.AgeAtScan > 14.44 & ' filename '.AgeAtScan <= 19.71;']);
    eval(['Sub_agegroup4 = ' filename '.AgeAtScan > 19.71;']);
    eval(['Sub_agegroup1to2 = ' filename '.AgeAtScan > 6 & ' filename '.AgeAtScan <= 14.44;']);
    eval(['Sub_agegroup1to3 = ' filename '.AgeAtScan > 6 & ' filename '.AgeAtScan <= 19.71;']);
        
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
    
    

disp('==========Finish=========');


