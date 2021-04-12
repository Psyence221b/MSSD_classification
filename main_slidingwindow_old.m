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
    Fol_subj = BasicInfo.subID(nSubj);
    cd(Fol_subj{1});
    File_BLA = dir('*_BLA_*'); File_BLA = File_BLA.name;
    File_mPFC = dir('*_mpfc_*'); File_mPFC = File_mPFC.name;
    File_VTA = dir('*_vta_*'); File_VTA = File_VTA.name;
    File_whole = dir('*aal120_epi_*'); File_whole = File_whole.name;
    
    % Load data
    formatSpec = '%f';
    fileID = fopen(File_mPFC,'r');
    Data_mPFC = fscanf(fileID, formatSpec);
    fileID = fopen(File_BLA,'r');
    Data_BLA = fscanf(fileID, formatSpec);
    fileID = fopen(File_VTA,'r');
    Data_VTA = fscanf(fileID, formatSpec);    
    fileID = fopen(File_whole,'r'); 
    Data_whole = fscanf(fileID, formatSpec);
    fclose(fileID);
    Data_whole = reshape(Data_whole,120,[])';
    Data_Amyg = mean(Data_whole(:,[45 46]),2); % 45&46 is Amyg (left & right)
    Data_Insular = mean(Data_whole(:,[33,34]),2); % 33&34 is Insular (left & right)
    
    Data_withBLA = [Data_mPFC Data_BLA];
    Data_withAmyg = [Data_mPFC Data_Amyg];
    clear File_BLA File_mPFC File_whole formatSpec fileID
    
    %% Calculate Measures with raw signals
    % mPFC
    [mPFC_mean(nSubj,1), mPFC_MSSD1(nSubj,1), mPFC_MSSD2(nSubj,1), mPFC_SD(nSubj,1), mPFC_VSD(nSubj,1)] = calculate_5variables(Data_mPFC,1);
    % BLA
    [BLA_mean(nSubj,1), BLA_MSSD1(nSubj,1), BLA_MSSD2(nSubj,1), BLA_SD(nSubj,1), BLA_VSD(nSubj,1)] = calculate_5variables(Data_BLA,1);
    % Amyg
    [Amyg_mean(nSubj,1), Amyg_MSSD1(nSubj,1), Amyg_MSSD2(nSubj,1), Amyg_SD(nSubj,1), Amyg_VSD(nSubj,1)] = calculate_5variables(Data_Amyg,1);  
    
    %% Calculate Dynamic connectivty using sliding window (with BLA)
    windowsize = 50; % follow the manuscript of ASD
    [TMP_Dfc] = sliding_window(Data_withBLA, windowsize);
    Dfc_BLA_50TR = squeeze(TMP_Dfc(1,2,:));
    Dfc_BLA_50TR = rmmissing(Dfc_BLA_50TR);
    Dfc_BLA_50TR = atanh(Dfc_BLA_50TR); % Fisher normalization r-to-z
    [Dfc_BLA_50TR_mean(nSubj,1), Dfc_BLA_50TR_MSSD1(nSubj,1), Dfc_BLA_50TR_MSSD2(nSubj,1), Dfc_BLA_50TR_SD(nSubj,1), Dfc_BLA_50TR_VSD(nSubj,1)] = calculate_5variables(Dfc_BLA_50TR,2);  
    clear TMP_Dfc
    
    windowsize = 40; % follow the manuscript of ASD
    [TMP_Dfc] = sliding_window(Data_withBLA, windowsize);
    Dfc_BLA_40TR = squeeze(TMP_Dfc(1,2,:));
    Dfc_BLA_40TR = rmmissing(Dfc_BLA_40TR);
    Dfc_BLA_40TR = atanh(Dfc_BLA_40TR); % Fisher normalization r-to-z
    [Dfc_BLA_40TR_mean(nSubj,1), Dfc_BLA_40TR_MSSD1(nSubj,1), Dfc_BLA_40TR_MSSD2(nSubj,1), Dfc_BLA_40TR_SD(nSubj,1), Dfc_BLA_40TR_VSD(nSubj,1)] = calculate_5variables(Dfc_BLA_40TR,2);  
    clear TMP_Dfc

    windowsize = 30; % follow the manuscript of ASD
    [TMP_Dfc] = sliding_window(Data_withBLA, windowsize);
    Dfc_BLA_30TR = squeeze(TMP_Dfc(1,2,:));
    Dfc_BLA_30TR = rmmissing(Dfc_BLA_30TR);
    Dfc_BLA_30TR = atanh(Dfc_BLA_30TR); % Fisher normalization r-to-z
    [Dfc_BLA_30TR_mean(nSubj,1), Dfc_BLA_30TR_MSSD1(nSubj,1), Dfc_BLA_30TR_MSSD2(nSubj,1), Dfc_BLA_30TR_SD(nSubj,1), Dfc_BLA_30TR_VSD(nSubj,1)] = calculate_5variables(Dfc_BLA_30TR,2);  
    clear TMP_Dfc
    
    %% Calculate Dynamic connectivty using sliding window (with Amyg)
    windowsize = 50; % follow the manuscript of ASD
    [TMP_Dfc] = sliding_window(Data_withAmyg, windowsize);
    Dfc_Amyg_50TR = squeeze(TMP_Dfc(1,2,:));
    Dfc_Amyg_50TR = rmmissing(Dfc_Amyg_50TR);
    Dfc_Amyg_50TR = atanh(Dfc_Amyg_50TR); % Fisher normalization r-to-z
    [Dfc_Amyg_50TR_mean(nSubj,1), Dfc_Amyg_50TR_MSSD1(nSubj,1), Dfc_Amyg_50TR_MSSD2(nSubj,1), Dfc_Amyg_50TR_SD(nSubj,1), Dfc_Amyg_50TR_VSD(nSubj,1)] = calculate_5variables(Dfc_Amyg_50TR,2);  
    clear TMP_Dfc
    
    windowsize = 40; % follow the manuscript of ASD
    [TMP_Dfc] = sliding_window(Data_withAmyg, windowsize);
    Dfc_Amyg_40TR = squeeze(TMP_Dfc(1,2,:));
    Dfc_Amyg_40TR = rmmissing(Dfc_Amyg_40TR);
    Dfc_Amyg_40TR = atanh(Dfc_Amyg_40TR); % Fisher normalization r-to-z
    [Dfc_Amyg_40TR_mean(nSubj,1), Dfc_Amyg_40TR_MSSD1(nSubj,1), Dfc_Amyg_40TR_MSSD2(nSubj,1), Dfc_Amyg_40TR_SD(nSubj,1), Dfc_Amyg_40TR_VSD(nSubj,1)] = calculate_5variables(Dfc_Amyg_40TR,2);  
    clear TMP_Dfc

    windowsize = 30; % follow the manuscript of ASD
    [TMP_Dfc] = sliding_window(Data_withAmyg, windowsize);
    Dfc_Amyg_30TR = squeeze(TMP_Dfc(1,2,:));
    Dfc_Amyg_30TR = rmmissing(Dfc_Amyg_30TR);
    Dfc_Amyg_30TR = atanh(Dfc_Amyg_30TR); % Fisher normalization r-to-z
    [Dfc_Amyg_30TR_mean(nSubj,1), Dfc_Amyg_30TR_MSSD1(nSubj,1), Dfc_Amyg_30TR_MSSD2(nSubj,1), Dfc_Amyg_30TR_SD(nSubj,1), Dfc_Amyg_30TR_VSD(nSubj,1)] = calculate_5variables(Dfc_Amyg_30TR,2);  
    clear TMP_Dfc
    
    
    cd(Dirdata);
end

%% Add results to the Table
% add variable of age group
Sub_child = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 12;
Sub_adolescent = BasicInfo.AgeAtScan > 12 & BasicInfo.AgeAtScan <= 19;
Sub_younger = BasicInfo.AgeAtScan > 19 & BasicInfo.AgeAtScan <= 25;
Sub_middleadult = BasicInfo.AgeAtScan > 25 & BasicInfo.AgeAtScan <= 55;
Sub_olderadult = BasicInfo.AgeAtScan > 55;
Sub_childtoadolescent = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 19;
Sub_childtoyounger = BasicInfo.AgeAtScan <= 25;
Sub_agegroup = Sub_child + Sub_adolescent.*2 + Sub_younger.*3 + Sub_middleadult.*4 + Sub_olderadult.*5;
fprintf('Child: %d / Adolescent: %d / Younger: %d / MiddleAdult: %d / OlderAdult: %d',sum(Sub_child),sum(Sub_adolescent),sum(Sub_younger),sum(Sub_middleadult),sum(Sub_olderadult));

Sub_Con = strcmp(BasicInfo.SubjectType,'CONTROL');
Sub_Pat = strcmp(BasicInfo.SubjectType,'PATIENT');
Sub_group = Sub_Con + Sub_Pat.*2;

BasicInfo = addvars(BasicInfo, Sub_agegroup, 'After','AgeAtScan');
BasicInfo = addvars(BasicInfo, Sub_group, 'After', 'SubjectType');

% Make Output Table
ResultTable_Raw_mPFC  = addvars(BasicInfo, mPFC_mean, mPFC_MSSD1, mPFC_MSSD2, mPFC_SD, mPFC_VSD);
ResultTable_Raw_BLA   = addvars(BasicInfo, BLA_mean, BLA_MSSD1, BLA_MSSD2, BLA_SD, BLA_VSD);
ResultTable_Raw_Amyg  = addvars(BasicInfo, Amyg_mean, Amyg_MSSD1, Amyg_MSSD2, Amyg_SD, Amyg_VSD);
ResultTable_BLA_50TR  = addvars(BasicInfo, Dfc_BLA_50TR_mean, Dfc_BLA_50TR_MSSD1, Dfc_BLA_50TR_MSSD2, Dfc_BLA_50TR_SD, Dfc_BLA_50TR_VSD);
ResultTable_BLA_40TR  = addvars(BasicInfo, Dfc_BLA_40TR_mean, Dfc_BLA_40TR_MSSD1, Dfc_BLA_40TR_MSSD2, Dfc_BLA_40TR_SD, Dfc_BLA_40TR_VSD);
ResultTable_BLA_30TR  = addvars(BasicInfo, Dfc_BLA_30TR_mean, Dfc_BLA_30TR_MSSD1, Dfc_BLA_30TR_MSSD2, Dfc_BLA_30TR_SD, Dfc_BLA_30TR_VSD);
ResultTable_Amyg_50TR = addvars(BasicInfo, Dfc_Amyg_50TR_mean, Dfc_Amyg_50TR_MSSD1, Dfc_Amyg_50TR_MSSD2, Dfc_Amyg_50TR_SD, Dfc_Amyg_50TR_VSD);
ResultTable_Amyg_40TR = addvars(BasicInfo, Dfc_Amyg_40TR_mean, Dfc_Amyg_40TR_MSSD1, Dfc_Amyg_40TR_MSSD2, Dfc_Amyg_40TR_SD, Dfc_Amyg_40TR_VSD);
ResultTable_Amyg_30TR = addvars(BasicInfo, Dfc_Amyg_30TR_mean, Dfc_Amyg_30TR_MSSD1, Dfc_Amyg_30TR_MSSD2, Dfc_Amyg_30TR_SD, Dfc_Amyg_30TR_VSD);

%% Exclude outliers (Criteria : MSSD1, 2.5 SD)
Conditions = {'_Raw_mPFC','_Raw_BLA','_Raw_Amyg','_BLA_50TR','_BLA_40TR','_BLA_30TR','_Amyg_50TR','_Amyg_40TR','_Amyg_30TR'};
Var_MSSD = {'mPFC_MSSD1','BLA_MSSD1','Amyg_MSSD1','Dfc_BLA_50TR_MSSD1','Dfc_BLA_40TR_MSSD1','Dfc_BLA_30TR_MSSD1','Dfc_Amyg_50TR_MSSD1','Dfc_Amyg_40TR_MSSD1','Dfc_Amyg_30TR_MSSD1'};
for ncondition = 1:length(Conditions)
    filename = ['ResultTable' Conditions{ncondition}];
    
    eval(['TMP_std = std(' filename '.' Var_MSSD{ncondition} ');']);
    eval(['TMP_mean = mean(' filename '.' Var_MSSD{ncondition} ');']);
    eval(['Outlier = ' filename '.' Var_MSSD{ncondition} '> TMP_mean+2.5*TMP_std | ' filename '.' Var_MSSD{ncondition} '< TMP_mean-2.5*TMP_std;']);
    eval([filename '(Outlier,:) = [];']);
    clear Outlier TMP_mean TMP_std
end

%% Save Results
cd(Dirdata); cd('..');
mkdir('MSSD_mPFC');
cd('MSSD_mPFC');

Conditions = {'_Raw_mPFC','_Raw_BLA','_Raw_Amyg','_BLA_50TR','_BLA_40TR','_BLA_30TR','_Amyg_50TR','_Amyg_40TR','_Amyg_30TR'};
Idx_variable = {'mPFC','BLA','Amyg','Dfc_BLA_50TR','Dfc_BLA_40TR','Dfc_BLA_30TR','Dfc_Amyg_50TR','Dfc_Amyg_40TR','Dfc_Amyg_30TR'};
for ncondition = 1:length(Conditions)
    filename = ['ResultTable' Conditions{ncondition}];
    % Find index of age group again because we exclude outliers
    eval(['Sub_child = ' filename '.AgeAtScan > 6 & ' filename '.AgeAtScan <= 12;']);
    eval(['Sub_adolescent = ' filename '.AgeAtScan > 12 & ' filename '.AgeAtScan <= 19;']);
    eval(['Sub_younger = ' filename '.AgeAtScan > 19 & ' filename '.AgeAtScan <= 25;']);
    eval(['Sub_middleadult = ' filename '.AgeAtScan > 25 & ' filename '.AgeAtScan <= 55;']);
    eval(['Sub_olderadult = ' filename '.AgeAtScan > 55;']);
    eval(['Sub_childtoadolescent = ' filename '.AgeAtScan > 6 & ' filename '.AgeAtScan <= 19;']);
    eval(['Sub_childtoyounger = ' filename '.AgeAtScan <= 25;']);
    
    eval(['mkdir(''' filename ''');']);
    eval(['cd(''' filename ''');']);
    
    % Divide age groups
    eval([filename '_child = ' filename '(Sub_child,:);']);
    eval([filename '_adolescent = ' filename '(Sub_adolescent,:);']);
    eval([filename '_younger = ' filename '(Sub_younger,:);']);
    eval([filename '_middleadult = ' filename '(Sub_middleadult,:);']);
    eval([filename '_olderadult = ' filename '(Sub_olderadult,:);']);
    eval([filename '_childtoadolescent = ' filename '(Sub_childtoadolescent,:);']);
    eval([filename '_childtoyounger = ' filename '(Sub_childtoyounger,:);']);
    
    eval(['writetable(' filename ');']); % save as text file
    eval(['writetable(' filename '_child);']);
    eval(['writetable(' filename '_adolescent);']);
    eval(['writetable(' filename '_younger);']);
    eval(['writetable(' filename '_middleadult);']);
    eval(['writetable(' filename '_olderadult);']);
    eval(['writetable(' filename '_childtoadolescent);']);
    eval(['writetable(' filename '_childtoyounger);']);
    
    % Plotting and save
    if ncondition < 4, option1 = 'Raw'; else option1 = 'Dfc'; end
    
    mkdir('whole'); cd('whole');
    eval(['plot_5variables(' filename ', ''' option1 ''', ''' Idx_variable{ncondition} ''');']);
    cd('..');
    mkdir('child'); cd('child');
    eval(['plot_5variables(' filename '_child' ', ''' option1 ''', ''' Idx_variable{ncondition} ''');']);
    cd('..');
    mkdir('adolescent'); cd('adolescent');
    eval(['plot_5variables(' filename '_adolescent' ', ''' option1 ''', ''' Idx_variable{ncondition} ''');']);
    cd('..');        
    mkdir('younger'); cd('younger');
    eval(['plot_5variables(' filename '_younger' ', ''' option1 ''', ''' Idx_variable{ncondition} ''');']);
    cd('..');        
    mkdir('middleadult'); cd('middleadult');
    eval(['plot_5variables(' filename '_middleadult' ', ''' option1 ''', ''' Idx_variable{ncondition} ''');']);
    cd('..');  
    mkdir('olderadult'); cd('olderadult');
    eval(['plot_5variables(' filename '_olderadult' ', ''' option1 ''', ''' Idx_variable{ncondition} ''');']);
    cd('..');  
    mkdir('childtoadolescent'); cd('childtoadolescent');
    eval(['plot_5variables(' filename '_childtoadolescent' ', ''' option1 ''', ''' Idx_variable{ncondition} ''');']);
    cd('..');       
    mkdir('childtoyounger'); cd('childtoyounger');
    eval(['plot_5variables(' filename '_childtoyounger' ', ''' option1 ''', ''' Idx_variable{ncondition} ''');']);
    cd('..');      
    
    
    
    cd('..');
end



