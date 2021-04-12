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

%% Make null variable
Whole_mPFC_Amyg = [];
Whole_mPFC_VTA = [];
Whole_mPFC_Insular = [];
Whole_VTA_Insular = [];

%% Get Data from 'analysis' folder
cd(Dirdata);

for nSubj = 1:size(BasicInfo,1)
    nSubj
    Fol_subj = BasicInfo.subID(nSubj);
    cd(Fol_subj{1});
%     File_BLA = dir('*_BLA_*'); File_BLA = File_BLA.name;
    File_mPFC = dir('*_mpfc_*'); File_mPFC = File_mPFC.name;
    File_BLA = dir('*_BLA_*'); File_BLA = File_BLA.name;
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
    
    Data_areas = [Data_mPFC Data_BLA Data_Amyg Data_VTA Data_Insular];
    clear File_BLA File_mPFC File_whole formatSpec fileID
    
    %% Calculate Dynamic connectivty using sliding window (varying window size)
    Var_window = [40];
    Var_methods = {'sw'};
    
    for nWindow = 1:1
    windowsize = Var_window(nWindow);
    [TMP_Dfc_sw] = sliding_window(Data_areas, windowsize); % sliding window
    Idx_nan = isnan(TMP_Dfc_sw);
    TMP_Dfc_sw(Idx_nan) = [];
    TMP_Dfc_sw = reshape(TMP_Dfc_sw,5,5,[]);
    TMP_Dfc_sw = atanh(TMP_Dfc_sw); % Fisher normalization r-to-z
    clear Idx_nan
    end
    Whole_mPFC_Amyg(nSubj,:) = resample(TMP_Dfc_sw(1,3,:),200,size(TMP_Dfc_sw,3));
    Whole_mPFC_VTA(nSubj,:) = resample(TMP_Dfc_sw(1,4,:),200,size(TMP_Dfc_sw,3));
    Whole_mPFC_Insular(nSubj,:) = resample(TMP_Dfc_sw(1,5,:),200,size(TMP_Dfc_sw,3));
    Whole_VTA_Insular(nSubj,:) = resample(TMP_Dfc_sw(4,5,:),200,size(TMP_Dfc_sw,3));
    
    
   cd(Dirdata); 
end

%% Select only child (just for plot)
Sub_child = BasicInfo.AgeAtScan > 6 & BasicInfo.AgeAtScan <= 12;
Sub_Con = strcmp(BasicInfo.SubjectType,'CONTROL');
Sub_Pat = strcmp(BasicInfo.SubjectType,'PATIENT');

Sub_child_Con = Sub_child.*Sub_Con;
Sub_child_Pat = Sub_child.*Sub_Pat;

Whole_mPFC_Amyg_Con = Whole_mPFC_Amyg(Sub_child_Con==1,:);
Whole_mPFC_Amyg_Pat = Whole_mPFC_Amyg(Sub_child_Pat==1,:);
Whole_mPFC_VTA_Con = Whole_mPFC_VTA(Sub_child_Con==1,:);
Whole_mPFC_VTA_Pat = Whole_mPFC_VTA(Sub_child_Pat==1,:);
Whole_mPFC_Insular_Con = Whole_mPFC_Insular(Sub_child_Con==1,:);
Whole_mPFC_Insular_Pat = Whole_mPFC_Insular(Sub_child_Pat==1,:);
Whole_VTA_Insular_Con = Whole_VTA_Insular(Sub_child_Con==1,:);
Whole_VTA_Insular_Pat = Whole_VTA_Insular(Sub_child_Pat==1,:);

%% Plotting
plot(Whole_mPFC_Amyg_Con','-k');
hold on; plot(Whole_mPFC_Amyg_Pat','-r');
xlim(gca,[1 200]);
ylabel(gca,'Dynamic Connectivity Strength','FontSize', 15, 'FontName', 'Calibri');
xlabel(gca,'Time(resampled)','FontSize', 15, 'FontName', 'Calibri');
set(gcf,'Position',[10 10 500 300]);

plot(Whole_mPFC_VTA_Con','-k');
hold on; plot(Whole_mPFC_VTA_Pat','-r');
xlim(gca,[1 200]);
ylabel(gca,'Dynamic Connectivity Strength','FontSize', 15, 'FontName', 'Calibri');
xlabel(gca,'Time(resampled)','FontSize', 15, 'FontName', 'Calibri');
set(gcf,'Position',[10 10 500 300]);

plot(Whole_mPFC_Insular_Con','-k');
hold on; plot(Whole_mPFC_Insular_Pat','-r');
xlim(gca,[1 200]);
ylabel(gca,'Dynamic Connectivity Strength','FontSize', 15, 'FontName', 'Calibri');
xlabel(gca,'Time(resampled)','FontSize', 15, 'FontName', 'Calibri');
set(gcf,'Position',[10 10 500 300]);

plot(Whole_VTA_Insular_Con','-k');
hold on; plot(Whole_VTA_Insular_Pat','-r');
xlim(gca,[1 200]);
ylabel(gca,'Dynamic Connectivity Strength','FontSize', 15, 'FontName', 'Calibri');
xlabel(gca,'Time(resampled)','FontSize', 15, 'FontName', 'Calibri');
set(gcf,'Position',[10 10 500 300]);