cd(Dirlog);
MRIInfo = readtable('BasicInfo_onlyRunData.csv');
% Change wrong variable name
MRIInfo.Properties.VariableNames{'SamlpingRate'} = 'SamplingRate';

FDInfo = readtable('final_FDmotion_before.csv');

cd(Dirdocu);
score = readtable('Andlab_ABIDE_summary.xlsx');
score(1,:) = [];

% Trim SamplingRate
List06 = MRIInfo.SamplingRate >= 0.59 & MRIInfo.SamplingRate <= 0.61;
MRIInfo.SamplingRate(List06) = 0.6;
List05 = MRIInfo.SamplingRate >= 0.49 & MRIInfo.SamplingRate <= 0.55;
MRIInfo.SamplingRate(List05) = 0.5;

% Remove Last row
MRIInfo(end,:) = [];

%% Load demographic information
cd(Dirdocu);
DemoInfo = readtable('Andlab_ABIDE_basic_subject_info.csv');

TMP_table = table();
for nTable1 = 1:size(MRIInfo,1)
    for nTable2 = 1:size(DemoInfo,1)
        if strcmp(MRIInfo.subID(nTable1), DemoInfo.x___Lab_ID(nTable2))
            TMP_table(nTable1,:) = DemoInfo(nTable2,:);
        end
    end
end
TMP_table = removevars(TMP_table, 'x___Lab_ID');

BasicInfo = [MRIInfo TMP_table FDInfo(:,2)];
BasicInfo.SubjectType{131} = 'CONTROL';
BasicInfo.Properties.VariableNames{'Var2'} = 'FD';

%% Calculate quartile of age to choose a criterion of age group
Age_quartile = quantile(DemoInfo.AgeAtScan,3);

%% Remove subjects
% Remove subjects if SamplingRate is too low
ListDelete = BasicInfo.VolumeN <= 100;
BasicInfo(ListDelete,:) = [];

% Remove subjects if other information is missing
ListDelete = BasicInfo.Sex ~= 1 & BasicInfo.Sex ~= 2;
BasicInfo(ListDelete,:) = [];

% Remove subjects if he/she has lots of movement
ListDelete = BasicInfo.FD >= 0.25;
BasicInfo(ListDelete,:) = [];


%% Add ADOS, SRS, SQT scores
% find subject indices corresponding to score variable
Idx_subj = [];
parfor nSubj = 1:size(BasicInfo,1)
    Idx_subj(nSubj) = find(strcmp(BasicInfo.subID(nSubj),score.Var3));    
end

L_ADOS = score(Idx_subj,[34:36,44:46,66:69]);
for i = 1:size(L_ADOS,2), L_ADOS.(i) = str2double(L_ADOS{:,i}); end % change to number instead of character

L_ADI = score(Idx_subj,[37:41]);
for i = 1:size(L_ADI,2), L_ADI.(i) = str2double(L_ADI{:,i}); end

L_SRS = score(Idx_subj,[48:49,73:77]);
for i = 1:size(L_SRS,2), L_SRS.(i) = str2double(L_SRS{:,i}); end

L_VABS = score(Idx_subj,[51:65]);
for i = 1:size(L_VABS,2), L_VABS.(i) = str2double(L_VABS{:,i}); end

% remove -9999 values
Table_score = [L_ADOS L_ADI L_SRS L_VABS];
for i = 1:size(Table_score,2)
    Idx_NaN = find(Table_score{:,i} == -9999);
    Table_score(Idx_NaN,i) = {NaN};
end

% Add to the BasicInfo table

BasicInfo = [BasicInfo Table_score];

%%
clear DemoInfo MRIInfo TMP_table nTable1 nTable2 List05 List06 FDInfo i Idx_subj L_ADOS L_ADI L_SRS L_VABS Table_score
