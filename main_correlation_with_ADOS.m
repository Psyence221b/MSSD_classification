Data = ResultTable_TR50sw_agegroup1;
Data = Data([1:72],:); % 1 = CON
Data = Data([73:end],:); % 2 = PAT


% Array_connect = Data.mPFC_Amyg;
Array_connect = Data.TR50sw_mPFC_Amyg;
Array_connect = Array_connect(:,2);

Array_ADOS_Total = Data.ADOS_Total;
Array_ADOS_Comm = Data.ADOS_Comm;
Array_ADOS_Social = Data.ADOS_Social;
Array_ADOS_StBeh = Data.ADOS_StBeh;
Array_ADOS_SocAffect = Data.ADOSGotham_SocAffect;
Array_ADOS_Rep_RepBeh = Data.ADOSGotham_Rest_RepBeh;
Array_ADOS_Gotham_Total = Data.ADOSGotham_Total;
Array_ADOS_Gotham_Severity = Data.ADOSGotham_Severity;
Array_ADI_Social_Total = Data.ADI_R_SocialTotal;
Array_ADI_Verbal_Total = Data.ADI_R_VerbalTotal;
Array_ADI_RRB = Data.ADI_R_RRB;
Array_SRS = Data.SRS_RawTotal;


%%
x = Array_connect;
y = Array_SRS;
Idx_non_y = ~isnan(y);
x = x(Idx_non_y);
y = y(Idx_non_y);

[rob,~,rob_corrw] = andlab_robustfit( x, y);
r = rob_corrw
p = rob.stats.p(2)
n = length(x)