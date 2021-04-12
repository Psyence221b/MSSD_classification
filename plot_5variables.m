function [TMP] = plot_5variables(Table,option1,option2)
% option 1 => 'Raw' or 'Dfc'
% option 2 => variable name in the table ex) 'mPFC', 'Dfc_BLA_50TR'

List_Con = Table([Table.Sub_group == 1],:);
List_Pat = Table([Table.Sub_group == 2],:);

TMP_index = {'_mean', '_MSSD1', '_MSSD2', '_SD', '_VSD'};
for i = 1:length(TMP_index)
    List_index{i} = [option2 TMP_index{i}];
end

x_con = ones([1 size(List_Con,1)])*2;
x_pat = ones([1 size(List_Pat,1)])*3;

%% Plotting
if strcmp(option1, 'Raw')
% Raw values
figure;
eval(['scatter(x_con,List_Con.' List_index{1} '(:));']);
xlim(gca, [1 4]); 
title(gca, 'Average of Raw signal', 'FontSize', 15, 'FontName', 'Calibri');
ylabel(gca,'Ave(raw)','FontSize', 15, 'FontName', 'Calibri');
xlabel(gca,'Group','FontSize', 15, 'FontName', 'Calibri');
hold on;
eval(['scatter(x_pat,List_Pat.' List_index{1} '(:));']);
eval(['mean_con = mean(List_Con.' List_index{1} ');']); 
eval(['mean_pat = mean(List_Pat.' List_index{1} ');']);
hline_con = refline([0 mean_con]); hline_con.Color = 'b';
hline_pat = refline([0 mean_pat]); hline_pat.Color = 'r';
legend('Control','Patient','Ave(Control)','Ave(Patient)');
saveas(gcf,'AveRaw.png');

elseif strcmp(option1,'Dfc')
% Normalized (Fisher) Dfc
figure;
eval(['scatter(x_con,List_Con.' List_index{1} '(:));']);
xlim(gca, [1 4]); 
title(gca, 'Average of Dfc', 'FontSize', 15, 'FontName', 'Calibri');
ylabel(gca,'Ave(Dfc) (pearson-r based)','FontSize', 15, 'FontName', 'Calibri');
xlabel(gca,'Group','FontSize', 15, 'FontName', 'Calibri');
hold on;
eval(['scatter(x_pat,List_Pat.' List_index{1} '(:));']);
eval(['mean_con = mean(List_Con.' List_index{1} ');']); 
eval(['mean_pat = mean(List_Pat.' List_index{1} ');']);
hline_con = refline([0 mean_con]); hline_con.Color = 'b';
hline_pat = refline([0 mean_pat]); hline_pat.Color = 'r';
legend('Control','Patient','Ave(Control)','Ave(Patient)');
saveas(gcf,'AveDfc.png');
end

% MSSD1
figure;
eval(['scatter(x_con,List_Con.' List_index{2} '(:));']);
xlim(gca, [1 4]); 
title(gca, 'MSSD1(manuscript): mPFC and BLA', 'FontSize', 15, 'FontName', 'Calibri');
ylabel(gca,'MSSD1 (pearson-r based)','FontSize', 15, 'FontName', 'Calibri');
xlabel(gca,'Group','FontSize', 15, 'FontName', 'Calibri');
hold on;
eval(['scatter(x_pat,List_Pat.' List_index{2} '(:));']);
eval(['mean_con = mean(List_Con.' List_index{2} ');']); 
eval(['mean_pat = mean(List_Pat.' List_index{2} ');']);
hline_con = refline([0 mean_con]); hline_con.Color = 'b';
hline_pat = refline([0 mean_pat]); hline_pat.Color = 'r';
legend('Control','Patient','Ave(Control)','Ave(Patient)');
saveas(gcf,'MSSD1.png');

% MSSD2
figure;
eval(['scatter(x_con,List_Con.' List_index{3} '(:));']);
xlim(gca, [1 4]); 
title(gca, 'MSSD2(article): mPFC and BLA', 'FontSize', 15, 'FontName', 'Calibri');
ylabel(gca,'MSSD2 (pearson-r based)','FontSize', 15, 'FontName', 'Calibri');
xlabel(gca,'Group','FontSize', 15, 'FontName', 'Calibri');
hold on;
eval(['scatter(x_pat,List_Pat.' List_index{3} '(:));']);
eval(['mean_con = mean(List_Con.' List_index{3} ');']); 
eval(['mean_pat = mean(List_Pat.' List_index{3} ');']);
hline_con = refline([0 mean_con]); hline_con.Color = 'b';
hline_pat = refline([0 mean_pat]); hline_pat.Color = 'r';
legend('Control','Patient','Ave(Control)','Ave(Patient)');
saveas(gcf,'MSSD2.png');

% SD
figure;
eval(['scatter(x_con,List_Con.' List_index{4} '(:));']);
xlim(gca, [1 4]);  
title(gca, 'SD', 'FontSize', 15, 'FontName', 'Calibri');
ylabel(gca,'SD (pearson-r based)','FontSize', 15, 'FontName', 'Calibri');
xlabel(gca,'Group','FontSize', 15, 'FontName', 'Calibri');
hold on;
eval(['scatter(x_pat,List_Pat.' List_index{4} '(:));']);
eval(['mean_con = mean(List_Con.' List_index{4} ');']); 
eval(['mean_pat = mean(List_Pat.' List_index{4} ');']);
hline_con = refline([0 mean_con]); hline_con.Color = 'b';
hline_pat = refline([0 mean_pat]); hline_pat.Color = 'r';
legend('Control','Patient','Ave(Control)','Ave(Patient)');
saveas(gcf,'SD.png');

% VSD
figure;
eval(['scatter(x_con,List_Con.' List_index{5} '(:));']);
xlim(gca, [1 4]);  
title(gca, 'VSD(manuscript): mPFC and BLA', 'FontSize', 15, 'FontName', 'Calibri');
ylabel(gca,'VSD (pearson-r based)','FontSize', 15, 'FontName', 'Calibri');
xlabel(gca,'Group','FontSize', 15, 'FontName', 'Calibri');
hold on;
eval(['scatter(x_pat,List_Pat.' List_index{5} '(:));']);
eval(['mean_con = mean(List_Con.' List_index{5} ');']); 
eval(['mean_pat = mean(List_Pat.' List_index{5} ');']);
hline_con = refline([0 mean_con]); hline_con.Color = 'b';
hline_pat = refline([0 mean_pat]); hline_pat.Color = 'r';
legend('Control','Patient','Ave(Control)','Ave(Patient)');
saveas(gcf,'VSD.png');

clear mean_con mean_pat hline_con hline_pat
close all; clear List_Con List_Pat

end

