%% statistical analysis of the lognormal parameters you gathered in the "logfit" script
% first step is to read the data into the script
% repeat / copy this section 

% change the folder to where you saved all the parameters and values 
% cd('C:\Users\Nina R��ler\Google Drive\PhD\Paper\12_2021_new figures\gof_sigler')


dd0 = dir('*.xlsx'); % this lists all your csv files in the folder 
 
gof_stats0 = {dd0.name}; % takes the file names from the data 
data0 = cell(numel(gof_stats0)); % makes a cell array with the number of elements in the folder
data0(:,1) = regexprep(gof_stats0, '.xlsx',''); % replaces the file names with the info of each file, filling up the cell arrays
 
%%%% important: check the csv file and look which column the wanted
%%%% information is; then change nc accordingly

% read the data of parameters that you need
for h = 1:numel(gof_stats0)    
    data0{h} = xlsread(gof_stats0{h}); % read and load in all the files   
%     gofs0{h} = ([data0{h}]); % if you want to read every row, use this
%     line
    gofs0{h} = ([data0{h}(:,2)]); % specify which row you'd need, default is second row   
     
end

% specify now which row in the data belongs to which 
param07g0 = gofs0{1}; % the number here specifies which row in the document it is
param07g1 = gofs0{2};
param14g0 = gofs0{3};
param14g1 = gofs0{4};
param21g0 = gofs0{5};
param21g1 = gofs0{6}; 

% copy these if necessary and there are more than one parameter you want to
% check (change name of the variable accordingly)

%% statistical comparisons between different parameters / conditions / time points
% statistical comparison, doing a simple pairwise comparison between two
% conditions at different ages; or only use one comparison 

[rs07, h07] = ranksum(param07g0, param07g1);
ph07 = [rs07, h07];
xlswrite('Ranksum_gof_07_mushroom.xlsx', ph07)

[rs14, h14] = ranksum(param14g0, param14g1);
ph14 = [rs14, h14];
xlswrite('Ranksum_gof_14_mushroom.xlsx', ph14)

[rs21, h21] = ranksum(param21g0, param21g1);
ph21 = [rs21, h21];
xlswrite('Ranksum_sigma_21_mushroom.xlsx', ph21)

% in case of multiple comparisons (three conditions or time points): ranksum 
% for pairwise comparisons with correction for multiple tests

[pm, hm] = ranksum(param07g1, param14g1);
pmhm0714 = [pm, hm]; 
xlswrite('gof_comp_0714.xlsx',pmhm0714)

[pm1,hm1] = ranksum(param07g1, param21g1);
pmhm0721 = [pm1,hm1]; 
xlswrite('gof_comp_0721.xlsx', pmhm0721)

[pm2, hm2] = ranksum(param14g1, param21g1); 
pmhm1421 = [pm2,hm2]; 
xlswrite('gof_comp_1421.xlsx', pmhm1421)

pval = [pm, pm1, pm2];
[cor_p, ha] = bonf_holm(pval, 0.05);
cor_pn = array2table(cor_p); 
writetable(cor_pn, 'gof_corrected_p_g1.xlsx')

%% visually compare statistical parameters: compare two parameters between ages or conditions

% define the groups / conditions you want to compare
me0 = param21g0; 
me1 = param21g1; 


si = [(mean(me0)),(mean(me1))];
sir = [((std(me0))/(sqrt(length(me0)))), ((std(me1))/(sqrt(length(me1))))];

x = [1,2];

plot(x(1), me0, 'r.')
hold on 
errorbar(x(1), si(1),sir(1),'k');
plot(x(2),me1, 'b.')
errorbar(x(2), si(2),sir(2),'k');

xlim([0 3])
ylim([0 1])

pagesize = [3 2]; %(x, y) size in cm

legend('boxoff'); 
legend('off');

set (gca, ...
    'ActivePositionProperty', 'position', ...
    'position', [0.2 0.2 0.75 0.75], ...
    'ytick', 0 : 0.1 : 1, ...
    'xtick', 0 : 1 : 3, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ticklength', [0.096 0.24] ./ max (pagesize), ...
    'tickdir', 'out', ...
    'linewidth', 0.5, ...
    'fontsize', 6, ...
    'fontname', 'arial', ... 
    'box', 'off');
tprint ( ...
    './sigma_comp_21_alt', ...
    '-HR -jpg -eps', ...
    pagesize);


%% look at graphical and statistical comparisons of your parameters
% starting with group comparison 
cd('C:\Users\Nina R��ler\Google Drive\PhD\Paper\12_2021_new figures')

% define which parameters you want to compare
group1 = param07g0; 
group2 = param07g1;

xm = [1,2];

ym = [(mean(group1)), (mean(group2))]; 
errm = [(std(group1) / sqrt(length(group1))), (std(group2) / sqrt(length(group2)))]; 

plot(xm(1), group1, 'r.')
hold on 
plot(xm(2), group2, 'b.')
errorbar(ym, errm, 'k')
xlim([0 3])
ylim([-1 4])
               
pagesize = [3 2]; %(x, y) size in cm

legend('boxoff'); 
legend('off');

set (gca, ...
    'ActivePositionProperty', 'position', ...
    'position', [0.2 0.2 0.75 0.75], ...
    'ytick', -1 : 1 : 4, ...
    'xtick', 0 : 1 : 3, ...
    'xticklabel', {'', '','',''}, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ticklength', [0.096 0.24] ./ max (pagesize), ...
    'tickdir', 'out', ...
    'linewidth', 0.5, ...
    'fontsize', 6, ...
    'fontname', 'arial', ... 
    'box', 'off');
tprint ( ...
    './skew_comp_group_all times', ...
    '-HR -jpg -eps', ...
    pagesize);

% make a comparison between groups for each time point - plot all
% comparisons in one graph
% first conidition
ym0_n = [(mean(param07g0)), (mean(param14g0)), (mean(param21g0))]; 
erm0_n = [((std(param07g0))/(sqrt(length(param07g0)))), ((std(param14g0))/(sqrt(length(param14g0)))), ...
    ((std(param21g0))/(sqrt(length(param21g0))))]; 

% second condition
ym1_n = [(mean(param07g1)), (mean(param14g1)), (mean(param21g1))]; 
erm1_n = [((std(param07g1))/(sqrt(length(param07g1)))), ((std(param14g1))/(sqrt(length(param14g1)))), ...
    ((std(param21g1))/(sqrt(length(param21g1))))]; 

figure
x = [1 2 3 4 5 6 7 8 9];

figure3 = plot(x(1), param07g0,'r.');
hold on
errorbar(x(1),ym0_n(1), erm0_n(1), 'k')
plot(x(2),param07g1, '.b')
errorbar(x(2),ym1_n(1), erm1_n(1), 'k')
plot(x(4),param14g0, 'r.')
errorbar(x(4),ym0_n(2), erm0_n(2), 'k')
plot(x(5),param14g1, 'b.') 
errorbar(x(5),ym1_n(2), erm1_n(2), 'k')
plot(x(7),param21g0, 'r.')
errorbar(x(7),ym0_n(3), erm0_n(3), 'k')
plot(x(8),param21g1, 'b.')
errorbar(x(8),ym1_n(3), erm1_n(3), 'k')

% figure2 = errorbar(yy_ipsi,err_ipsi,'k', 'LineWidth',2);

xlim([0 9]);
% ylim([-1 4]);
ylim([0 1.2]);               
pagesize = [3 2]; %(x, y) size in cm

legend('boxoff'); 
legend('off');

set (gca, ...
    'ActivePositionProperty', 'position', ...
    'position', [0.2 0.2 0.75 0.75], ...
    'ytick', 0 : 0.2 : 1.2, ...
    'xtick', 0 : 3 : 9, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'xticklabels', {'','','',''}, ...
    'ticklength', [0.096 0.24] ./ max (pagesize), ...
    'tickdir', 'out', ...
    'linewidth', 0.5, ...
    'fontsize', 6, ...
    'fontname', 'arial', ... 
    'box', 'off');
tprint ( ...
    './sigma_comp_groups and times_allsp', ...
    '-HR -jpg -eps', ...
    pagesize);


%% if applicable: take a closer look at the differences at different time points
% first for condition 1 

ym0 = [(mean(param07g0)), (mean(param14g0)), (mean(param21g0))]; 
erm0 = [((std(param07g0))/(sqrt(length(param07g0)))), ((std(param14g0))/(sqrt(length(param14g0)))), ...
    ((std(param21g0))/(sqrt(length(param21g0))))]; 

% now compare the spines in different times
figure
x = [1,2,3]; 
gmush_time = plot(x(1), param07g0, 'r.'); 
hold on
plot(x(2), param14g0, 'r.')
plot(x(3), param21g0, 'r.')
errorbar(ym0, erm0, 'k')
hold off

xlim([0 4]);
ylim([-1 4]);
               
pagesize = [3 2]; %(x, y) size in cm

legend('boxoff'); 
legend('off');

set (gca, ...
    'ActivePositionProperty', 'position', ...
    'position', [0.2 0.2 0.75 0.75], ...
    'ytick', -1 : 1 : 4, ...
    'xtick', 0 : 1 : 4, ...
    'xticklabels', {'' '' '' ''}, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ticklength', [0.096 0.24] ./ max (pagesize), ...
    'tickdir', 'out', ...
    'linewidth', 0.5, ...
    'fontsize', 6, ...
    'fontname', 'arial', ... 
    'box', 'off');
tprint ( ...
    './skew_comp_g0_all times', ...
    '-HR -jpg -eps', ...
    pagesize);

% next step is looking at time differences in condition 2
ym1 = [(mean(param07g1)), (mean(param14g1)), (mean(param21g1))]; 
erm1 = [((std(param07g1))/(sqrt(length(param07g1)))), ((std(param14g1))/(sqrt(length(param14g1)))), ...
    ((std(param21g1))/(sqrt(length(param21g1))))]; 

% now compare the spines in different times
figure
x = [1,2,3]; 
gmush_time1 = plot(x(1), param07g1, 'b.'); 
hold on
plot(x(2), param14g1, 'b.')
plot(x(3), param21g1, 'b.')
errorbar(ym1, erm1, 'k')
hold off

xlim([0 4]);
ylim([-1 4]);
               
pagesize = [3 2]; %(x, y) size in cm

legend('boxoff'); 
legend('off');

set (gca, ...
    'ActivePositionProperty', 'position', ...
    'position', [0.2 0.2 0.75 0.75], ...
    'ytick', -1 : 1 : 4, ...
    'xtick', 0 : 1 : 4, ...
    'xticklabels', {'' '' '' ''}, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ticklength', [0.096 0.24] ./ max (pagesize), ...
    'tickdir', 'out', ...
    'linewidth', 0.5, ...
    'fontsize', 6, ...
    'fontname', 'arial', ... 
    'box', 'off');
tprint ( ...
    './skew_comp_g1_all times', ...
    '-HR -jpg -eps', ...
    pagesize);

%% aic analysis for each condition and time point seperately

cd('C:\Users\Nina R��ler\Google Drive\PhD\Paper\12_2021_new figures')

% define which condition and time point you want to compare
aicnew = param07g0s'; 

figure
for j = 1:length(aicnew)
    
    y = [(mean(aicnew(j,1))),(mean(aicnew(j,2))),(mean(aicnew(j,3)))];
    er = [((std(aicnew(j,1)))/(sqrt(length(aicnew(j,1))))), ((std(aicnew(j,2)))/(sqrt(length(aicnew(j,2))))), ...
        ((std(aicnew(j,3)))/(sqrt(length(aicnew(j,3)))))];
    
    x = [1 2 3];
    
    f1 = plot(x(1), aicnew(j,1), '.r');
    hold on 
    plot(x(2), aicnew(j,2), '.r')
    plot(x(3), aicnew(j,3), '.r')
    f2 = errorbar(y, er, 'r'); 
    
    xlim([0 4])
    
end 

yy = [(mean(aicnew(:,1))),(mean(aicnew(:,2))),(mean(aicnew(:,3)))];
err = [((std(aicnew(:,1)))/(sqrt(length(aicnew(:,1))))), ((std(aicnew(:,2)))/(sqrt(length(aicnew(:,2))))), ...
    ((std(aicnew(:,3)))/(sqrt(length(aicnew(:,3)))))];


f3 = errorbar(yy,err,'k');
hold off
xlim([0 4]);
% ylim([-1 4]);
               
pagesize = [3 2]; %(x, y) size in cm

legend('boxoff'); 
legend('off');

set (gca, ...
    'ActivePositionProperty', 'position', ...
    'position', [0.2 0.2 0.75 0.75], ...
    'xtick', 0 : 1 : 4, ...
    'xticklabels', {'' '' '' ''}, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ticklength', [0.096 0.24] ./ max (pagesize), ...
    'tickdir', 'out', ...
    'linewidth', 0.5, ...
    'fontsize', 6, ...
    'fontname', 'arial', ... 
    'box', 'off');
tprint ( ...
    './aic_comp_d07g0_stub', ...
    '-HR -jpg -eps', ...
    pagesize);



