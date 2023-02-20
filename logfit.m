%% testing for lognormal fit
% this script fits your data visually to a lognormal distribution, and
% gathers statistical parameters that will be used for a statistical
% analysis in a different script 
% for the statisticial analysis, refer to "stat_analysis_lognormal_data"
 
% clear all variables and spaces before starting the analysis
clc
clear all
% close all

%  for PC in Gie�en 
% cd('C:\Users\Nina R��ler\Documents\Master project\Spine data analysis\spine-size_information-storage_20180604\new_dpi21\IML\Ipsi');
cd('G:\My Drive\PhD\Data Analysis\TNF data\TNF alpha\SP-\KO')
% load data into matlab
dd = dir('*.xlsx'); % this lists all your csv files in the folder 
 
spineposition = {dd.name}; % takes the file names from the data 
data = cell(numel(spineposition)); % makes a cell array with the number of elements in the folder
data(:,1) = regexprep(spineposition, '.xlsx',''); % replaces the file names with the info of each file, filling up the cell arrays
 

%%%% important: check the csv file and look which column the wanted
%%%% information is; then change nc accordingly

% nc is number of the row of your data, in my example it is row number 2; for dpi21 it is 6;
% please change accordingly 
% nt is number of row in which you want to start (sometimes text in the
% first rows is not accepted, so you have to cut it off; 21-21 is 0, 21 is
% 5
currfold = pwd;

% if currfold(111:115) == 'dpi21'
%     nt = 5; 
%     nc = 6;
% else
%     nt = 0;
%     nc = 2; 
% end
spineposition2 = spineposition(~contains(spineposition,'$'));
% read the data, and only take e.g. spinesize out of the data 

for h = 1:numel(spineposition2)    
%     data{h} = csvread(spineposition{h},1); % read and load in all the files   
    data{h} = readcell(spineposition2{h}, 'Range',2);
%     data2{h} = string(data{h}(:,2));
    spinearea{h} = ([data{h,1}]);
%     for m = 1:length(data2{h})
%         if data2{h}(m) == 'ipsi'
%            spinearea{h}(m) = ([data{h,1}]); % only include the area row    
%         end
%     end
end

% create matrix of spine sizes for easier further analysis
spinea = nonzeros(cell2mat(spinearea{1}'));
% spinea = spinea*100; 

% calculate the skewness of the data 
skew_data = skewness(spinea); 
kurt_data = kurtosis(spinea); 

%% fit the data to a lognormal distribution
% change folder to your analysis destination
% cd('C:\Users\Nina R��ler\Google Drive\PhD\Paper\12_2021_new figures')
cd('G:\My Drive\PhD\Data Analysis\TNF data\Auswertung')

% fit the data with a lognormal fit
% t = sprintf('Logfit IML21dpi ipsilateral'); % the title of your graph - not necessary
binsize = 0.05; % size of bins data will be grouped
mulow = -2;
silow = 0;
muup = 2;
siup = 2; 
[log_plot, gof_stats, coeff, fiti, xa, yaf] = logfitting(spinea, binsize, mulow, silow, muup, siup);

% save the figure of the fit and the statistical values (goodness of fit (gof) and
% parameters (mu and sigma))
% savefig('name of figure.fig'); % not necessary
writetable((struct2table(gof_stats)), '23_02_14_gof_log_TNF alpha KO SP-.xlsx') %%% you can change the name of the table to whatever you prefer) 
writetable((table(coeff)), '23_02_14_param_log_TNF alpha KO SP-.xlsx') %%% name can be changed as well

% collect the skewnesses in a table 
writetable((table(skew_data)), '23_02_14_skew_comp_TNF alpha KO SP-.xlsx')

% as another goodness of fit measure, calculate the aic of the different
% distributions to compare and see, which distributions best fits the data 

% [aic_l, aic_g, aic_w] = aiccomp(spinea);
% writetable(table(aic_g, aic_l, aic_w, 'VariableNames', {'AIC_Gamma', 'AIC_Lognormal', 'AIC_Weibull'}), '23_01_12_aic_comp_IMLday0lesion.xlsx');


%% if applicable: control data or any data point you want to compare against the first one
% change folder to data storage
% cd('C:\Users\Nina R��ler\Documents\Master project\Spine data analysis\spine-size_information-storage_20180604\new_dpi21\IML\Contra');
cd('G:\My Drive\PhD\Data Analysis\TNF data\TNF alpha\SP-\Control')

% load data into matlab
dd_c = dir('*.xlsx'); % this lists all your csv files in the folder 

spineposition_c = {dd_c.name}; % takes the file names from the data 
data_c = cell(numel(spineposition_c)); % makes a cell array with the number of elements in the folder
data_c(:,1) = regexprep(spineposition_c, '.xlsx',''); % replaces the file names with the info of each file, filling up the cell arrays
 

%%%% important: check the csv file and look which column the wanted
%%%% information is; then change nc accordingly

% nc is number of the row of your data, in my example it is row number 2; for dpi21 it is 6;
% please change accordingly 
% nt is number of row in which you want to start (sometimes text in the
% first rows is not accepted, so you have to cut it off; 21-21 is 0, 21 is
% 5
spineposition_c2 = spineposition_c(~contains(spineposition_c,'$'));

% read the data, and only take e.g. spinesize out of the data 
for hc = 1:numel(spineposition_c2)    
%     data_c{hc} = csvread(spineposition_c{hc},1,nt); % read and load in all the files   
    data_c{hc} = readcell(spineposition_c2{hc},'Range', 2); 
%     spinearea_c{hc} = ([data_c{hc,1}]'); % only include the area row    
%     data_c2{hc} = string(data_c{hc}(:,2));
    
%     for mc = 1:length(data_c2{hc})
%         if data_c2{hc}(mc) == 'contra'
           spinearea_c{hc} = ([data_c{hc,1}]); % only include the area row    
%         end
%     end
end

% create matrix of spine sizes for easier further analysis
spinea_c = nonzeros(cell2mat(spinearea_c{1})');
% spinea = spinea*100; 

% calculate the skewness of the data 
skew_data_control = skewness(spinea_c); 
kurt_data_control = kurtosis(spinea_c); 

%% fit lognormal distribution to your data
% change folder again to where analysis is stored
% cd('C:\Users\Nina R��ler\Google Drive\PhD\Paper\12_2021_new figures')
cd('G:\My Drive\PhD\Data Analysis\TNF data\Auswertung')

% fit the data with a lognormal fit
% tc = sprintf('Logfit IML21dpi contralateral'); % not necessary
binsize_con = 0.05; 
[eac, log_plotc, gof_statsc, coeffc, fitc, xa_c, yaf_c] = logfitting(spinea_c, binsize_con, mulow, silow, muup, siup);

% if aaplicable: normalize the data to properly compare goodness of fit and
% skewness (as seen in Hazan & Ziv, 2019(?), figure ?)
% for normalization, calculate plasticity / control ratio of spine sizes
a = mean(spinea); 
b = mean(spinea_c); 
ratio = a/b; 

% now normalise all spine sizes 
new_spinea_c = spinea_c .* ratio; 
skew_cn = skewness(new_spinea_c); 

% tcn = sprintf('Normalized Logfit IML21dpi contralateral'); 
[log_plotc_n, gof_statsc_n, coeffc_n, fitc_n, xa_c_n, yaf_c_n] = logfitting(new_spinea_c, binsize);

% save the figure of the fit and the statistical values (gof and
% parameters)
% savefig('logfitting_IML21contra.fig');
writetable((struct2table(gof_statsc_n)), '23_02_14_norm_gof_log_TNF alpha control SP-.xlsx') %%% you can change the name of the table to whatever you prefer) 
writetable((struct2table(gof_statsc)), '23_02_14_gof_log_TNF alpha control SP-.xlsx') %%% you can change the name of the table to whatever you prefer) 
writetable((table(coeffc)), '23_02_14_param_log_TNF alpha control SP+.xlsx') %%% name can be changed as well

% save the skewness both normalised (if applicable) and not normalised
writetable(table(skew_data_control), '23_02_14_skew_comp_TNF alpha control SP-.xlsx')
% writetable(table(skew_cn, gof_statsc_n.rsquare), 'norm_skew_gof_IML21contra.xlsx')


% as another goodness of fit measure, calculate the aic of the different
% distributions to compare and see, which distributions best fits the data 
% [aic_lc, aic_gc, aic_wc] = aiccomp(spinea_c);
% writetable(table(aic_gc, aic_lc, aic_wc, 'VariableNames', {'AIC_Gamma', 'AIC_Lognormal', 'AIC_Weibull'}), 'aic_comp_IML21contra.xlsx');

%% generate comparison graphs: plot both fits / conditions in one graph
     
figure
figure3 = plot(fiti, '-g', xa, yaf, '-g');
set(figure3, 'LineWidth', 0.2);
hold on 
plot(fitc, '-b', xa_c, yaf_c, '-b'); 
set(figure3, 'LineWidth', 0.2); 
hold off

xlim([0 0.75]);
ylim([0 6]);        
pagesize = [3 2]; %(x, y) size in cm

savefig('23_02_14_logfit comparison_TNF alpha SP-.fig')
% you can either add the legend manually later or add it here
% legend('','ipsilateral', '', 'contralateral', ... 
%     'box', 'off', ...
%     'location', 'northeast');

legend('boxoff'); 
legend('off');

% check all parameters on the graph, see that axes and linewidth etc. are
% set to what you need
set (gca, ...
    'ActivePositionProperty', 'position', ...
    'position', [0.2 0.2 0.75 0.75], ...
    'ytick', 0 : 1 : 6, ...
    'xtick', 0 : 0.25 : 0.75, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ticklength', [0.096 0.24] ./ max (pagesize), ...
    'tickdir', 'out', ...
    'linewidth', 0.5, ...
    'fontsize', 6, ...
    'fontname', 'arial', ... 
    'box', 'off');
tprint ( ...
    './23_02_14_logfit_comp_TNF alpha SP-', ...
    '-HR -jpg -eps', ...
    pagesize);

%% cut off x axis at 0.5 to see differences in the peak a bit better (if
% applicable)
figure
figure3 = plot(fiti, '-g', xa, yaf, '-g');
set(figure3, 'LineWidth', 0.2);
hold on

figure3 = plot(fitc, '-b', xa_c, yaf_c, '-b'); 
set(figure3, 'LineWidth', 0.2); 
hold off
    

xlim([0 0.5]);
ylim([0 15]);        
pagesize = [3 2]; %(x, y) size in cm

legend('boxoff'); 
legend('off');

set (gca, ...
    'ActivePositionProperty', 'position', ...
    'position', [0.2 0.2 0.75 0.75], ...
    'ytick', 0 : 5 : 15, ...
    'xtick', 0 : 0.025 : 0.5, ...
    'xticklabel', {'0' '' '' '' '' '' '' '' '' '' '0.25' '' '' '' '' '' '' '' '' '' '0.5'}, ... 
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ticklength', [0.096 0.24] ./ max (pagesize), ...
    'tickdir', 'out', ...
    'linewidth', 0.5, ...
    'fontsize', 6, ...
    'fontname', 'arial', ... 
    'box', 'off');
tprint ( ...
    './logfit_comp_IML21_cutoff', ...
    '-HR -jpg -eps', ...
    pagesize);

%% look at lognbormal fits for each individual cell in each condition

% cd('C:\Users\Nina R��ler\Google Drive\PhD\Paper\02_2022_new figures')
% fit inidiviual cells and plot each individually
for k = 1:numel(spinearea_ind) 
    
    ske{k} = skewness(spinearea_ind{k});
%     [aic_l{k}, aic_g{k}, aic_w{k}] = aiccomp(nonzeros(spinearea{k}));

    % print out the lognormal fit for each individual cell and put out
    % plots for each
%     ti = sprintf('Logfit IML21dpi ipsilateral Single Cell'); 
    [l_plot{k}, gof_s{k}, co{k}, lfit{k}, x_a{k}, y_af{k}] = logfitting(spinearea_ind{k}, binsize, mulow, silow, muup, siup); 
    
    close 
    
end

% save all the needed parameters in tables in the folder of your analysis
% goodness of fit statistics
coef_tab_indcells = cell2mat(gof_s); 
writetable(struct2table(coef_tab_indcells), '23_01_31_gof_stats_indcells_AISipsi_synpo cluster density.xlsx');
% aic comparison
% writetable(table(aic_g', aic_l', aic_w', 'VariableNames', {'AIC_Gamma', 'AIC_Lognormal', 'AIC_Weibull'}), 'aic_comp_IML21ipsi.xlsx');
% skewness
% writetable(table(ske), 'skew_indcells_ipsiIML21.xlsx')
% parameters mu and sigma
param_tab_indcell = cell2mat(co'); 
writetable(table(param_tab_indcell), '23_01_31_param_indcells_AISipsi_synpo cluster density.xlsx')


% put all individual cells plus the average in the same graph 
figure()
for i = 1:numel(spinearea_ind)
    % bin the spine sizes together, you can change bin size
    ni{i} = 0:0.001:max(spinearea_ind{i});
    % take parameters mu and sigma from the individual fits above and feed
    % them into these fits 
    mi{i} = lfit{i}.mu_all;
    si{i} = lfit{i}.s_all;
    % pdf of lognormal distribution to fit to data
    y_l{i} = (1./(ni{i}.*si{i}.*sqrt(2.*pi))) .* exp(-((log(ni{i})-mi{i}).^2) ./ (2.*si{i}.^2));
    % plot each individual fit into the same graph
    figure26 = plot(ni{i}, y_l{i}, 'g');
    figure26.LineWidth = 0.2;
    
    hold on
    
end 

% plot the average that was calculated as a first step into the individual
% plots as well
p1 = plot(fiti,'r');
p1(1).LineWidth = 0.2; 

xlim([0 0.4]);
ylim([0 16]);        
pagesize = [3 2]; %(x, y) size in cm

legend('boxoff'); 
legend('off');

savefig('23_01_31_indcells logfit AISipsi_synpo cluster density.fig')
set (gca, ...
    'ActivePositionProperty', 'position', ...
    'position', [0.2 0.2 0.75 0.75], ...
    'ytick', 0 : 2 : 16, ...
    'xtick', 0 : 0.1 : 0.4, ...
    'xticklabel', {'0' '' '' '' '0.4'}, ... 
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ticklength', [0.096 0.24] ./ max (pagesize), ...
    'tickdir', 'out', ...
    'linewidth', 0.5, ...
    'fontsize', 6, ...
    'fontname', 'arial', ... 
    'box', 'off');
tprint ( ...
    './23_01_18_logfit_allcells_comp_AISipsi_synpo cluster density', ...
    '-HR -jpg -eps', ...
    pagesize);

    
%% look at lognormal fits for each individual cell - control group (if applicable) 

for kc = 1:numel(spinearea_indc) 
        
    skec{kc} = skewness(spinearea_indc{kc});
%     [aic_lc{kc}, aic_gc{kc}, aic_wc{kc}] = aiccomp(nonzeros(spinearea_c{kc}));
    
    % print out the lognormal fit for each individual cell and put out
    % plots for each
%     t_c = sprintf('Logfit IML21dpi contralateral Single Cell'); 
    [l_plotc{kc}, gof_sc{kc}, coc{kc}, lfitc{kc}, x_ac{kc}, y_afc{kc}] = logfitting(spinearea_indc{kc}, binsize, mulow, silow, muup, siup); 
    
    close 
    
end

coef_tab_indcells_c = cell2mat(gof_sc); 
writetable(struct2table(coef_tab_indcells_c), '23_01_31_gof_stats_indcells_AIScontra_synpo cluster density.xlsx');

param_tab_indcell_c = cell2mat(coc');
writetable(table(param_tab_indcell_c), '23_01_31_param_indcells_AIScontra_synpo cluster density.xlsx')

% writetable(table(aic_gc', aic_lc', aic_wc', 'VariableNames', {'AIC_Gamma', 'AIC_Lognormal', 'AIC_Weibull'}), 'aic_comp_IML21contra.xlsx');
writetable(table(skec), '23_01_31_skew_indcells_AIScontra_synpo cluster density.xlsx')

figure
% put all individual cells and the average in the same graph 
for ic = 1:numel(spinearea_indc)
    
    nic{ic} = 0:0.001:max(spinearea_indc{ic});
    mci{ic} = lfitc{ic}.mu_all;
    sc{ic} = lfitc{ic}.s_all;
    
    y_lc{ic} = (1./(nic{ic}.*sc{ic}.*sqrt(2.*pi))) .* exp(-((log(nic{ic})-mci{ic}).^2) ./ (2.*sc{ic}.^2));
    
    figure27 = plot(nic{ic}, y_lc{ic}, 'b'); 
    figure27.LineWidth = 0.2; 
    
    hold on
end 
    
p2 = plot(fitc,'r');
p2(1).LineWidth = 0.2; 
    
xlim([0 0.4]);
ylim([0 16]);        
pagesize = [3 2]; %(x, y) size in cm

legend('boxoff'); 
legend('off');

savefig('23_01_31_indcells logfit AIScontra_synpo cluster density.fig')

set (gca, ...
    'ActivePositionProperty', 'position', ...
    'position', [0.2 0.2 0.75 0.75], ...
    'ytick', 0 : 2 : 16, ...
    'xtick', 0 : 0.1 : 0.4, ...
    'xticklabel', {'0' '' '' '' '0.4'}, ... 
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ticklength', [0.096 0.24] ./ max (pagesize), ...
    'tickdir', 'out', ...
    'linewidth', 0.5, ...
    'fontsize', 6, ...
    'fontname', 'arial', ... 
    'box', 'off');
tprint ( ...
    './23_01_31_logfit_allcells_comp_AIScontra_synpo cluster density', ...
    '-HR -jpg -eps', ...
    pagesize);


%% plotting the logarithm to check whether it is normal and thus a true lognormal distribution 

%take the logarithm of the data, set beginning and end points
norm_d = log10(spinea); 
start_analysis = min(norm_d);
end_analysis = max(norm_d); % the endpoint is your biggest data point 

% creating size categories or whatever else categories you need for the
% analysis
binsize = 0.1;
ed = start_analysis:binsize:end_analysis; % the first number is your starting point, the middle number is the steps / size of the bins
% the last number is the end point for the bins; 
%%% change your start and endpoints and steps in between accordingly

[num_s, eds] = histcounts(norm_d,ed); % this function gives you the number of bins with the specified parameter above
% and the number of data points in each bin

ed_f = ed(1:length(num_s));
       
xa = ed_f'; 
ya = num_s'; 
    
Aa = trapz(xa,ya); 
yaf = ya / Aa;
       
[nfit, ngof] = fit(xa, ya, 'gauss1');
    
writetable(struct2table(ngof), '23_02_14_gof_loglog_TNF alpha KO SP-.xlsx');
figure
plot(nfit, xa,ya)
 
norm_dc = log10(spinea_c); 
start_analysisc = min(norm_dc); 
end_analysisc = max(norm_dc); % the endpoint is your biggest data point 

% creating size categories or whatever else categories you need for the
% analysis
edc = start_analysisc:binsize:end_analysisc; % the first number is your starting point, the middle number is the steps / size of the bins
% the last number is the end point for the bins; 
%%% change your start and endpoints and steps in between accordingly

[num_sc, edsc] = histcounts(norm_dc,edc); % this function gives you the number of bins with the specified parameter above
% and the number of data points in each bin
   

ed_fc = edc(1:length(num_sc));
  
xac = ed_fc'; 
yac = num_sc'; 
    
Aac = trapz(xac,yac); 
yafc = yac / Aac;
[nfitc, ngofc] = fit(xac, yac, 'gauss1'); 
   
writetable(struct2table(ngofc), '23_02_14_gof_loglog_TNF alpha control SP-.xlsx');
    
%     plot(nfitc, xac,yac)
 
my = nfit(xa);
myc = nfitc(xac); 
figure
p2 = plot(nfit, '-g', xa,ya, 'g-'); 
set(p2, 'LineWidth', 0.2);
hold on
pa = patch( ...
     [xa; fliplr(xa);xa(1)], ...
     [ya; fliplr(my);my(1)], 'green', 'EdgeColor', 'none');
pa.FaceAlpha = 0.2;
    
%     hold on 
p3 = plot(nfitc, '-b', xac, yac, 'b-');
set(p3, 'LineWidth', 0.2); 
pe = patch(...
     [xac; fliplr(xac); xac(1)], ... 
     [max([yac,myc],[],2); fliplr(min([yac,myc],[],2));min([yac(1),myc(1)],[],2);], 'blue', 'EdgeColor', 'none'); 
pe.FaceAlpha = 0.2; 
    

ylim([0 260])
xlim([-2.5 0])        
pagesize = [3 2]; %(x, y) size in cm

legend('boxoff'); 
legend('off');

set (gca, ...
    'ActivePositionProperty', 'position', ...
    'position', [0.2 0.2 0.75 0.75], ...
    'ytick', 0 : 20 : 260, ...
    'xtick', -2.5 : 0.55 : 0, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ticklength', [0.096 0.24] ./ max (pagesize), ...
    'tickdir', 'out', ...
    'linewidth', 0.5, ...
    'fontsize', 6, ...
    'fontname', 'arial', ... 
    'box', 'off');
tprint ( ...
    './23_02_14_loglogpatch_comp_TNF alpha SP-', ...
    '-HR -jpg -eps', ...
    pagesize);

% loglogplot but without the patch

figure 
     
p2 = plot(nfit, '-g', xa,ya, 'g-'); 
set(p2, 'LineWidth', 0.2);
hold on
p3 = plot(nfitc, '-b', xac, yac, 'b-');
set(p3, 'LineWidth', 0.2); 
   
ylim([0 260])
xlim([-2.5 0])        
pagesize = [3 2]; %(x, y) size in cm

legend('boxoff'); 
legend('off');

savefig('23_02_14_loglog comparison TNF alpha SP-.fig')

set (gca, ...
    'ActivePositionProperty', 'position', ...
    'position', [0.2 0.2 0.75 0.75], ...
    'ytick', 0 : 20 : 260, ...
    'xtick', -2.5 : 0.5 : 0, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ticklength', [0.096 0.24] ./ max (pagesize), ...
    'tickdir', 'out', ...
    'linewidth', 0.5, ...
    'fontsize', 6, ...
    'fontname', 'arial', ... 
    'box', 'off');
tprint ( ...
    './23_02_14_loglog_comp_TNF alpha SP-', ...
    '-HR -jpg -eps', ...
    pagesize);



%% plot the logarithm of the data for each individual cell 

 
for no = 1:numel(spinearea) 
        
%     spinearea{no} = nonzeros(spinearea{no});
    norm_d{no} = log10(spinearea{no}); 
    start_analysis{no} = min(norm_d{no});
    end_analysis{no} = max(norm_d{no}); % the endpoint is your biggest data point 

    % creating size categories or whatever else categories you need for the
    % analysis
    ed{no} = start_analysis{no}:binsize:end_analysis{no}; % the first number is your starting point, the middle number is the steps / size of the bins
    % the last number is the end point for the bins; 
    %%% change your start and endpoints and steps in between accordingly

    [num_s{no}, eds{no}] = histcounts(norm_d{no},ed{no}); % this function gives you the number of bins with the specified parameter above
    % and the number of data points in each bin

    ed_f{no} = ed{no}(1:length(num_s{no}));
        
    xa{no} = ed_f{no}'; 
    ya{no} = num_s{no}'; 
       
    [nfit{no}, ngof{no}] = fit(xa{no}, ya{no}, 'gauss1');
    
     close 
end

ngofs = cell2mat(ngof);
writetable(struct2table(ngofs), 'gof_logarithm_IML21ipsi_indcells.xlsx');

% contralateral side

for noc = 1:numel(spinearea_c) 
        
    norm_dc{noc} = log10(spinearea_c{noc}); 
    start_analysisc{noc} = min(norm_dc{noc});
    end_analysisc{noc} = max(norm_dc{noc}); % the endpoint is your biggest data point 

    % creating size categories or whatever else categories you need for the
    % analysis
    edc{noc} = start_analysisc{noc}:binsize:end_analysisc{noc}; % the first number is your starting point, the middle number is the steps / size of the bins
    % the last number is the end point for the bins; 
    %%% change your start and endpoints and steps in between accordingly

    [num_sc{noc}, edsc{noc}] = histcounts(norm_dc{noc},edc{noc}); % this function gives you the number of bins with the specified parameter above
    % and the number of data points in each bin

    ed_fc{noc} = edc{noc}(1:length(num_sc{noc}));
        
    xac{noc} = ed_fc{noc}'; 
    yac{noc} = num_sc{noc}'; 
       
    [nfitc{noc}, ngofc{noc}] = fit(xac{noc}, yac{noc}, 'gauss1');

    close 
end

ngofsc = cell2mat(ngofc);
writetable(struct2table(ngofsc), 'gof_loglog_IML21con_indcells.xlsx');


