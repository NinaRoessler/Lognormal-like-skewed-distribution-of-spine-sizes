% lognormal fit to data

function [log_plot, gof_stats, coeff, logfit, xa, yaf] = logfitting(spinea, binsize, mulow, silow, muup, siup) 
      
    % fit all spines from all cells into lognormal fit 
%     start_analysis = min(spinea);
    start_analysis = 0.001; 

    end_analysis = max(spinea); % the endpoint is your biggest data point 

    % creating size categories or whatever else categories you need for the
    % analysis
    ed = start_analysis:binsize:end_analysis; % the first number is your starting point, the middle number is the steps / size of the bins
    % the last number is the end point for the bins; 
    %%% change your start and endpoints and steps in between accordingly

    [num_s, eds] = histcounts(spinea,ed); % this function gives you the number of bins with the specified parameter above
    % and the number of data points in each bin

    ed_f = ed(1:length(num_s));

    % ex_all = nonzeros(spinea); % excludes values that are 0 
    % ex_all = ex_all(~isnan(ex_all)); % excludes all NaN variables in script
    mu_all = mean(log(spinea)); % calculates the mean of your data
    s_all = std(log(spinea)); % calculates the std of your data


    xa = ed_f'; 
    ya = num_s';

    % normalizing the y-axis / the number of cells in each bin otherwise the
    % fit and the raw data could not be plotted in the same figure
%     tot_n = numel(spinea); 
    Aa = trapz(xa,ya); 
    yaf = ya / Aa; % normalizing by dividing the data points with the area calculated above
      
    % writing the function for the fit you want to do - in this case,
    % probability density function for the lognormal distribution
    % define your variables and the parameters / coefficients you want to use
    logfittype_all = fittype('(1./(xa.*s_all.*sqrt(2.*pi))) .* exp(-((log(xa)-mu_all).^2) ./ (2.*s_all.^2))',...
        'dependent',{'ya'},'independent',{'xa'},...
        'coefficients',{'mu_all','s_all'});


    % fitting the fittype over the data and as last step plotting both
    % put in the data, your fit-function (above), borders for the analysis and
    % fitting and start and endpoints for the parameters

    %%% you can change Lower and Upper and Start; these are the values the
    %%% fit-algorithm takes into account to make the fitting faster; the two
    %%% numbers are for the coefficients [mu s] 
    
%     [logfit, gof] = fit(xa,yaf, logfittype_all, 'Start', [mu_all, s_all]);
      
    [logfit, gof] = fit(xa,yaf, logfittype_all, 'Lower', [mulow silow], 'Upper', [muup siup], ...
            'Start', [mu_all, s_all]);
  
    % getting the goodness of fit statistics and your mu and sigma for the fitted curve 
    % you can save them if needed 
    gof_stats = [gof];
    coeff = coeffvalues(logfit);
    
    figure
    log_plot = plot(logfit, xa, yaf);
    
       
         %%% change title for each category you have %%%
         % Create title of plot
%          title(sprintf(t));
               
         % Create name of y-axis 
%          ylabel('Normalized Unit');

         % Create name of x-axis (change accordingly)
%          xlabel('Spine sizes (µm3)');
    
         %%% change the legend labels, location of legend
         legend('Data', 'Fitted curve', 'Location', 'northeast');
         

end