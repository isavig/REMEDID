%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%
% This function is a modification of the original remedid.m function,
% which implemented de REMEDID algorithm. 
% REMEDID: Retrospective Methodology to Estimate Daily Infections from 
%          Deaths  
%
% This new version applies a Monte Carlo simulation to estimate errors in
% the estimated infections.
%
% Written by 
%                       David Garcia-Garcia
%                       University of Alicante, Spain
%                       d.garcia@ua.es
%
%                                                           September 2023
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% 
% This algorithm has been designed to estimate infections during the 
% COVID-19 pandemic, although it can be used for any disease causing death.
%
% This function is freely distributed without any warranty. It has been
% tested for Matlab R2021b
%
%-------------------------------------------------------------------------
% If you use this function, please cite the following publications:
% - Garcia-Garcia, D., I. Vigo, E. S. Fonfria, Z. Herrador, M. Navarro, and
% C. Bordehore. Retrospective Methodology to Estimate Daily Infections from
% Deaths (REMEDID) in COVID-19: the Spain case study. Scientific Reports, 
% 11:11274, 2021. https://doi.org/10.1038/s41598-021-90051-7
%
% - Marquez, J., D. Garcia-Garcia, M. I. Vigo, and C. Bordehore. Estimacion 
% retrospectiva de los casos iniciales de COVID-19 en Santiago Region 
% Metropolitana en Chile. Gaceta Sanitaria. 2024.
%-------------------------------------------------------------------------
%
% INPUTS:
%   - deaths: time series of daily deaths. It must be a row or a column.
%   - dates_deaths: dates associated to time series of daily deaths. It 
%         must be a row or a column.
%
%   - Incubation period (IP) distribution:
%       - IP_distribution: 'Lognormal'
%       - IP_parameter_1_lower_error
%       - IP_parameter_1_upper_error
%       - IP_parameter_2_lower_error
%       - IP_parameter_2_upper_error
%       
%   - Illness onset to death (IOD)  distribution:
%       - IOD_distribution: 'Lognormal'
%       - IOD_parameter_1_lower_error
%       - IOD_parameter_1_upper_error
%       - IOD_parameter_2_lower_error
%       - IOD_parameter_2_upper_error
%       
%   - CI_percentage: percentage of the confidence interval error of the
%                    estimate. Value from 0 to 100
%
%   - N_mc: number of Monte Carlo simulations. Usually it is a large value
%        (several thousands)
%
%   - N_days: number of days where the Probability Density Functions (PDF) 
%             of the infection to death (I2D) will be estimated.
%             For COVID it is usually 60.
%
%   - min_percentage: Infections close to the end of the time series are 
%                     inferred only with a percentage of their associated 
%                     deaths. The infection time series is truncated in 
%                     such way that infections are estimated at least with 
%                     the minimum of percentage defined by the parameter. 
%                     Usually, min_percentage = 95 (%)
%
%   - CFR: Case Fatality Ratio in percentage. It can be a number (then it
%          is constant along the time series) or a vector (then there is a
%          CFR for each day) with the same dimensions than deaths.
%
%   - plot_option: If plot_option=1, the probability distribution functions 
%                  are plotted
%
% Parameters 1 and 2 are different for each distribution. See Matlab 
% documentation for further details. In Matlab R2021b:
% - Lognormal: 
%       - Parameter_1: Mean 
%       - Parameter_2: Median 
%       
%
% OUTPUTS:
%   - infections_sim: all simulated time series 
%
%   - infections_mean: averaged time series of all simulated daily infections
%
%   - dates_infections: dates associated to time series of daily infections.
% 
%   - infections_CI_lower_error
%
%   - infections_CI_upper_error



function [x_days, y_pdf_N_days, ...
          infections_sim, infections_mean, dates_infections,...
          infections_CI_lower_error, infections_CI_upper_error]...
        = remedid_monte_carlo_error(deaths, dates_deaths,...
                              IP_distribution,  ...
                              IP_parameter_1_lower_error,  IP_parameter_1_upper_error,...
                              IP_parameter_2_lower_error,  IP_parameter_2_upper_error,...
                              IOD_distribution, ...
                              IOD_parameter_1_lower_error,  IOD_parameter_1_upper_error,...
                              IOD_parameter_2_lower_error,  IOD_parameter_2_upper_error,...
                              CI_percentage, N_mc, N_days,...
                              min_percentage, CFR,...
                              plot_option)



%Check the CI_percentage value:
if CI_percentage < 0
    error('CI_percentage must be larger than 0')
elseif CI_percentage > 100
    error('CI_percentage must be lower than 100')
end



% Generation of the PDF of infection to death (I2D):
[x_days, y_pdf_N_days] = ...
    pdf_infection2death_monte_carlo(N_mc, N_days,...
                              IP_distribution,  ...
                              IP_parameter_1_lower_error,  IP_parameter_1_upper_error,...
                              IP_parameter_2_lower_error,  IP_parameter_2_upper_error,...
                              IOD_distribution, ...
                              IOD_parameter_1_lower_error,  IOD_parameter_1_upper_error,...
                              IOD_parameter_2_lower_error,  IOD_parameter_2_upper_error);




%% --------------------------------------
%% REMEDID errores PDF. Monte Carlo
%----------------------------------------

%[N_mc, N_days] = size(y_convo_all_daily);

N_truncation = 60; % Days to evaluate the PDF. 60 is enough for COVID-19. 
                   % You may like change it

for i=1:N_mc

    N_mc - i

    pdf = y_pdf_N_days(i,:);
    
    [infect_remedid, Nmin, Nextra] = remedid_pdf(deaths, ...
                                  pdf, N_truncation, ...
                                  min_percentage, CFR,...
                                  0);
    
    infect_remedid_monte_carlo_aux{i} = infect_remedid;
    Nextra_all(i) = Nextra;
    Nmin_all(i) = Nmin;

end


Nextra_max = max(Nextra_all);
Nmin_max   = max(Nmin_all);
Nmin = 33; %Value from the PDF estimated with the main values of the parameters.

num_datos = length(deaths);
infections_sim = zeros(N_mc, num_datos+Nextra_max-Nmin) *NaN; 
dates_remedid_monte_carlo = datetime(zeros(N_mc, num_datos+Nextra_max-Nmin), zeros(N_mc, num_datos+Nextra_max-Nmin), zeros(N_mc, num_datos+Nextra_max-Nmin), 'format','yyyy-MM-dd');



% We extend the dates Nextra days before the start of the series of deaths.
%
% We add zeros to the beginning of the series so that they all start on
% the same day: Nextra_max days before dates_deaths(1):
for i=1:N_mc

    Nextra = Nextra_all(i);

    if Nextra==Nextra_max
        dates_aux = dates_deaths(1)-Nextra : dates_deaths(end);
        infect_remedid_aux = infect_remedid_monte_carlo_aux{i};
    elseif Nextra<Nextra_max
        k=Nextra_max - Nextra;
        dates_aux = dates_deaths(1)-(Nextra+k) : dates_deaths(end);
        infect_remedid_aux = [zeros(1,k), infect_remedid_monte_carlo_aux{i}];

    end
        
     infections_sim(i,:) = infect_remedid_aux(1:end-Nmin);
     dates_remedid_monte_carlo(i,:)  = dates_aux(1:end-Nmin);

clear k
clear infect_remedid_aux
clear dates_aux

end



%% Mean and confidence interval

dates_infections = dates_remedid_monte_carlo(1,:);

b = length(dates_infections);

infections_mean           = zeros(1, b) *NaN;
infections_CI_lower_error = zeros(1, b) *NaN;
infections_CI_upper_error = zeros(1, b) *NaN;

quartile_lower = (100 - CI_percentage) /2 /100;
quartile_upper = 1 - quartile_lower;

for i=1:b
    
    aux = infections_sim(:,i);
    infections_mean(i)           = round(mean(aux));
    infections_CI_lower_error(i) = quantile(aux, quartile_lower);
    infections_CI_upper_error(i) = quantile(aux, quartile_upper);
    
end



%% Plot
if plot_option==1

    figure('Renderer', 'painters', 'Position', [10 10 700 350])
    hold on, grid on

    p=fill([dates_infections, fliplr(dates_infections)], [infections_CI_upper_error, fliplr(infections_CI_lower_error)], 'r', 'FaceAlpha', 0.2);
    p.EdgeColor = rgb('lightsalmon');

    h_infect = plot(dates_infections, infections_mean, '-', 'color', 'r', 'linewidth',0.5);

end





