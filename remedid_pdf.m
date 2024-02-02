%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%
% This function implements de REMEDID algorithm
% REMEDID: Retrospective Methodology to Estimate Daily Infections from 
%          Deaths  
%
% In this version the probability density function (PDF) of the 
% "Infection to death period" is an input.
%
% Written by 
%                       David Garcia-Garcia
%                       University of Alicante, Spain
%                       d.garcia@ua.es
%
%                                                                July 2023
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
% If you use this function, please cite the following publication:
% Garcia-Garcia, D., I. Vigo, E. S. Fonfria, Z. Herrador, M. Navarro, and
% C. Bordehore. Retrospective Methodology to Estimate Daily Infections from
% Deaths (REMEDID) in COVID-19: the Spain case study. Scientific Reports, 
% 11:11274, 2021. https://doi.org/10.1038/s41598-021-90051-7
%-------------------------------------------------------------------------
%
% INPUTS:
%   - deaths: time series of daily deaths. It must be a row or a column.
%
%   - pdf: Probability Density Function (PDF) of the "Infection to death
%             period". Each value represents a day. It must be a row or a 
%             column.
%
%   - N_truncation: Days of the PDF used
% 
%   - min_percentage: Infections close to the end of the time series are 
%                     inferred only with a percentage of their associated 
%                     deaths. The infection time series is truncated in 
%                     such way that infections are estimated at least with 
%                     the minimum of percentage defined by the parameter. 
%                     Usually, min_percentage = 95 (%)
%   - CFR: Case Fatality Ratio in percentage. It can be a number (then it
%          is constant along the time series) or a vector (then there is a
%          CFR for each day) with the same dimensions than deaths.
%   - plot_option: If plot_option=1, the probability distribution functions 
%                  are plotted
%
% OUTPUTS:
%   - infections: time series of daily infections
%   - Nmin: Last Nmin days of the "infections" time series must be 
%           eliminated to accomplish with min_percentage  
%   - Nextra: "infections" time series is longer than "deaths" time series
%             in "Nextra" days.
%             Nextra is the number of days from first infection to first 
%             day in "deaths" time series. 
%             It may happend that first death in "deaths" is not in death(1)
%             and first infection is later than date of death(1). In that
%             case, Nextra = 0.




function [infections, Nmin, Nextra, infections_aux] = remedid_pdf(deaths, ...
                              pdf, N_truncation, min_percentage, CFR,...
                              plot_option)


                          
N=N_truncation;     %Days to evaluate the PDF. 60 is enough for COVID-19. You may like change it
pdf_ID = pdf(1:N);


%% --------------
% Plot (optional)
% ---------------

% Only plotted if plot_option==1
if plot_option==1
    
    figure('Renderer', 'painters', 'Position', [10 10 700 350])
    
    hold on, grid on
    h_convolution = plot(1:N, pdf_ID,'linewidth',2);
    
    xlim([0,N+1])
    xlabel('Days ','FontSize',15)
    
    aa = get(gca,'XTickLabel');
    set(gca,'XTickLabel',aa,'fontsize',13)
    
    title('Probability Density Function ', 'fontweight','bold','FontSize',25);
    
end






%% ---------------------------------------------------------
% Cumulative Distribution Function (CDF) of pdf_convolution
%-----------------------------------------------------------

for i=1:length(pdf_ID)
    cdf_ID(i) = sum(pdf_ID(1:i));
end

ind = find(cdf_ID < min_percentage/100 );

%Minimal number of days needed accordingly to min_percentage
Nmin = ind(end);




%% ----------
% Infections 
%------------

% This number must be larger than the number of days between the first day 
% of inferred infections and the first day of deaths time series. It will
% be adjust later.
Nextra_aux = 200;

%Check dimension of deaths:
[rows,colu] = size(deaths);
if     rows==1  
elseif colu==1, deaths = deaths';
else, error('Check dimensions of deaths time series')
end

%Check dimension of CFR
[rows_cfr,colu_cfr] = size(CFR);
if     rows_cfr==1 && colu_cfr==1 
    CFR2 = CFR * ones(1,length(deaths) );
elseif rows_cfr==1 && colu_cfr==length(deaths)
    CFR2 = CFR;
elseif rows_cfr==length(deaths) && colu_cfr==1
    CFR2 = CFR';
else, error('Check dimensions of CFR')
end


% Nextra zeros are added before the beginning of deaths and CFR2 time series
deaths_extended = [zeros(1, Nextra_aux), deaths];
CFR_extended    = [ones(1, Nextra_aux)*CFR2(1), CFR2];

         
Nextended = length(deaths_extended);
infections_CFR100 = zeros(1, Nextended);


for i=1:Nextended
    
    aux = deaths_extended(i:end);
    
    if length(aux)>N %+1
        k = length(aux) - (N);  %+1);
        aux_pd = [pdf_ID, zeros(1,k)];
    else
        aux_pd = pdf_ID(1:length(aux));
    end
    
    % Equation 1:
    infections_CFR100(i) = sum(  aux .* aux_pd  );
    
end
    


% Equation 2:
infections_aux = round( infections_CFR100 ./CFR_extended *100 );


% Time series is cut to avoid zeros at the begining:
ind_positive = find(infections_aux > 0);

if ind_positive(1) <= Nextra_aux
    infections = infections_aux(ind_positive(1):end);
else
    infections = infections_aux(Nextra_aux+1:end);
end

Nextra = length(infections) - length(deaths);











