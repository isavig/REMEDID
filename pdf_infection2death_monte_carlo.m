%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%
% This function simulates N Probability Density Functions (PDF) of the
% random variable representing the period from infection to death (I2D). 
% The PDF are obtainted convolving N PDF representing:
% - Incubation Period (IP)
% - Illness Onset to Death (IOD) period
% which are obtained from N random combinations (within error estimates) of
% the parameters defining the PDF of IP and IOD.
% The N I2D PDF will be used to estimate errors applying the REMEDID 
% algorithm in a Monte Carlo analysis.
%
% REMEDID: Retrospective Methodology to Estimate Daily Infections from 
%          Deaths  
%
% Written by 
%                       David Garcia-Garcia
%                       University of Alicante, Spain
%                       d.garcia@ua.es
%
%                                                            September 2023
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
%   - N: number of I2D PDF to be estimated,
%
%   - N_days: number of days where the PDF will be estimated
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
% REMARK: Parameters 1 and 2 are different for each distribution. See 
% Matlab documentation for further details. In Matlab R2021b:
% - Lognormal: 
%       - Parameter_1: Mean 
%       - Parameter_2: Median 
%
% OUTPUTS:
%   - x_days: Days where the PDF is estimated. 
%             Row vector: 1 x N_days
%
%   - y_pdf_N_days: Daily values of the PDF.
%              Matrix: N x N_days (each row represents a PDF)




function [x_days, y_pdf_N_days] = ...
    pdf_infection2death_monte_carlo(N, N_days,...
                              IP_distribution,  ...
                              IP_parameter_1_lower_error,  IP_parameter_1_upper_error,...
                              IP_parameter_2_lower_error,  IP_parameter_2_upper_error,...
                              IOD_distribution, ...
                              IOD_parameter_1_lower_error,  IOD_parameter_1_upper_error,...
                              IOD_parameter_2_lower_error,  IOD_parameter_2_upper_error)



% Log-Normal
% Y ~ log-N(mu, sigma^2)
% ln(Y) ~ N(mu, sigma^2)
%
% median(Y) = exp(mu)
% mean(Y) = exp( mu + (sigma^2)/2 )
%
% Then, mu and sigma can be estimated from the mean and the median as:
% mu = log( median(Y) )
% sigma^2 = 2 * ln( mean(Y) / median(Y) ) 



h_paso = 0.1; %Longitud de paso
x_ip = 0:h_paso:N_days;
x_i2d = 0:h_paso:N_days;

num_datos_convo = 2*length(x_ip)-1;

y_pdf_N  = zeros(N, num_datos_convo) *NaN;



if strcmp('Lognormal', IP_distribution) && strcmp('Lognormal', IOD_distribution)
    
    for i=1:N
        
        N - i
    
        %% Incubation Period
        control_ip = 0;
        while control_ip==0
            mean_Y_ip   = random('uniform', IP_parameter_1_lower_error, IP_parameter_1_upper_error);
            median_Y_ip = random('uniform', IP_parameter_2_lower_error, IP_parameter_2_upper_error);
            
            %mu y sigma
            mu_ip     =  log( median_Y_ip );
            sigma2_ip =  2 *  log( mean_Y_ip / median_Y_ip ) ;
            
            %sigma must be positive
            if sigma2_ip>0
                %PDF
                pd_ip = makedist('Lognormal','mu', mu_ip  ,'sigma', sqrt(sigma2_ip)  );
                y_ip = pdf(pd_ip, x_ip);
                control_ip = 1;
            end
    
        end %while
    
    
    
        %% Illnes onset to death
        control_i2d = 0;
        while control_i2d==0
            mean_Y_i2d   = random('uniform', IOD_parameter_1_lower_error, IOD_parameter_1_upper_error);
            median_Y_i2d = random('uniform', IOD_parameter_2_lower_error, IOD_parameter_2_upper_error);
            
            mu_i2d     =  log( median_Y_i2d );
            sigma2_i2d =  2 *  log( mean_Y_i2d / median_Y_i2d ) ;
            
            %sigma must be positive
            if sigma2_i2d>0
                %PDF
                pd_i2d = makedist('Lognormal','mu', mu_i2d  ,'sigma', sqrt(sigma2_i2d)  );
                y_i2d = pdf(pd_i2d, x_i2d);
                control_i2d = 1;
            end
    
        end %while
    
        %% ----------------------------------------
        %% Convolution
        %------------------------------------------
        % We generate a density function by convolving the two density 
        % functions.
        % For each x of ip we add the curve of i2d multiplied by the 
        % probability of x. Then we add all the resulting curves together. 
        
        prob_all = zeros(length(x_ip), num_datos_convo )*NaN;
        
        for x=1:length(x_i2d)
            
            prob_all(x, x:x+length(x_ip)-1 ) = y_ip(x) * y_i2d * h_paso;
            
        end
        
        y_pdf_N(i,:) = nansum(prob_all, 1);
    
    
    
    end %for

end %if



%% PDF, daily values

% Daily version of the convolution: Infection to Death (ID)
x_days = 0: ((size(y_pdf_N,2)-1) *h_paso) -1;

y_pdf_N_days = zeros(N, N_days) *NaN;

k=0;
for i=1 : 1/h_paso : size(y_pdf_N,2)-1
    k=k+1;
    y_pdf_N_days(:,k) = mean(y_pdf_N(:,i:i+9), 2);
end









