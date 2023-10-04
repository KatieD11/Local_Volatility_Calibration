% Calibration (basic version)
% First run SPX_analysis file and SPX_bid_ask_spread (in Data_prep folder)
% to get filtered option data (with BS implied vols),
% discount data, and bid-ask spreads
clear; clc;

% Select data set
dataset = "spx_20220401";

% Filtered dataset with BS implied volatilities
spx_df=readtable("../Data_prep/Data/"+dataset+"_filtered_optionDataWithImplVol.csv");
% Data set with discount factors (BT, QT) and ATM total implied variance
discountData_df=readtable("../Data_prep/Data/"+dataset+"_discountData.csv");
% Implied vol bid-ask spread data
option_df =readtable("../Data_prep/Data/"+dataset+"_option_bid_ask_spread.csv");

% Intial value for spx_20220401
S0 = 4545.86;
%% 
f = @(params) obj_fnc(discountData_df, option_df, params(1), params(2));
%eps = 0.8467; rho = -0.6887;
opt_params = fminsearch(f, [0.5, 0.5]);
eps_opt = opt_params(1)
rho_opt = opt_params(2)
cost = obj_fnc(discountData_df, option_df, eps_opt, rho_opt)
% Save params in csv file
calibration_params = table;
calibration_params.eps = eps_opt;
calibration_params.rho = rho_opt;
% Store constants used for calibration
gamma1 = 0.238; gamma2 = 0.253; 
beta1 = exp(5.18); beta2 = exp(-3);
calibration_params.gamma1 = gamma1;
calibration_params.gamma2 = gamma2;
calibration_params.beta1 = beta1;
calibration_params.beta2 = beta2;
% Save cost
calibration_params.cost = cost;
writetable(calibration_params, "Calibration_results/"+dataset+"_calibration_params_without_weights.csv")
%% 
% Define objective function for optimisation problem (calibration)
% Columns of option_df: TimeToExpiration, Strike, logStrike, sigma_ask, sigma_bid
% Columns of discountData_df: T, BT, QT, TotImplVar, rT
function summation = obj_fnc(discountData_df, option_df, eps, rho) 
    % Define constants (S3)
    gamma1 = 0.238; gamma2 = 0.253; 
    beta1 = exp(5.18); beta2 = exp(-3);
    eta =  2.016048*exp(eps);
    
    % Define functions phi and w
    phi = @(theta) eta/(theta^gamma1*(1+beta1*theta)^gamma2* ...
        (1+beta2*theta)^(1-gamma1-gamma2));
    w = @(thetaT, T, k) thetaT/2*(1+rho*phi(thetaT)*k + ...
        sqrt((phi(thetaT)*k + rho)^2 + (1-rho^2)));

    summation = 0;
    
    % [Start with nested loops, then try to make more efficient]
    %for Ti=discountData_df.T
    for i = 1:length(discountData_df.T)
        Ti = discountData_df.T(i);
        thetaTi = discountData_df.TotImplVar(i);
        Ks = option_df.Strike(option_df.TimeToExpiration == Ti);
        %for Kj = spx_df.Strike(spx_df.TimeToExpiration == Ti)
        for j = 1:length(Ks)
            Kj = Ks(j);
            % Filter by Ti and Kj
            filter = option_df.TimeToExpiration == Ti & ...
                option_df.Strike == Kj;
            kj = option_df.logStrike(filter);
            sigma_ask = option_df.sigma_ask(filter);
            sigma_bid = option_df.sigma_bid(filter);

            summation = summation + (max(0, w(thetaTi, Ti, kj) ...
                - sigma_ask^2*Ti)^2 + min(0, w(thetaTi, Ti, kj) ...
                - sigma_bid^2*Ti)^2)/w(thetaTi, Ti, kj)^2;
        end
    end

end