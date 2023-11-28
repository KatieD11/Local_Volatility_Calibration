% Calibration (using power-law phi fnc)
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
%% Solve optimisation problem
f = @(params) obj_fnc(discountData_df, option_df, params(1), ...
    params(2),params(3));
%rho = -0.6887;
%opt_params = fminsearch(f, [0.5, 0.5]);
lb = [-Inf, -1, 0];
ub = [Inf, 1, 0.5];
opt_params = fmincon(f,[0.5, 0.5, 0.2],[],[],[],[],lb,ub);
eta_opt = opt_params(1)
rho_opt = opt_params(2)
gamma_opt = opt_params(3)
cost = cost_fnc(discountData_df, option_df, eta_opt, rho_opt, gamma_opt)
cost_w = obj_fnc(discountData_df, option_df, eta_opt, rho_opt, gamma_opt)
%% Save params in csv file
calibration_params = table;
calibration_params.eta = eta_opt;
calibration_params.rho = rho_opt;
calibration_params.gamma = gamma_opt;
% Save cost (without and with weights)
calibration_params.cost = cost;
calibration_params.cost_w = cost_w;
writetable(calibration_params, "Calibration_results/"+dataset+"_calibration_params_powerLaw.csv")
%% 
% Define objective function for optimisation problem (calibration)
% Columns of option_df: TimeToExpiration, Strike, logStrike, sigma_ask, sigma_bid
% Columns of discountData_df: T, BT, QT, TotImplVar, rT
function summation = obj_fnc(discountData_df, option_df, eta, rho, gamma) 
    % Check that static arbitrage free condition is met
    if (eta*(1+abs(rho)) > 2)
        summation = 10^6; % Set cost high
        return
    end

    % Define functions phi (power-law) and w
    phi = @(theta) eta/(theta^gamma*(1+theta)^(1-gamma));
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
            weight = option_df.open_interest_weight(filter);

            summation = summation + weight*((max(0, w(thetaTi, Ti, kj) ...
                - sigma_ask^2*Ti)^2 + min(0, w(thetaTi, Ti, kj) ...
                - sigma_bid^2*Ti)^2)/w(thetaTi, Ti, kj)^2);
        end
    end

end

% Compute the cost without the open weights 
% (to compare across calibration techniques)
function summation = cost_fnc(discountData_df, option_df, eta, rho, gamma) 
    % Check that static arbitrage free condition is met
    if (eta*(1+abs(rho)) > 2)
        disp("Static arbitrage exists")
        return
    end
    
    % Define functions phi (power-law) and w
    phi = @(theta) eta/(theta^gamma*(1+theta)^(1-gamma));
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