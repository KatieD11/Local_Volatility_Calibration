% Implied volatilities
% Compute implied volatilies using SSVI formulation 
% with parameters from calibration
clear; clc;
addpath('./Data_prep');

spx_df=readtable("Data_prep/Data/spx_quotedata20220401_filtered_optionDataWithImplVol.csv");
discountData_df=readtable("Data_prep/Data/spx_quotedata20220401_discountData.csv");
calibration_params=readtable("Calibration_results/spx_20220401_calibration_params.csv");
%calibration_params=readtable("Calibration_results/spx_20220401_calibration_params_without_weights.csv");
bid_ask_spread=readtable("Calibration_results/spx_20220401_calibration_bid_ask_spread.csv");

S0 = 4545.86; % spx_20220401
%% Compute implied vols for all maturities in dataset
SSVI_vols_cellArray = {}; % Store arrays of different lengths
ks_cellArray = {};

T_maturities = discountData_df.T; 
for i = 1:length(T_maturities)
    T_i = T_maturities(i);
    ks = spx_df.logStrike(spx_df.TimeToExpiration == T_i);
    thetaT = discountData_df.TotImplVar(discountData_df.T == T_i);
    SSVI_vols = SSVIimpliedVolatility(thetaT, T_i, ks, ...
            calibration_params.rho, calibration_params.eps);
    % Store arrays
    ks_cellArray{end+1} = ks;
    SSVI_vols_cellArray{end+1} = SSVI_vols;
end

%% Plot results for a particular maturity T
% Choose a time to maturity in days
Tn_days = 90;
%Tn_days = 105; % [17, 19, 21, 24, 26, 28, 31, 35, 42, 49 ..., ]

% Get the time to expiration in days (approx)
spx_df.T_days = round(spx_df.TimeToExpiration*365);
discountData_df.T_days = round(discountData_df.T*365);
bid_ask_spread.T_days = round(bid_ask_spread.TimeToExpiration*365);

% Strike values for the maturity
ks = ks_cellArray{discountData_df.T_days == Tn_days};

filter = spx_df.T_days == Tn_days;
figure(1)
plot(ks, spx_df.callBid_BSvol(filter), ".b");
hold on
plot(ks, spx_df.callAsk_BSvol(filter), ".r");
hold on
plot(ks, spx_df.putBid_BSvol(filter), ".",'Color',"#0072BD");
hold on
plot(ks, spx_df.putAsk_BSvol(filter), ".",'Color',"#A2142F");
hold on
%plot(ks, SSVI_vols, "-g", "LineWidth",1);
plot(ks, SSVI_vols_cellArray{discountData_df.T_days == Tn_days}, "-g", "LineWidth",1);
hold on
plot(bid_ask_spread.logStrike(bid_ask_spread.T_days == Tn_days), ...
    bid_ask_spread.sigma_target(bid_ask_spread.T_days== Tn_days), ...
    "-c", "LineWidth",1);
hold off
xlabel("Log strike")
ylabel("BS implied vol")
legend(["Call bid", "Call ask", "Put bid", "Put ask", "SSVI", "Target"])
title("SPX: maturity " +Tn_days+ " days (2022-4-1)")
%% Implied variance plot
figure(2)
for i = 1:length(T_maturities)
    T_i = T_maturities(i);
    ks = ks_cellArray{discountData_df.T == T_i};
    implied_var = (SSVI_vols_cellArray{discountData_df.T == T_i}).^2*T_i;
    plot(ks, implied_var, "-", "LineWidth",1);
    hold on
end
hold off
legend(string(T_maturities))
title("Total implied variance plot, SPX (2022-4-1)")
xlabel("Log strike")
ylabel("Total implied variance")
%% Check subset of the smaller maturities
T_small = T_maturities(T_maturities < 0.1);
figure(3)
for i = 1:length(T_small)
    T_i = T_small(i);
    ks = ks_cellArray{discountData_df.T == T_i};
    implied_var = (SSVI_vols_cellArray{discountData_df.T == T_i}).^2*T_i;
    plot(ks, implied_var, "-", "LineWidth",1);
    hold on
end
hold off
legend(string(T_small))
title("Subset: Total implied variance plot, SPX (2022-4-1)")
xlabel("Log strike")
ylabel("Total implied variance")
