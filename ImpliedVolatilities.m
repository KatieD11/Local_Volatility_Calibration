% Implied volatilities
% Test if plots in report outline can be replicated
clear; clc;
addpath('./Data_prep');

spx_df=readtable("Data_prep/Data/spx_quotedata20220401_filtered_optionDataWithImplVol.csv");
discountData_df=readtable("Data_prep/Data/spx_quotedata20220401_discountData.csv");
calibration_params=readtable("spx_20220401_calibration_params.csv");

S0 = 4545.86;

% Choose a time to maturity in days
Tn_days = 90;

% Get the time to expiration in days (approx)
spx_df.T_days = round(spx_df.TimeToExpiration*365);
discountData_df.T_days = round(discountData_df.T*365);

% Compute implied vols using SVI framework for Tn_days
ks = spx_df.logStrike(spx_df.T_days == Tn_days);
SSVI_vols = zeros(length(ks),1);
thetaT = discountData_df.TotImplVar(discountData_df.T_days == Tn_days);
T_i = discountData_df.T(discountData_df.T_days == Tn_days);
for j=1:length(ks)
    SSVI_vols(j) = SSVIimpliedVolatility(thetaT, T_i, ks(j), ...
        calibration_params.rho, calibration_params.eps);
end

filter = spx_df.T_days == Tn_days;
figure(1)
plot(ks, spx_df.callBid_BSvol(filter), ".");
hold on
plot(ks, spx_df.callAsk_BSvol(filter), ".");
hold on
plot(ks, spx_df.putBid_BSvol(filter), ".");
hold on
plot(ks, spx_df.putAsk_BSvol(filter), ".");
hold on
plot(ks, SSVI_vols, "-", "LineWidth",2);
hold off
xlabel("Log strike")
ylabel("BS implied vol")
legend(["Call bid", "Call ask", "Put bid", "Put ask", "SSVI"])
title("SPX: maturity " +Tn_days+ " days (2022-4-1)")
%% 
% Try create a surface plot
