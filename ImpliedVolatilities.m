% Implied volatilities
% Test if plots in report outline can be replicated
clear; clc;
addpath('./Data_prep');

spx_df=readtable("Data_prep/Data/spx_quotedata20220401_filtered_optionDataWithImplVol.csv");
discountData_df=readtable("Data_prep/Data/spx_quotedata20220401_discountData.csv");

S0 = 4545.86;

% Find time to expiration corresponding to 90 days
T_vals = unique(spx_df.TimeToExpiration); % Maturities from option data
T_below = max(T_vals(T_vals <= 90/365));
T_above = min(T_vals(T_vals >= 90/365));

filter90 = (spx_df.TimeToExpiration >= T_below) & (spx_df.TimeToExpiration < T_above);
figure(1)
plot(spx_df.logStrike(filter90), spx_df.callBid_BSvol(filter90), ".");
hold on
plot(spx_df.logStrike(filter90), spx_df.callAsk_BSvol(filter90), ".");
hold on
plot(spx_df.logStrike(filter90), spx_df.putBid_BSvol(filter90), ".");
hold on
plot(spx_df.logStrike(filter90), spx_df.putAsk_BSvol(filter90), ".");
hold off
xlabel("Log strike")
ylabel("BS implied vol")
legend(["Call bid", "Call ask", "Put bid", "Put ask"])
title("SPX: maturity 90 days (2022-4-1)")