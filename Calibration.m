% Test calibration (S5 of project outline)
clear; clc;
spx_df=readtable("Data_prep/Data/spx_quotedata20220401_filtered_optionData.csv");
%discountData_df=readtable("Data_prep/Data/spx_quotedata20220401_discountData.csv");

S0 = 4545.86;

% Black-Scholes implied volatilities (using interpolated prices)
% Note: blsimpv(Price,Strike,Rate,Time,Value, [Limit], [Yield], [Class])
spx_df.callBid_BSvol = zeros(length(spx_df.TimeToExpiration),1);
spx_df.callAsk_BSvol = zeros(length(spx_df.TimeToExpiration),1);
for i = 1: length(spx_df.TimeToExpiration)
    spx_df.callBid_BSvol(i) = ...
        blsimpv(S0,spx_df.Strike(i),spx_df.rT(i),spx_df.TimeToExpiration(i), ... 
        spx_df.CallBidPrice(i),'Class', {'Call'});
    spx_df.callAsk_BSvol(i) = ...
        blsimpv(S0,spx_df.Strike(i),spx_df.rT(i),spx_df.TimeToExpiration(i), ... 
        spx_df.CallAskPrice(i),'Class', {'Call'});
end

% Find time to expiration corresponding to 90 days
T_vals = unique(spx_df.TimeToExpiration); % Maturities from option data
T_below = max(T_vals(T_vals <= 90/365));
T_above = min(T_vals(T_vals >= 90/365));

filter90 = (spx_df.TimeToExpiration >= T_below) & (spx_df.TimeToExpiration < T_above);
figure(1)
plot(spx_df.logStrike(filter90), spx_df.callBid_BSvol(filter90), ".");
hold on
plot(spx_df.logStrike(filter90), spx_df.callAsk_BSvol(filter90), ".");
hold off
xlabel("Log strike")
ylabel("BS implied vol")
legend(["Call bid", "Call ask"])

%% 

% Define objective function for optimisation problem
% Columns of filteredOptionData:
% TimeToExpiration, Strike, CallMktPrice, PutMktPrice, TotalImplVar, 
% logStrike, CallBidPrice, CallAskPrice, PutBidPrice, PutAskPrice
%function res = CalibrationObjFnc(filteredOptionData, ) 