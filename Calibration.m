% Test calibration (S5 of project outline)
clear; clc;
addpath('./Data_prep');

% Filtered dataset with BS implied volatilities
spx_df=readtable("Data_prep/Data/spx_quotedata20220401_filtered_optionDataWithImplVol.csv");
% Data set with discount factors (BT, QT) and ATM total implied variance
discountData_df=readtable("Data_prep/Data/spx_quotedata20220401_discountData.csv");

S0 = 4545.86;

spx_df.Properties.VariableNames
discountData_df.Properties.VariableNames

% Find sigma_ask and sigma_bid for the option dataset
% σask (Ti , Kj ) is the smallest Black-Scholes implied volatility for the call and put ask prices at maturity Ti and strike Kj
% σbid(Ti,Kj) is the largest Black-Scholes implied volatility for the call and put bid prices at maturity Ti and strike Kj
col_names = {'TimeToExpiration', 'Strike', 'sigma_ask', 'sigma_bid'};
col_types = {'double', 'double','double', 'double'}; 

% Create an empty table with specified column names and data types
option_df = table('Size', [0, length(col_names)], 'VariableNames', ...
    col_names, 'VariableTypes', col_types);

% [Start with nested loops, then try to make more efficient]
for i = 1:length(discountData_df.T)
    Ti = discountData_df.T(i);
    Ks = spx_df.Strike(spx_df.TimeToExpiration == Ti);
    %for Kj = spx_df.Strike(spx_df.TimeToExpiration == Ti)
    for j = length(Ks)
        Kj = Ks(j);
        % Find sigma_ask from call and put ask vols
        callAsk_BSvol = spx_df.callAsk_BSvol(spx_df.TimeToExpiration == Ti & ...
            spx_df.Strike == Kj);
        putAsk_BSvol = spx_df.putAsk_BSvol(spx_df.TimeToExpiration == Ti & ...
            spx_df.Strike == Kj);
        sigma_ask = min(callAsk_BSvol, putAsk_BSvol);
        % Find sigma_bid from call and put bid vols
        callBid_BSvol = spx_df.callBid_BSvol(spx_df.TimeToExpiration == Ti & ...
            spx_df.Strike == Kj);
        putBid_BSvol = spx_df.putBid_BSvol(spx_df.TimeToExpiration == Ti & ...
            spx_df.Strike == Kj);
        sigma_bid = max(callBid_BSvol, putBid_BSvol);
        % Add the values to the table
        new_row = {Ti, Kj, sigma_ask,sigma_bid};
        option_df = [option_df; new_row];

    end
end


%% 
% Define objective function for optimisation problem (calibration)
% Columns of option_df: TimeToExpiration, Strike, sigma_ask, sigma_bid
function f = obj_fnc(discountData_df, option_df, eta, rho) 
    % Define constants (S3)
    gamma1 = 0.238; gamma2 = 0.253; 
    beta1 = exp(5.18); beta2 = exp(-3);
    eta = @(eps) 2.016048*exp(eps);
    
    % Define functions phi and w
    phi = @(eta, theta) eta/(theta^gamma1*(1+beta1*theta)^gamma2* ...
        (1+beta2*theta)^(1-gamma1-gamma2));
    w = @(eta, thetaT, rho, T, k) thetaT/2*(1+rho*phi(thetaT)*k + ...
        sqrt((phi(thetaT)*k + rho)^2 + (1-rho^2)));

    min_f = 0;
    
    % [Start with nested loops, then try to make more efficient]
    %for Ti=discountData_df.T
    for i = 1:length(discountData_df.T)
        Ti = discountData_df.T(i);
        Ks = option_df.Strike(option_df.TimeToExpiration == Ti);
        %for Kj = spx_df.Strike(spx_df.TimeToExpiration == Ti)
        for j = length(Ks)
            Kj = Ks(j);
            kj = option_df.logStrike(option_df.TimeToExpiration == Ti & ...
                option_df.Strike == Kj);
            w(eta, thetaT, rho, T, k)
        end
    end

end

% Columns of filteredOptionData:
% TimeToExpiration, Strike, CallMktPrice, PutMktPrice, TotalImplVar, 
% logStrike, CallBidPrice, CallAskPrice, PutBidPrice, PutAskPrice
%function res = CalibrationObjFnc(filteredOptionData, ) 