% Test calibration (S5 of project outline)
clear; clc;
addpath('./Data_prep');

spx_df=readtable("Data_prep/Data/spx_quotedata20220401_filtered_optionData.csv");
discountData_df=readtable("Data_prep/Data/spx_quotedata20220401_discountData.csv");

S0 = 4545.86;

% Find BS implied vols for the bid/ask calls and puts
spx_df.callBid_BSvol = zeros(length(spx_df.TimeToExpiration),1);
spx_df.callAsk_BSvol = zeros(length(spx_df.TimeToExpiration),1);
spx_df.putBid_BSvol = zeros(length(spx_df.TimeToExpiration),1);
spx_df.putAsk_BSvol = zeros(length(spx_df.TimeToExpiration),1);
for i = 1: length(spx_df.TimeToExpiration)
    % Find BS implied vols for the bid/ask calls and puts
    T = spx_df.TimeToExpiration(i);
    K = spx_df.Strike(i);
    QT = discountData_df.QT(discountData_df.T == T);
    BT = discountData_df.BT(discountData_df.T == T);
%     % Call bid
%     spx_df.callBid_BSvol(i) = fzero(@(BSvol) BScall(T,K,S0,BSvol,QT, BT) ...
%         - spx_df.CallBidPrice(i),0.1);
%     % Call ask
%     spx_df.callAsk_BSvol(i) = fzero(@(BSvol) BScall(T,K,S0,BSvol,QT, BT) ...
%         - spx_df.CallAskPrice(i),0.1);
%     % Put bid
%     spx_df.putBid_BSvol(i) = fzero(@(BSvol) BSput(T,K,S0,BSvol,QT, BT) ...
%         - spx_df.PutBidPrice(i),0.1);
%     % Put ask
%     spx_df.putAsk_BSvol(i) = fzero(@(BSvol) BSput(T,K,S0,BSvol,QT, BT) ...
%         - spx_df.PutAskPrice(i),0.1); 
    % fzero with tolerance set:
%     options = optimset;
%     opt = optimset(options,'TolX',1e-25);
%     % Call bid
%     spx_df.callBid_BSvol(i) = fzero(@(BSvol) BScall(T,K,S0,BSvol,QT, BT) ...
%         - spx_df.CallBidPrice(i),0.1,opt);
%     % Call ask
%     spx_df.callAsk_BSvol(i) = fzero(@(BSvol) BScall(T,K,S0,BSvol,QT, BT) ...
%         - spx_df.CallAskPrice(i),0.1,opt);
%     % Put bid
%     spx_df.putBid_BSvol(i) = fzero(@(BSvol) BSput(T,K,S0,BSvol,QT, BT) ...
%         - spx_df.PutBidPrice(i),0.1,opt);
%     % Put ask
%     spx_df.putAsk_BSvol(i) = fzero(@(BSvol) BSput(T,K,S0,BSvol,QT, BT) ...
%         - spx_df.PutAskPrice(i),0.1,opt); 
    % fminsearch:
    options = optimset('TolX', 1e-10);
    % Call bid
    spx_df.callBid_BSvol(i) = fminsearch(@(BSvol) abs(BScall(T,K,S0,BSvol,QT, BT) ...
        - spx_df.CallBidPrice(i)),0.1,options);
    % Call ask
    spx_df.callAsk_BSvol(i) = fminsearch(@(BSvol) abs(BScall(T,K,S0,BSvol,QT, BT) ...
        - spx_df.CallAskPrice(i)),0.1,options);
    % Put bid
    spx_df.putBid_BSvol(i) = fminsearch(@(BSvol) abs(BSput(T,K,S0,BSvol,QT, BT) ...
        - spx_df.PutBidPrice(i)),0.1,options);
    % Put ask
    spx_df.putAsk_BSvol(i) = fminsearch(@(BSvol) abs(BSput(T,K,S0,BSvol,QT, BT) ...
        - spx_df.PutAskPrice(i)),0.1,options);  

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
hold on
plot(spx_df.logStrike(filter90), spx_df.putBid_BSvol(filter90), ".");
hold on
plot(spx_df.logStrike(filter90), spx_df.putAsk_BSvol(filter90), ".");
hold off
xlabel("Log strike")
ylabel("BS implied vol")
legend(["Call bid", "Call ask", "Put bid", "Put ask"])
title("SPX: maturity 90 days (2022-4-1)")

%% 

% Define objective function for optimisation problem
% Columns of filteredOptionData:
% TimeToExpiration, Strike, CallMktPrice, PutMktPrice, TotalImplVar, 
% logStrike, CallBidPrice, CallAskPrice, PutBidPrice, PutAskPrice
%function res = CalibrationObjFnc(filteredOptionData, ) 