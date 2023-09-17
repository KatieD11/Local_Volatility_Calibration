% S&P500 data analysis 2022–04–01
clear; clc;
% Load data: 
% Closing prices for all European call and put options 
% for various strikes and maturities
spx_df=readtable("Data/spx_quotedata20220401_all.xlsx");
%% 
% Spot price at closing on 2022–04–01 
S0 = 4545.86;
t0 = datetime('2022-04-01', 'Format', 'yyyy-MM-dd');

% Only consider options with a time to maturity
% >= 14 days and <= 400 days
dateFormat = 'eee MMM dd yyyy';
ExpirationDate = datetime(spx_df.ExpirationDate, 'InputFormat', dateFormat);
timeFilter = days(ExpirationDate - t0)>=14 & ...
    days(ExpirationDate - t0)<=400;
spx_df_trim = spx_df(timeFilter,:);

% ExpirationDate and Strike: common to both calls and puts
% Bid, Ask and OpenInterest: for call options
% Bid_1, Ask_1 and OpenInterest_1: for put options

% Calculate time to expiration from ExpirationDate (trimmed version)
ExpirationDate_trim = datetime(spx_df_trim.ExpirationDate, 'InputFormat', dateFormat);
TimeToExpiration = years(ExpirationDate_trim - t0);

% Option data
optionData = table;
optionData.TimeToExpiration = TimeToExpiration;
optionData.Strike = spx_df_trim.Strike;
% Estimate market price as midpoint of bid and ask
optionData.CallMktPrice = (spx_df_trim.Bid+spx_df_trim.Ask)/2;
optionData.PutMktPrice = (spx_df_trim.Bid_1+spx_df_trim.Ask_1)/2;

% Store option data
%writetable(optionData, "Data/spx_quotedata20220401_optionData.csv")
%% 
% Find discount factors from dataset
DFs = DiscountFactors(optionData, S0);
BT = DFs(:,1);
QT = DFs(:,2);
% Maturity times
T_vals = unique(optionData.TimeToExpiration);

% Plot discount factors
figure(1)
plot(T_vals, BT, "-b")
hold on
plot(T_vals, QT, "-r")
hold off
legend(["B(T)", "Q(T)"])
xlabel("Maturity (T)")
title("Discount factors")
%% 
% Find the ATM total implied variance θTi for each maturity

% Forward prices 
FT = S0*QT./BT; 
%plot(T_vals, FT)

% Find discount rates from BT 
rT = -1./T_vals.*log(BT);

% Iterate through each time to maturity Ti and 
% find the ATM total implied variance for FTi
% θTi =σ^2(Ti,FTi)* Ti
total_impl_vars = zeros(length(T_vals),1);
impl_vols = zeros(length(T_vals),1); % BS implied volatilities
optionData.TotalImplVar = zeros(length(optionData.TimeToExpiration),1);
for i=1:length(T_vals) 
    [total_impl_vars(i), impl_vols(i)] = TotalImpliedVariance(FT(i), ...
        T_vals(i), QT(i), BT(i), S0, optionData);
    % Add total implied variance to option data table
    optionData.TotalImplVar(optionData.TimeToExpiration == T_vals(i)) ...
    = total_impl_vars(i);

end

figure(2)
plot(T_vals, total_impl_vars)
xlabel("Maturity (T)")
ylabel("θT")
title("ATM total implied variance θT for different maturities")

figure(3)
plot(T_vals, impl_vols)
xlabel("Maturity (T)")
ylabel("BS implied volatility")
title("BS implied volatilities for different maturities")
%% 
% Compute values of kj,  k(T,K)=log(K/FT) 
k_vals = [];
for i = 1:length(T_vals)
    Ks = optionData.Strike(optionData.TimeToExpiration == T_vals(i));
    kj = log(Ks/FT(i));
    k_vals = [k_vals; kj];
end
optionData.logStrike = k_vals;

% Add bid and ask prices back before filtering data set
optionData.CallBidPrice = spx_df_trim.Bid;
optionData.CallAskPrice = spx_df_trim.Ask;
optionData.PutBidPrice = spx_df_trim.Bid_1;
optionData.PutAskPrice = spx_df_trim.Ask_1;

% Filter option data:
% All options with |kj|/sqrt(θTi) > 3.5 are censored from the data set
filter = ((abs(optionData.logStrike)./sqrt(optionData.TotalImplVar)) > 3.5);
filtered_optionData = optionData(~filter,:);

%writetable(filtered_optionData, "Data/spx_quotedata20220401_filtered_optionData.csv")
%% 
% Create a table with discount factors, TotalImplVar, and discount rates
discountData = table;
discountData.T = T_vals;
discountData.BT = BT;
discountData.QT = QT;
discountData.TotImplVar = total_impl_vars;
discountData.rT = rT;
%writetable(discountData, "Data/spx_quotedata20220401_discountData.csv")
