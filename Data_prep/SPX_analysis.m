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

% Iterate through each time to maturity Ti and 
% find the ATM total implied variance for FTi
% θTi =σ^2(T,FTi)* Ti
total_impl_vars = zeros(length(T_vals),1);
impl_vols = zeros(length(T_vals),1); % BS implied volatilities
r=5/100; % risk-free rate
for i=1:length(T_vals) 
    % var = TotalImpliedVariance(Ki, Ti, r, S0, optionData);
    [total_impl_vars(i), impl_vols(i)] = TotalImpliedVariance(FT(i), T_vals(i), r, S0, optionData);
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
optionData.k = k_vals;

% Filter option data:
% All options with |kj|/sqrt(θTi) > 3.5 are censored from the data set
filtered_optionData = optionData;
for i = 1:length(T_vals)
    filter = (filtered_optionData.TimeToExpiration == T_vals(i)) & ...
        (abs(filtered_optionData.k)/sqrt(total_impl_vars(i)) > 3.5);
    filtered_optionData(filter,:) =[];
end
%writetable(filtered_optionData, "Data/spx_quotedata20220401_filtered_optionData.csv")
