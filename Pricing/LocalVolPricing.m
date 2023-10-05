% Pricing using the local volatility fnc
clear; clc;
addpath('../Data_prep');
% Select files
dataset = "spx_20220401";

% Select parameter set
calibration_set = "with_free_params";

spx_df=readtable("../Data_prep/Data/"+dataset+"_filtered_optionDataWithImplVol.csv");
discountData_df=readtable("../Data_prep/Data/"+dataset+"_discountData.csv");
calibration_params=readtable("../Calibration/Calibration_results/"+dataset+"_calibration_params_"+calibration_set+".csv");
bid_ask_spread=readtable("../Data_prep/Data/"+dataset+"_option_bid_ask_spread.csv");

S0 = 4545.86; % spx_20220401
%% 
%rT = -log(discountData_df.BT)./discountData_df.T
T_min = min(discountData_df.T)
T_max = max(discountData_df.T)

% Set up w as a function of T and k
w = @(k,T) SSVItotalImpliedVariance(discountData_df, T, k, ...
            calibration_params.rho, calibration_params.eps, ...
            calibration_params.gamma1, calibration_params.gamma2, ...
            calibration_params.beta1, calibration_params.beta2);

% Central finite difference estimates for the partial derivatives
dT=0.001; dk=0.001;
%dT=0.0001; dk=0.0001;
dwdT = @(k, T) (w(k, T+dT) - w(k, T-dT))/(2*dT);
dwdk = @(k, T) (w(k+dk, T) - w(k-dk, T))/(2*dk);
d2wdk2 = @(k, T) (w(k+2*dk, T) - 2*w(k, T) + w(k-2*dk, T))/(4*dk^2);

% Local volatility function (implemented as vol^2)
local_var = @(k,T) dwdT(k, T)./(1 - k./w(k,T) .* dwdk(k, T) + ...
    1/4 * (-1/4 - 1./w(k,T) + k.^2./w(k,T).^2).*(dwdk(k, T)).^2 ...
    + 1/2* d2wdk2(k,T));

%% Test: Plot results for a particular maturity T
% Choose a time to maturity in days
Tn_days = 90; T = (Tn_days-0)/365;
%Tn_days = 259; % [17, 19, 21, 24, 26, 28, 31, 35, 42, 49 ..., ]

% Calculate SSVI implied vols & local vols
k_set = -0.4:0.01:0.2;
SSVI_vols = zeros(length(k_set),1);
local_vols = zeros(length(k_set),1);
for i = 1:length(k_set)
    wi = w(k_set(i), T);
    SSVI_vols(i) = sqrt(wi/T);
    local_vols(i) = sqrt(local_var(k_set(i),T));
end

% Get the time to expiration in days (approx)
spx_df.T_days = round(spx_df.TimeToExpiration*365);
discountData_df.T_days = round(discountData_df.T*365);
bid_ask_spread.T_days = round(bid_ask_spread.TimeToExpiration*365);

% Strike values for the maturity
ks = bid_ask_spread.logStrike(bid_ask_spread.T_days == Tn_days);

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
plot(k_set, SSVI_vols, "-g", "LineWidth",1);
hold on
plot(bid_ask_spread.logStrike(bid_ask_spread.T_days == Tn_days), ...
    bid_ask_spread.sigma_target(bid_ask_spread.T_days== Tn_days), ...
    "-c", "LineWidth",1);
hold on
plot(k_set, local_vols, "-k", "LineWidth",1);
hold off
xlabel("Log strike")
ylabel("BS implied vol")
legend(["Call bid", "Call ask", "Put bid", "Put ask", "SSVI", "Target", "Local vol"])
title("Implied volatilities for maturity " +Tn_days+ " days, "+ regexprep(dataset,'_',' '))
%% Pricing
% Choose a time to maturity in days
Tn_days = 90; 
%Tn_days = 140; 
T = Tn_days/365;

% Approximate qT (from QT) as average dividend yield over the period
discountData_df.qT = -log(discountData_df.QT)./discountData_df.T;
% Compute FT in the dataset
discountData_df.FT = S0*(discountData_df.QT)./(discountData_df.BT);

% Use interpolation to estimate r(T), q(T), FT, BT for maturities not in the data set
r = @(T) interp1(discountData_df.T,discountData_df.rT,T);
q = @(T) interp1(discountData_df.T,discountData_df.qT,T);
F = @(T) interp1(discountData_df.T,discountData_df.FT,T);
B = @(T) interp1(discountData_df.T,discountData_df.BT,T); %DF

% Compute K given a log-strike k(T) and k given K
K = @(k, T) F(T).*exp(k);
k = @(K,T) log(K./F(T));

% Use MC simulation to estimate put option prices

% MC parameters
Ks = spx_df.Strike(spx_df.T_days == Tn_days); % strikes in data set
%M = 100; % # time-steps
M = 5;
dt = T/M;
t = dt:dt:T;
n = 50000; % # MC simulations

% Estimate prices for each strike at maturity
p_hat = zeros(length(Ks),1);
p_mkt = zeros(length(Ks),1);
for i = 1: length(Ks)
    ki = k(Ks(i),T); % log-strike
    % Approximate time-average of local_vol^2
    local_vari = @(T) local_var(ki,T);
    %ave_local_var = 1/T * integral(local_vari,0,T)
    N = 2000;
    %ave_local_var = AveLocalVar(local_vari,T,N)
    ave_local_var = local_var(ki,T); % Test just using terminal vol

    % Terminal stock values
    Z = randn(1,n);
    Si = S0*exp((r(T)-q(T)-0.5*ave_local_var)'*T + ...
        sqrt(ave_local_var)'.*sqrt(T).*Z);
    
    % Discounted payoff function
    f_put = B(T)*max(Ks(i)-Si,0);
    % Price estimate
    p_hat(i) = mean(f_put); 

    % Store corresponding market price
    p_mkt(i) = mean(spx_df.PutMktPrice(spx_df.Strike == Ks(i) & spx_df.T_days == Tn_days));
end

figure(2)
plot(Ks, p_hat, ".")
hold on
plot(Ks, p_mkt, ".")
hold off
title("MC put price estimates for T ="+Tn_days+" days")
xlabel("Strike (K)")
ylabel("Put price")
legend(["MC estimates", "Market values"])
%% Reverse pricing to get BS implied vols
% Find Black-Scholes implied volatilities (using interpolated prices)
%call_BSvol = fzero(@(BSvol) BScall(Ti,Ki,S0,BSvol,QT, BT) - call_i,0.1);
QT = discountData_df.QT(discountData_df.T_days == Tn_days);
BT = discountData_df.BT(discountData_df.T_days == Tn_days);

p_BSvol = zeros(length(p_hat),1);
for i = 1:length(p_hat)
    p_BSvol(i) = fzero(@(BSvol) BSput(T,Ks(i),S0,BSvol,QT, BT) - p_hat(i),0.1); 
end

k_set2 = k(Ks,T);

filter = spx_df.T_days == Tn_days;
figure(3)
% plot(ks, spx_df.callBid_BSvol(filter), ".b");
% hold on
% plot(ks, spx_df.callAsk_BSvol(filter), ".r");
% hold on
plot(ks, spx_df.putBid_BSvol(filter), ".",'Color',"#0072BD");
hold on
plot(ks, spx_df.putAsk_BSvol(filter), ".",'Color',"#A2142F");
hold on
plot(k_set, SSVI_vols, "-g", "LineWidth",1);
hold on
plot(bid_ask_spread.logStrike(bid_ask_spread.T_days == Tn_days), ...
    bid_ask_spread.sigma_target(bid_ask_spread.T_days== Tn_days), ...
    "-c", "LineWidth",1);
hold on
plot(k_set2, p_BSvol, "x", "LineWidth",1);
hold off
xlabel("Log strike")
ylabel("BS implied vol")
legend(["Put bid", "Put ask", "SSVI", "Target", "Put estimates"])
title("Implied volatilities for maturity " +Tn_days+ " days, "+ regexprep(dataset,'_',' '))

%% 
% Numerical integration to estimate average local volatility
function ave_local_var = AveLocalVar(local_var_fnc,T,N)
    dt = T / N;
    t = linspace(0, T, N+1);
    integral = 0;
    
    for i = 1:N
        t_i = t(i);
        t_next = t(i+1);
        
        % Use mid-pt of function
        var_t = local_var_fnc((t_i+t_next)/2);
        if isnan(var_t)
            continue
        end
        % Integrate σ(t, St) over the time interval [t_i, t_next]
        integral = integral + var_t * dt;
    end
    
    % Calculate the average local volatility
    ave_local_var = integral / T;

end
