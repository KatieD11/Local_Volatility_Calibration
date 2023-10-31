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

% Set up total implied variance w as a function of T and k
w = @(k,T) SSVItotalImpliedVariance(discountData_df, T, k, ...
            calibration_params.rho, calibration_params.eps, ...
            calibration_params.gamma1, calibration_params.gamma2, ...
            calibration_params.beta1, calibration_params.beta2);

% % Central finite difference estimates for the partial derivatives
% dT=0.001; dk=0.001;
% %dT=0.0001; dk=0.0001;
% dwdT = @(k, T) (w(k, T+dT) - w(k, T-dT))/(2*dT);
% dwdk = @(k, T) (w(k+dk, T) - w(k-dk, T))/(2*dk);
% d2wdk2 = @(k, T) (w(k+2*dk, T) - 2*w(k, T) + w(k-2*dk, T))/(4*dk^2);
% 
% % Local volatility function (implemented as vol^2)
% local_var = @(k,T) dwdT(k, T)./(1 - k./w(k,T) .* dwdk(k, T) + ...
%     1/4 * (-1/4 - 1./w(k,T) + k.^2./w(k,T).^2).*(dwdk(k, T)).^2 ...
%     + 1/2* d2wdk2(k,T));

%% Test local vol function
dT=0.001; dk=0.001;

%T = T_min-20*dT; %0.25;
%T = 0.022; dT=0.001;
dT=0.0001; T = 0.021;
k1=-0.3;
k2=0.2;

vol = LocalVolFD(dT, dk, w, k1, T)
vol = LocalVolFD(dT, dk, w, k2, T)
%% 
% Find vols for a range of times to maturity
% k=0.2;
% T = discountData_df.T;
% vols = LocalVolFD(dT, dk, w, k, T)
% 
% % Extrapolate values for a T < T_min
% Te = 0:0.01:T_max;
% vol_extrap = @(Te) interp1(T,vols,Te, 'linear', 'extrap');
% 
% plot(T,vols, "-b")
% hold on
% plot(T,vols, ".b")
% hold on
% plot(Te,vol_extrap(Te), "or")
% hold off
% title("Local vol as a function of T for log-strike: "+k)
% xlabel("Time to maturity")
% ylabel("Local vol")

%% Test: local vol surface plot
%[X,Y] = meshgrid(-0.3:0.005:0.3, T_min:0.005:T_max);
[X,Y] = meshgrid(-0.3:0.005:0.3, 0.021:0.005:T_max);
%Z = sqrt(local_var(X,Y));
Z = sqrt(LocalVolFD(dT, dk, w, X, Y));
figure(1)
surf(X,Y,Z)
xlabel("log-strike")
ylabel("T")
zlabel("Local volatility")
title("Local vol surface")
%% Test: Plot results for a particular maturity T
% Choose a time to maturity in days
Tn_days = 90; % 17; % 60; % 35; %21;
%Tn_days = 17; % [17, 19, 21, 24, 26, 28, 31, 35, 42, 49 ..., 259, 273]
T = (Tn_days-0)/365;

% Calculate SSVI implied vols & local vols
%k_set = -0.4:0.01:0.2;
k_set = -0.2:0.01:0.2;
SSVI_vols = zeros(length(k_set),1);
local_vols = zeros(length(k_set),1);
for i = 1:length(k_set)
    wi = w(k_set(i), T);
    SSVI_vols(i) = sqrt(wi/T);
    local_vols(i) = LocalVolFD(dT, dk, w, k_set(i), T); 
end

% Get the time to expiration in days (approx)
spx_df.T_days = round(spx_df.TimeToExpiration*365);
discountData_df.T_days = round(discountData_df.T*365);
bid_ask_spread.T_days = round(bid_ask_spread.TimeToExpiration*365);

% Strike values for the maturity
ks = bid_ask_spread.logStrike(bid_ask_spread.T_days == Tn_days);

filter = spx_df.T_days == Tn_days;
figure(2)
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
%% MC Pricing
rng(0);

% Approximate qT (from QT) as average dividend yield over the period
discountData_df.qT = -log(discountData_df.QT)./discountData_df.T;
% Compute FT in the dataset
discountData_df.FT = S0*(discountData_df.QT)./(discountData_df.BT);

% Use interpolation to estimate r(T), q(T), FT, BT for maturities not in the data set
r = @(T) interp1(discountData_df.T,discountData_df.rT,T, 'linear', 'extrap');
q = @(T) interp1(discountData_df.T,discountData_df.qT,T, 'linear', 'extrap');
F = @(T) interp1(discountData_df.T,discountData_df.FT,T, 'linear', 'extrap');
B = @(T) interp1(discountData_df.T,discountData_df.BT,T, 'linear', 'extrap'); %DF
Q = @(T) interp1(discountData_df.T,discountData_df.QT,T, 'linear', 'extrap'); 

% Compute K given a log-strike k(T) and k given K
K = @(k, T) F(T).*exp(k);
k = @(K,T) log(K./F(T));

% Use MC simulation to estimate put option prices

% MC parameters
Ks = spx_df.Strike(spx_df.T_days == Tn_days); % strikes in data set
%n = 50000; % # MC simulations
n = 1000000;
t_start = 0.01; dT=0.0001;
%t_start = 0.001;
M = 100;
m=0:(M-1);

t = t_start + m*(T-t_start)/M;
dt = t(1:end) - [0, t(1:end-1)];

r_ave = r(T); q_ave = q(T);

% Generate Sobol sequences
q = qrandstream('sobol',M);
X = rand(q,n,M);

% Simulate stock price paths (Euler-Maruyama)
St = S0;
Xt = log(St);
for j = 1:M
    dtj = dt(j);
%     Z1 = randn(1,n);
%     Z2 = -Z1;
%     Z = [Z1, Z2];
%     Wt = sqrt(dtj) * Z;
    
    % Using sobol sequences
    Z1 = norminv(X(:,j));
    Z2 = -Z1;
    Z = [Z1; Z2]';
    Wt = sqrt(dtj) * Z;

    % Local_vol at t, St 
    % set initial t0 time to 0.001 (to approx 0, since vol_t at 0 is undefined)
    if (j==1)
        rj = -1/dtj*(log(B(t(j))) - log(B(0.001)));
        qj = -1/dtj*(log(Q(t(j))) - log(Q(0.001)));
        vol_t = LocalVolFD(dT, dk, w, k(St,0.001), 0.001);
    else
        rj = -1/dtj*(log(B(t(j))) - log(B(t(j-1))));
        qj = -1/dtj*(log(Q(t(j))) - log(Q(t(j-1))));            
        vol_t = LocalVolFD(dT, dk, w, k(St,t(j-1)), t(j-1));
    end  

    Xt = Xt + (r_ave - q_ave - 0.5 * vol_t.^2) * dtj + vol_t .* Wt;
    %Xt = Xt + (rj - qj - 0.5 * vol_t.^2) * dtj + vol_t .* Wt;
    St = exp(Xt);
end 
 
% Estimate prices for each strike at maturity
p_hat = zeros(length(Ks),1);
p_mkt = zeros(length(Ks),1);
for i = 1: length(Ks)
    % Discounted payoff function
    f_put = B(T)*max(Ks(i)-St,0);
    % Price estimate
    p_hat(i) = mean(f_put); 

    % Store corresponding market price
    p_mkt(i) = mean(spx_df.PutMktPrice(spx_df.Strike == Ks(i) & spx_df.T_days == Tn_days));
end

figure(3)
plot(Ks, p_hat, ".")
hold on
plot(Ks, p_mkt, ".")
hold off
title("MC put price estimates for T ="+Tn_days+" days")
xlabel("Strike (K)")
ylabel("Put price")
legend(["MC estimates", "Market values"])
%% Reverse pricing to get BS implied vols
% Find Black-Scholes implied volatilities (from MC price estimates)
QT = discountData_df.QT(discountData_df.T_days == Tn_days);
BT = discountData_df.BT(discountData_df.T_days == Tn_days);

p_BSvol = zeros(length(p_hat),1);
for i = 1:length(p_hat)
    %call_BSvol = fzero(@(BSvol) BScall(Ti,Ki,S0,BSvol,QT, BT) - call_i,0.1);
    %p_BSvol(i) = fzero(@(BSvol) BSput(T,Ks(i),S0,BSvol,QT, BT) - p_hat(i),0.1); 
    options = optimset('TolFun', 1e-30);
    p_BSvol(i) = fminsearch(@(BSvol) abs(BSput(T,Ks(i),S0,BSvol,QT, BT) - p_hat(i)),0.1, options); 
end

ks = k(Ks,T);

filter = spx_df.T_days == Tn_days;
figure(4)
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
plot(ks, p_BSvol, "xk", "LineWidth",1);
hold off
xlabel("Log strike")
ylabel("BS implied vol")
legend(["Put bid", "Put ask", "SSVI", "Target", "Put MC estimates"])
title("Implied volatilities for maturity " +Tn_days+ " days, "+ regexprep(dataset,'_',' '))
