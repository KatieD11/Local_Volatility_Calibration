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
%dT=0.001; dk=0.001;
dT=0.001/100; dk=0.001*10;
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
%% Pricing using theta method
% Choose a time to maturity in days
Tn_days = 90; 
%Tn_days = 140; 
T = Tn_days/365;

% Approximate qT (from QT) as average dividend yield over the period
discountData_df.qT = -log(discountData_df.QT)./discountData_df.T;
% Compute FT in the dataset
discountData_df.FT = S0*(discountData_df.QT)./(discountData_df.BT);

% Use interpolation to estimate r(T), q(T), FT, BT for maturities not in the data set
r_fnc = @(T) interp1(discountData_df.T,discountData_df.rT,T);
q_fnc = @(T) interp1(discountData_df.T,discountData_df.qT,T);
F_fnc = @(T) interp1(discountData_df.T,discountData_df.FT,T);
B_fnc = @(T) interp1(discountData_df.T,discountData_df.BT,T); %DF

% Compute K given a log-strike k(T) and k given K
K_fnc = @(k, T) F_fnc(T).*exp(k);
k_fnc = @(K,T) log(K./F_fnc(T));

% Use theta-method to estimate put option prices

% Option related parameters
Ks = spx_df.Strike(spx_df.T_days == Tn_days); % strikes in data set
K = Ks(40); % Select a strike to compute option value for
%K = Ks(150);
r_ave = r_fnc(T);
q_ave = q_fnc(T);

% Mesh related parameters
Smin=S0*0.8; Smax=S0*1.2; 
%M=80; N=160; 
M=5; N=100; % local var doesn't work for too many time steps
theta=0.5;
dS=(Smax - Smin)/N;
n=1:(N-1);
dTau=T/M;
m=0:M;

% tstart = 0.05;
% M = 5;
% dTau=(T-tstart)/M;
% m=0:M;

% Compute the log equivalents of S-values (to use in local vol fnc)
s_min = k_fnc(Smin,T); s_max = k_fnc(Smax,T);
ds = k_fnc(dS,T);

% Initial condition
Uinit = @(S) max(K-S,0); % Terminal condition (reverse time transformation)
% Boundary conditions
U0 = @(S,tau) K*exp(-(r_ave-q_ave)*tau) - S;
Uinf = @(S,tau) 0; 

% Theta method:
% Create matrices
D1 = diag(Smin/dS + (1:N-1));
D2 = D1^2;
T1 = diag(ones(N-2,1), 1) - diag(ones(N-2,1), -1);
T2 = diag(ones(N-2,1), 1) + diag(ones(N-2,1), -1) - ...
    2*diag(ones(N-1,1));

% Boundary condition vectors (in (N-1)x(M+1) matrix)
tm = T - m*dTau;
%tm = tstart + T - m*dTau;
b = zeros(N-1, M+1);
b(1,:) = 0.5*dTau*(Smin/dS+1)*(local_var(s_min+ds,tm).*(Smin/dS+1) ...
    -(r_ave - q_ave)).*U0(Smin, m*dTau);
b(end,:) = 0.5*dTau*(Smax/dS-1)*(local_var(s_max-ds,tm).*(Smax/dS-1) ...
    +(r_ave-q_ave))*Uinf(Smax, m*dTau);
% Replace NaN values with 0
b(isnan(b))=0;

% Solution U
U = zeros(N-1, M+1);
U(:,1) = Uinit(Smin+n*dS)'; % Initial condition

for i=1:M
    t = tm(i);
    % Since local vol fnc is time dependent for this model,
    % sigma_matrix, Fm and Gm+1 are computed inside the loop
    %local_var(s_min+(n*ds),tm)
    sigma_matrix = diag(local_var(s_min+(n*ds),t));
    if (max(isnan(sigma_matrix),[],'all')==1)
        disp("Local var is nan at t="+t)
        %sigma_matrix(isnan(sigma_matrix))=0;
    end
    F = (1-(r_ave-q_ave)*dTau)*eye(N-1) + 0.5*(r_ave-q_ave)*dTau*D1*T1 + 0.5*sigma_matrix*dTau*D2*T2;
    G = (1+(r_ave-q_ave)*dTau)*eye(N-1) - 0.5*(r_ave-q_ave)*dTau*D1*T1 - 0.5*sigma_matrix*dTau*D2*T2;

    U(:,i+1) = (theta*G + (1-theta)*eye(N-1))\ ...
        (((1-theta)*F + theta*eye(N-1))*U(:,i) + ...
        (1-theta)*b(:,i) + theta*b(:,i+1));
end

% Surface plot of the result (in terms of S and tau)
[tau,Ss] = meshgrid(m*dTau,Smin+(n*dS));
figure(1);
surf(tau,Ss,U);
%t = T-tau;
%surf(t,Ss,U);
%xlabel("t")
xlabel("tau")
ylabel("S")
zlabel("V")
title("Put option pricing surface")

% Price estimate
U(Ss==S0 & tau==T)
%U(Ss==S0 & tau==T-tstart)
spx_df.PutMktPrice(spx_df.T_days == Tn_days & spx_df.Strike == K)

% Plots of prices at t=0
figure(2)
plot(Smin+(n*dS), U(tau==T), ".k")
title("Put price estimates at t=0")
xlabel("S_0")
ylabel("Price")
legend(["Closed form", "Estimate"])

%% 

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

