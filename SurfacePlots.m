% Volatility surface plots 
clear; clc;
addpath('./Data_prep');

% Select files
dataset = "spx_20220401";

spx_df=readtable("Data_prep/Data/"+dataset+"_filtered_optionDataWithImplVol.csv");
discountData_df=readtable("Data_prep/Data/"+dataset+"_discountData.csv");
bid_ask_spread=readtable("Data_prep/Data/"+dataset+"_option_bid_ask_spread.csv");
calibration_params_SPX=readtable("Calibration/Calibration_results/"+dataset+"_calibration_params_"+"with_free_params"+".csv");
calibration_params_HST=readtable("Calibration/Calibration_results/"+dataset+"_calibration_params_"+"heston"+".csv");
calibration_params_PWR=readtable("Calibration/Calibration_results/"+dataset+"_calibration_params_"+"powerLaw"+".csv");

S0 = 4545.86; % spx_20220401
%% Compute implied vols
ks=bid_ask_spread.logStrike;
T_maturities = discountData_df.T; 

T =  linspace(min(T_maturities), max(T_maturities), 100);
k = linspace(min(ks), max(ks), 100);

impliedVolsSPX = zeros(length(T), length(k));
impliedVolsHST = zeros(length(T), length(k));
impliedVolsPWR = zeros(length(T), length(k));

totImplVarHST = zeros(length(T), length(k)); % Plot for report (See ImpliedVolatilies file for testing the condition)
totImplVarPWR = zeros(length(T), length(k));
totImplVarSPX = zeros(length(T), length(k));

% Interpolate to get thetaT values at different maturities
thetaT_fnc = @(T) interp1(T_maturities,discountData_df.TotImplVar,T, 'linear');

for i = 1:length(T)
    T_i = T(i);
    %ks = spx_df.logStrike(spx_df.TimeToExpiration == T_i);
    thetaT = thetaT_fnc(T_i);
        impliedVolsHST(i,:) = SSVIimpliedVolatility_Heston(thetaT, T_i, k, ...
            calibration_params_HST.rho, calibration_params_HST.lambda);
        totImplVarHST(i,:) =impliedVolsHST(i,:).^2*T_i;

        impliedVolsPWR(i,:) = SSVIimpliedVolatility_PowerLaw(thetaT, T_i, k, ...
            calibration_params_PWR.rho, calibration_params_PWR.eta, ...
            calibration_params_PWR.gamma);  
        totImplVarPWR(i,:) =impliedVolsPWR(i,:).^2*T_i;

        impliedVolsSPX(i,:) = SSVIimpliedVolatility(thetaT, T_i, k, ...
                calibration_params_SPX.rho, calibration_params_SPX.eps, ...
                calibration_params_SPX.gamma1, calibration_params_SPX.gamma2, ...
                calibration_params_SPX.beta1, calibration_params_SPX.beta2);
        totImplVarSPX(i,:) =impliedVolsSPX(i,:).^2*T_i;
    
end

%% Plots
[k_axis, T_axis] = meshgrid(k,T); 

figure(1)
surf(k_axis, T_axis, impliedVolsSPX)
xlabel("Log-strike")
ylabel("Maturity")
zlabel("Implied volatility")
title("SPX-fit surface")
saveas(gcf,"Calibration/Calibration_results/"+dataset+"SPX_surf"+".png")

figure(2)
surf(k_axis, T_axis, impliedVolsPWR)
xlabel("Log-strike")
ylabel("Maturity")
zlabel("Implied volatility")
title("Power-law surface")
saveas(gcf,"Calibration/Calibration_results/"+dataset+"PWR_surf"+".png")

figure(3)
surf(k_axis, T_axis, impliedVolsHST)
xlabel("Log-strike")
ylabel("Maturity")
zlabel("Implied volatility")
title("Heston-like surface")
saveas(gcf,"Calibration/Calibration_results/"+dataset+"HST_surf"+".png")
%% Combined plot
figure(4);

% Plot the first surface (SPX)
surf(k_axis, T_axis, impliedVolsHST, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on
surf(k_axis, T_axis, impliedVolsPWR, 'FaceColor', 'magenta', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
surf(k_axis, T_axis, impliedVolsSPX, 'FaceColor', 'green', 'EdgeColor', 'none');
xlabel("Log-strike");
ylabel("Maturity");
zlabel("Implied volatility");
%legend(["Heston-like", "Power-law", "SPX-fit"])
%hold off; % Release the hold on the current axes

% Add a vertical plane
y_value = 0.2466;
%yline(y_value, 'Color', 'red', 'LineWidth', 2, 'Label', 'Y = 0.2466');
xPlane = [-0.8 0.6 0.6 -0.8];      % X coordinates of plane corners, ordered around the plane
%yPlane1 = xPlane+0.5;      % Corresponding y coordinates for plane 1
yPlane1 = [y_value y_value y_value y_value];
zPlane = [0.1 0.1 0.6 0.6];  % Z coordinates of plane corners
                 % Add to existing plot
patch(xPlane, yPlane1, zPlane,'k', 'FaceAlpha', 0.3);  % Plot plane 1
hold off
%% Implied variance plot (for HST surface)
%totImplVarHST
figure(5)
numLines = length(T);
colorMap = winter(numLines);
for i = 1:length(T)
    T_i = T(i);
    color = colorMap(i, :);
    plot(k, totImplVarHST(i,:), "-", "LineWidth",1, 'Color', color);
    hold on
end
hold off
%legend(string(T))
title("Total implied variance plot for the Heston-like surface")
xlabel("Log strike")
ylabel("Total implied variance")
exportgraphics(gcf,'Calibration/Calibration_results/'+dataset+'totImplVarHST.pdf','ContentType','vector')
%% Implied variance plot (for SPX-fit surface)
figure(6)
for i = 1:numLines
    T_i = T(i);
    color = colorMap(i, :); 
    plot(k, totImplVarSPX(i, :), "-", "LineWidth", 1, 'Color', color);
    hold on
end
hold off
%legend(string(T))
title("Total implied variance plot for the SPX-fit surface")
xlabel("Log strike")
ylabel("Total implied variance")
exportgraphics(gcf,'Calibration/Calibration_results/'+dataset+'totImplVarSPX.pdf','ContentType','vector')
%% Implied variance plot (for PWR surface)
figure(7)
for i = 1:length(T)
    T_i = T(i);
    color = colorMap(i, :);
    plot(k, totImplVarPWR(i,:), "-", "LineWidth",1, 'Color', color);
    hold on
end
hold off
%legend(string(T))
title("Total implied variance plot for the power-law surface")
xlabel("Log strike")
ylabel("Total implied variance")
exportgraphics(gcf,'Calibration/Calibration_results/'+dataset+'totImplVarPWR.pdf','ContentType','vector')
%% Check subset of the smaller maturities
figure(8)
colorMap = winter(18);
for i = 1:18
    T_i = T(i);
    color = colorMap(i, :);
    plot(k, totImplVarPWR(i,:), "-", "LineWidth",1, 'Color', color);
    hold on
end
hold off
%legend(string(T))
title("Total implied variance plot for the power-law surface, maturities 17 to 80 days")
xlabel("Log strike")
ylabel("Total implied variance")
exportgraphics(gcf,'Calibration/Calibration_results/'+dataset+'totImplVarPWR2.pdf','ContentType','vector')
%% Check conditions for absence of butterfly arbitrage:
    % i)  θϕ(θ) (1 + |ρ|)  <  4,  for  all  θ  >  0
    % ii) θϕ(θ)^2  (1 + |ρ|)  ≤  4,  for  all  θ  >  0
theta = discountData_df.TotImplVar;
%theta = (min(discountData_df.TotImplVar):1:20)'; % Test extreme values

% Power-law surface
phi_PWR = calibration_params_PWR.eta./(theta.^calibration_params_PWR.gamma.*(1+theta).^(1-calibration_params_PWR.gamma));
C1 = theta.*phi_PWR*(1+abs(calibration_params_PWR.rho));
C2 = theta.*phi_PWR.^2*(1+abs(calibration_params_PWR.rho));

figure(9)
yline(4)
hold on
plot(theta, C1)
plot(theta, C2)
hold off
title("Power-law surface: Absence of butterfly arbitrage conditions")
xlabel("θ")
ylabel("Value of C1/C2")
legend(["","C1: θ ϕ(θ) (1 + |ρ|)", "C2: θ ϕ(θ)^2  (1 + |ρ|)"])
%exportgraphics(gcf,'Calibration/Calibration_results/'+dataset+'butterflyPWR.png')

% Heston-like surface
phi_HST = 1./(calibration_params_HST.lambda*theta).* ... 
        (1-(1-exp(-calibration_params_HST.lambda*theta))./...
        (calibration_params_HST.lambda*theta));
C1 = theta.*phi_HST*(1+abs(calibration_params_HST.rho));
C2 = theta.*phi_HST.^2*(1+abs(calibration_params_HST.rho));

figure(10)
yline(4)
hold on
plot(theta, C1)
plot(theta, C2)
hold off
title("Heston-like surface: Absence of butterfly arbitrage conditions")
xlabel("θ")
ylabel("Value of C1/C2")
legend(["","C1: θ ϕ(θ) (1 + |ρ|)", "C2: θ ϕ(θ)^2  (1 + |ρ|)"])
%exportgraphics(gcf,'Calibration/Calibration_results/'+dataset+'butterflyHST.png')
%exportgraphics(gcf,'Calibration/Calibration_results/'+dataset+'butterflyHST3.png')

figure(11)
plot(theta, C1)
hold on
plot(theta, C2)
hold off
title("Heston-like surface: Absence of butterfly arbitrage conditions")
xlabel("θ")
ylabel("Value of C1/C2")
legend(["C1: θ ϕ(θ) (1 + |ρ|)", "C2: θ ϕ(θ)^2  (1 + |ρ|)"])
%exportgraphics(gcf,'Calibration/Calibration_results/'+dataset+'butterflyHST2.png')

% SPX-fit surface
eta = 2.016048*exp(calibration_params_SPX.eps);
phi_SPX = eta./(theta.^calibration_params_SPX.gamma1.*...
    (1+calibration_params_SPX.beta1*theta).^calibration_params_SPX.gamma2.* ...
        (1+calibration_params_SPX.beta2*theta).^...
        (1-calibration_params_SPX.gamma1-calibration_params_SPX.gamma2));
C1 = theta.*phi_SPX*(1+abs(calibration_params_SPX.rho));
C2 = theta.*phi_SPX.^2*(1+abs(calibration_params_SPX.rho));

figure(12)
yline(4)
hold on
plot(theta, C1)
plot(theta, C2)
hold off
title("SPX-fit surface: Absence of butterfly arbitrage conditions")
xlabel("θ")
ylabel("Value of C1/C2")
legend(["","C1: θ ϕ(θ) (1 + |ρ|)", "C2: θ ϕ(θ)^2  (1 + |ρ|)"])
%exportgraphics(gcf,'Calibration/Calibration_results/'+dataset+'butterflySPX.pdf','ContentType','vector')
%exportgraphics(gcf,'Calibration/Calibration_results/'+dataset+'butterflySPX.png')


