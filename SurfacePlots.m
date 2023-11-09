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

% Interpolate to get thetaT values at different maturities
thetaT_fnc = @(T) interp1(T_maturities,discountData_df.TotImplVar,T, 'linear');

for i = 1:length(T)
    T_i = T(i);
    %ks = spx_df.logStrike(spx_df.TimeToExpiration == T_i);
    thetaT = thetaT_fnc(T_i);
        impliedVolsHST(i,:) = SSVIimpliedVolatility_Heston(thetaT, T_i, k, ...
            calibration_params_HST.rho, calibration_params_HST.lambda);

        impliedVolsPWR(i,:) = SSVIimpliedVolatility_PowerLaw(thetaT, T_i, k, ...
            calibration_params_PWR.rho, calibration_params_PWR.eta, ...
            calibration_params_PWR.gamma);        

        impliedVolsSPX(i,:) = SSVIimpliedVolatility(thetaT, T_i, k, ...
                calibration_params_SPX.rho, calibration_params_SPX.eps, ...
                calibration_params_SPX.gamma1, calibration_params_SPX.gamma2, ...
                calibration_params_SPX.beta1, calibration_params_SPX.beta2);
    
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
hold off; % Release the hold on the current axes
