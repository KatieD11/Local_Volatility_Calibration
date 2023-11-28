% Model comparison 
% Uses results from the ImpliedVolatility script
clear; clc;

% Select dataset
dataset = "spx_20220401";

% Load option datasets
spx_df=readtable("Data_prep/Data/"+dataset+"_filtered_optionDataWithImplVol.csv");
discountData_df=readtable("Data_prep/Data/"+dataset+"_discountData.csv");
bid_ask_spread=readtable("Data_prep/Data/"+dataset+"_option_bid_ask_spread.csv");

% Load SSVI implied volatility results for each model
SSVIvols_withWeights_data = load('Calibration/Calibration_results/'+dataset+'_SSVIimpliedVols_with_weights.mat');
SSVIvols_withWeights_cellArray = SSVIvols_withWeights_data.SSVI_vols_cellArray;
SSVIvols_freeParams_data = load('Calibration/Calibration_results/'+dataset+'_SSVIimpliedVols_with_free_params.mat');
SSVIvols_freeParams_cellArray = SSVIvols_freeParams_data.SSVI_vols_cellArray;
SSVIvols_heston_data = load('Calibration/Calibration_results/'+dataset+'_SSVIimpliedVols_heston.mat');
SSVIvols_heston_cellArray = SSVIvols_heston_data.SSVI_vols_cellArray;
SSVIvols_powerLaw_data = load('Calibration/Calibration_results/'+dataset+'_SSVIimpliedVols_powerLaw.mat');
SSVIvols_powerLaw_cellArray = SSVIvols_powerLaw_data.SSVI_vols_cellArray;
%% Plot implied volatility results for a particular maturity T
% Choose a time to maturity in days
%Tn_days = 90;
Tn_days = 17; % [17, 19, 21, 24, 26, 28, 31, 35, 42, 49 ...,259, 273 ]

% Get the time to expiration in days (approx)
spx_df.T_days = round(spx_df.TimeToExpiration*365);
discountData_df.T_days = round(discountData_df.T*365);
bid_ask_spread.T_days = round(bid_ask_spread.TimeToExpiration*365);

% Strike values for the maturity
ks = bid_ask_spread.logStrike(bid_ask_spread.T_days == Tn_days);

filter = spx_df.T_days == Tn_days;
figure(1)
plot(ks, spx_df.callBid_BSvol(filter), ".", 'color', [.5 .5 .5]);
hold on
plot(ks, spx_df.callAsk_BSvol(filter), ".", 'color', [.5 .5 .5]);
hold on
plot(ks, spx_df.putBid_BSvol(filter), ".",'color', [.5 .5 .5]);
hold on
plot(ks, spx_df.putAsk_BSvol(filter), ".",'color', [.5 .5 .5]);
hold on
% plot(ks, SSVIvols_withWeights_cellArray{discountData_df.T_days == Tn_days}, "-g", "LineWidth",1);
% hold on
plot(ks, SSVIvols_heston_cellArray{discountData_df.T_days == Tn_days}, "-b", "LineWidth",1);
hold on
plot(ks, SSVIvols_powerLaw_cellArray{discountData_df.T_days == Tn_days}, "-m", "LineWidth",1);
hold on
plot(ks, SSVIvols_freeParams_cellArray{discountData_df.T_days == Tn_days}, "-", "LineWidth",1, 'color',"#00d400");
hold on
plot(bid_ask_spread.logStrike(bid_ask_spread.T_days == Tn_days), ...
    bid_ask_spread.sigma_target(bid_ask_spread.T_days== Tn_days), ...
    "-c", "LineWidth",1);
hold off
xlabel("Log strike")
ylabel("Implied volatility")
% legend(["Call bid", "Call ask", "Put bid", "Put ask", ...
%     "SSVI (with weights)", "SSVI (free parameters)", ...
%     "SSVI (heston-like)", "SSVI (power-law)","Target"])
legend(["Call bid", "Call ask", "Put bid", "Put ask", ...
    "SSVI: heston-like", "SSVI: power-law", ...
    "SSVI: SPX", "Target"])
title("Implied volatilities for maturity " +Tn_days+ " days")
%% Plot MAPEs (Mean absolute percentage errors) vs maturity
% Load MAPEs for each model
mapes_withWeights=readtable("Calibration/Calibration_results/"+dataset+"_mapes_with_weights.csv");
mapes_freeParams=readtable("Calibration/Calibration_results/"+dataset+"_mapes_with_free_params.csv");
mapes_heston=readtable("Calibration/Calibration_results/"+dataset+"_mapes_heston.csv");
mapes_powerLaw=readtable("Calibration/Calibration_results/"+dataset+"_mapes_powerLaw.csv");

T_maturities = mapes_withWeights.T;

% Compute the total open interest weight for each maturity
total_open_int_weight = zeros(length(T_maturities),1);
for i=1:length(T_maturities)
    Ti = T_maturities(i);
    total_open_int_weight(i) = sum(bid_ask_spread.open_interest_weight(...
        bid_ask_spread.TimeToExpiration == Ti));
end

figure(2)
% MAPEs plots
% plot(T_maturities*365, mapes_withWeights.mape, "xg", "LineWidth",2)
% hold on
plot(T_maturities*365, mapes_heston.mape, "xb", "LineWidth",2)
hold on
plot(T_maturities*365, mapes_powerLaw.mape, "xm", "LineWidth",2)
hold on
plot(T_maturities*365, mapes_freeParams.mape, "x", "LineWidth",2, "MarkerEdgeColor","#00d400")
hold on
ylabel("MAPE (%)")
% Open interest weight plot
yyaxis right;
bar(T_maturities*365, total_open_int_weight*100, 'FaceColor', ...
    [0.5, 0.5, 0.5], 'EdgeColor',[0.5, 0.5, 0.5], ...
    'FaceAlpha', 0.2, 'EdgeAlpha', 0.5)
hold off
title("Mean absolute percentage errors of SSVI implied volatilities")
xlabel("Maturity (in days)")
ylabel("Open interest weight (%)")
set(gca, 'ycolor', [0.5, 0.5, 0.5]);  % Set y-axis color to gray
% legend(["SSVI (with weights)", "SSVI (free parameters)", ...
%     "SSVI (heston-like)", "SSVI (power-law)"])
legend(["SSVI: heston-like", "SSVI: power-law", "SSVI: SPX"])
