% Implied volatilities
% Compute implied volatilies using SSVI formulation 
% with parameters from calibration
clear; clc;
addpath('./Data_prep');

% Select files
dataset = "spx_20220401";
%calibration_set = "with_weights";
%calibration_set = "without_weights";
calibration_set = "with_free_params";

spx_df=readtable("Data_prep/Data/"+dataset+"_filtered_optionDataWithImplVol.csv");
discountData_df=readtable("Data_prep/Data/"+dataset+"_discountData.csv");
calibration_params=readtable("Calibration_results/"+dataset+"_calibration_params_"+calibration_set+".csv");
bid_ask_spread=readtable("Calibration_results/"+dataset+"_calibration_bid_ask_spread.csv");

S0 = 4545.86; % spx_20220401
%% Compute implied vols for all maturities in dataset
SSVI_vols_cellArray = {}; % Store arrays of different lengths
ks_cellArray = {};

T_maturities = discountData_df.T; 
vol_mapes = zeros(length(T_maturities),1);
for i = 1:length(T_maturities)
    T_i = T_maturities(i);
    %ks = spx_df.logStrike(spx_df.TimeToExpiration == T_i);
    ks = bid_ask_spread.logStrike(bid_ask_spread.TimeToExpiration == T_i);
    thetaT = discountData_df.TotImplVar(discountData_df.T == T_i);
    SSVI_vols = SSVIimpliedVolatility(thetaT, T_i, ks, ...
            calibration_params.rho, calibration_params.eps, ...
            calibration_params.gamma1, calibration_params.gamma2, ...
            calibration_params.beta1, calibration_params.beta2);
    % Store arrays
    ks_cellArray{end+1} = ks;
    SSVI_vols_cellArray{end+1} = SSVI_vols;
    % Calculate RMSE
    target_vols = bid_ask_spread.sigma_target(bid_ask_spread.TimeToExpiration == T_i);
    vol_mapes(i) = mape(SSVI_vols, target_vols);
end
%% Plot MAPEs (Mean absolute percentage errors) vs maturity
figure(1)
plot(T_maturities, vol_mapes, "x", "LineWidth",2)
title("Mean absolute percentage errors of SSVI vols vs target implied vols, "+ regexprep(dataset,'_',' '))
xlabel("Maturity (T)")
ylabel("MAPE (%)")
results_errors = table;
results_errors.T = T_maturities;
results_errors.mape = vol_mapes;
writetable(results_errors, "Calibration_results/"+dataset+"_mapes_"+calibration_set+".csv")
exportgraphics(gcf,'Calibration_results/'+dataset+'mapes_'+calibration_set+'.pdf','ContentType','vector')
%% Plot results for a particular maturity T
% Choose a time to maturity in days
Tn_days = 90;
%Tn_days = 17; % [17, 19, 21, 24, 26, 28, 31, 35, 42, 49 ..., ]

% Get the time to expiration in days (approx)
spx_df.T_days = round(spx_df.TimeToExpiration*365);
discountData_df.T_days = round(discountData_df.T*365);
bid_ask_spread.T_days = round(bid_ask_spread.TimeToExpiration*365);

% Strike values for the maturity
ks = ks_cellArray{discountData_df.T_days == Tn_days};

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
plot(ks, SSVI_vols_cellArray{discountData_df.T_days == Tn_days}, "-g", "LineWidth",1);
hold on
plot(bid_ask_spread.logStrike(bid_ask_spread.T_days == Tn_days), ...
    bid_ask_spread.sigma_target(bid_ask_spread.T_days== Tn_days), ...
    "-c", "LineWidth",1);
hold off
xlabel("Log strike")
ylabel("BS implied vol")
legend(["Call bid", "Call ask", "Put bid", "Put ask", "SSVI", "Target"])
title("Implied volatilities for maturity " +Tn_days+ " days, "+ regexprep(dataset,'_',' '))
%% Implied variance plot
% For log-strikes in dataset
figure(3)
for i = 1:length(T_maturities)
    T_i = T_maturities(i);
    ks = ks_cellArray{discountData_df.T == T_i};
    implied_var = (SSVI_vols_cellArray{discountData_df.T == T_i}).^2*T_i;
    plot(ks, implied_var, "-", "LineWidth",1);
    hold on
end
hold off
legend(string(T_maturities))
title("Total implied variance plot, "+regexprep(dataset,'_',' '))
xlabel("Log strike")
ylabel("Total implied variance")
%% Check subset of the smaller maturities
T_small = T_maturities(T_maturities < 0.1);
figure(4)
for i = 1:length(T_small)
    T_i = T_small(i);
    ks = ks_cellArray{discountData_df.T == T_i};
    implied_var = (SSVI_vols_cellArray{discountData_df.T == T_i}).^2*T_i;
    plot(ks, implied_var, "-", "LineWidth",1);
    hold on
end
hold off
legend(string(T_small))
title("Subset: Total implied variance plot, "+ regexprep(dataset,'_',' '))
xlabel("Log strike")
ylabel("Total implied variance")
%% Implied variance plot
% For a fixed set of log-strikes (interpolation and extrapolation)
k_set = -0.4:0.01:0.4; % set of log-strikes to compute over
implied_var_est = zeros(length(T_maturities), length(k_set));
SSVI_vol_est = zeros(length(T_maturities), length(k_set));
figure(5)
for i = 1:length(T_maturities)
    T_i = T_maturities(i);
    thetaT = discountData_df.TotImplVar(discountData_df.T == T_i);
    SSVI_vol = SSVIimpliedVolatility(thetaT, T_i, k_set, ...
            calibration_params.rho, calibration_params.eps, ...
            calibration_params.gamma1, calibration_params.gamma2, ...
            calibration_params.beta1, calibration_params.beta2);
    implied_var = SSVI_vol.^2*T_i;
    implied_var_est(i,:) = implied_var;
    SSVI_vol_est(i,:) = SSVI_vol;
    plot(k_set, implied_var, "-", "LineWidth",1);
    hold on
end
hold off
legend(string(T_maturities))
title("Total implied variance plot (inter-&extrapolated), "+regexprep(dataset,'_',' '))
xlabel("Log strike")
ylabel("Total implied variance")

%% 
check_crossing(implied_var_est)
% check_crossing([[1,2,3]; [2,4,6]; [5,3,8]]) %test
%% Volatility surface
[T_axis,k_axis] = meshgrid(T_maturities,k_set);
surf(T_axis,k_axis,SSVI_vol_est')
%surf(T_axis,k_axis,implied_var_est')
xlabel("Time to maturity (T)")
ylabel("log-strike (k)")
zlabel("SSVI Volatility")
%% 
% Check for a crossing on the volatility surface (calendar spread arbitrage)
% implied_var_matrix: row_i -> set of implied variances for maturity Ti
% Note: this is just a simple test (only checks adjacent lines)
function check_crossing(implied_var_matrix)
    % Flag to check if values increase from row to row (line to line)
    values_increase = true;
    
    % Iterate through the rows (lines)
    for i = 1:size(implied_var_matrix, 1) - 1
        current_row = implied_var_matrix(i, :);
        next_row = implied_var_matrix(i + 1, :);
        
        % Check if values increase from row to row
        if any(current_row > next_row)
            values_increase = false;
            break; 
        end
    end
    
    % Display the result
    if values_increase
        disp('No crossings (no calendar spread arbitrage)');
    else
        disp('Lines cross (calendar spread arbitrage exists)');
    end

end




