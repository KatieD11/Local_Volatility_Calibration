% Test inter- & extrapolation for parameters used in pricing
% r(T), q(T), F(T), B(T), thetaT

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

% Approximate qT (from QT) as average dividend yield over the period
discountData_df.qT = -log(discountData_df.QT)./discountData_df.T;
% Compute FT in the dataset
discountData_df.FT = S0*(discountData_df.QT)./(discountData_df.BT);

% Use interpolation to estimate r(T), q(T), FT, BT, thetaT
% for maturities not in the data set
r = @(T) interp1(discountData_df.T,discountData_df.rT,T, 'linear', 'extrap');
q = @(T) interp1(discountData_df.T,discountData_df.qT,T, 'linear', 'extrap');
F = @(T) interp1(discountData_df.T,discountData_df.FT,T, 'linear', 'extrap');
B = @(T) interp1(discountData_df.T,discountData_df.BT,T, 'linear', 'extrap'); %DF
%thetaT = @(T) interp1(discountData_df.T,discountData_df.TotImplVar,T, 'linear', 'extrap');
T0 = [0; discountData_df.T]; TotImplVar0 = [0; discountData_df.TotImplVar];
thetaT = @(T) interp1(T0,TotImplVar0,T, 'linear', 'extrap');

T_maturities = discountData_df.T;

%% Plot actual vs inter-&extra-polated values
T_extrap = 0: 0.02: (max(T_maturities)+0.05);

% short rate r(T)
figure(1)
plot(T_maturities, discountData_df.rT, "-b");
hold on 
plot(T_maturities, discountData_df.rT, ".b");
hold on 
plot(T_extrap, r(T_extrap), "or");
hold off
title("Rate r_T")
xlabel("Maturity T")
legend(["", "Actual", "Inter/extra-polated"])

% short rate q(T)
figure(2)
plot(T_maturities, discountData_df.qT, "-b");
hold on 
plot(T_maturities, discountData_df.qT, ".b");
hold on 
plot(T_extrap, q(T_extrap), "or");
hold off
title("Yield q_T")
xlabel("Maturity T")
legend(["", "Actual", "Inter/extra-polated"])

% Discount factor B(T)
figure(3)
plot(T_maturities, discountData_df.BT, "-b");
hold on 
plot(T_maturities, discountData_df.BT, ".b");
hold on 
plot(T_extrap, B(T_extrap), "or");
hold off
title("Discount factor B(T)")
xlabel("Maturity T")
legend(["", "Actual", "Inter/extra-polated"])

% ThetaT 
figure(4)
plot(T_maturities, discountData_df.TotImplVar, "-b");
hold on 
plot(T_maturities, discountData_df.TotImplVar, ".b");
hold on 
plot(T_extrap, thetaT(T_extrap), "or");
hold off
title("ThetaT")
xlabel("Maturity T")
legend(["", "Actual", "Inter/extra-polated"])

% F(T) 
figure(5)
plot(T_maturities, discountData_df.FT, "-b");
hold on 
plot(T_maturities, discountData_df.FT, ".b");
hold on 
plot(T_extrap, F(T_extrap), "or");
hold off
title("F(T)")
xlabel("Maturity T")
legend(["", "Actual", "Inter/extra-polated"])

%% Average rates / yields between maturities

Ts = discountData_df.T;
dT = Ts(1:end) - [0;Ts(1:end-1)];
r_int = -1./dT.*(log(discountData_df.BT(1:end)) - [0;log(discountData_df.BT(1:end-1))]);
q_int = -1./dT.*(log(discountData_df.QT(1:end)) - [0;log(discountData_df.QT(1:end-1))]);

figure(6)
stairs([0;Ts]*365, [r_int; r_int(end)], "-r");
hold on
plot(Ts*365, r_int, ".r");
plot(Ts*365, discountData_df.rT, "-b");
plot(Ts*365, discountData_df.rT, ".b");
hold off
title("Short rate estimates")
xlabel("t_i (days)")
ylabel("r")
legend(["Average rate over interval", "", "Average rate to maturity"], "Location","southeast")

figure(7)
stairs([0;Ts]*365, [q_int; q_int(end)], "-r");
hold on
plot(Ts*365, q_int, ".r");
plot(Ts*365, discountData_df.qT, "-b");
plot(Ts*365, discountData_df.qT, ".b");
hold off
title("Dividend yield estimates")
xlabel("t_i (days)")
ylabel("q")
legend(["Average yield over interval", "", "Average yield to maturity"], "Location","southeast")
