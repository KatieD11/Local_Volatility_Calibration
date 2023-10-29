% Find the total implied variance using SVI framework
% Use interpolation to estimate thetaT for maturities not in the data set
% w = SSVItotalImpliedVariance(thetaT, T, k, rho, eps, gamma1, gamma2, beta1, beta2)
function w = SSVItotalImpliedVariance(discountData_df, T, k, rho, eps, ...
    gamma1, gamma2, beta1, beta2)
    % Define eta
    eta = 2.016048*exp(eps);

    % Find thetaT
    thetaT = interp1(discountData_df.T,discountData_df.TotImplVar,T, 'linear', 'extrap');
    if (thetaT <=0)
        % Extrapolate with the condition that thetaT at t=0 is zero
        T0 = [0; discountData_df.T]; TotImplVar0 = [0; discountData_df.TotImplVar];
        thetaT = interp1(T0,TotImplVar0,T, 'linear', 'extrap');
    end
    % Note: extrapolation from T0 works better for options with T <= 31
%     T0 = [0; discountData_df.T]; TotImplVar0 = [0; discountData_df.TotImplVar];
%     thetaT = interp1(T0,TotImplVar0,T, 'linear', 'extrap');

    % Define phi and w (for given thetaT, T, k)
    phi = @(eta, theta) eta./(theta.^gamma1.*(1+beta1.*theta).^gamma2.* ...
        (1+beta2.*theta).^(1-gamma1-gamma2));
    w = thetaT/2.*(1+rho.*phi(eta,thetaT).*k + ...
        sqrt((phi(eta,thetaT).*k + rho).^2 + (1-rho.^2)));

end