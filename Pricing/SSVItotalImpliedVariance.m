% Find the total implied variance using SVI framework
% Use interpolation to estimate thetaT for maturities not in the data set
% w = SSVItotalImpliedVariance(thetaT, T, k, rho, eps, gamma1, gamma2, beta1, beta2)
function w = SSVItotalImpliedVariance(discountData_df, T, k, rho, eps, ...
    gamma1, gamma2, beta1, beta2)
    % Define eta
    eta = 2.016048*exp(eps);

    % Find thetaT
%     if (T< min(discountData_df.T) || T > max(discountData_df.T))
%         thetaT = interp1(discountData_df.T,discountData_df.TotImplVar,T, 'linear', 'extrap');
%         thetaT = max(thetaT,0.01);
%     else
%         thetaT = interp1(discountData_df.T,discountData_df.TotImplVar,T);
%     end

%     thetaT = interp1(discountData_df.T,discountData_df.TotImplVar,T);
%     if isnan(thetaT)
%         thetaT = min(discountData_df.TotImplVar);
%     end

    if (T< min(discountData_df.T))
        thetaT = min(discountData_df.TotImplVar);
    else
        thetaT = interp1(discountData_df.T,discountData_df.TotImplVar,T);
    end

    % Define phi and w (for given thetaT, T, k)
    phi = @(eta, theta) eta./(theta.^gamma1.*(1+beta1.*theta).^gamma2.* ...
        (1+beta2.*theta).^(1-gamma1-gamma2));
    w = thetaT/2.*(1+rho.*phi(eta,thetaT).*k + ...
        sqrt((phi(eta,thetaT).*k + rho).^2 + (1-rho.^2)));

end