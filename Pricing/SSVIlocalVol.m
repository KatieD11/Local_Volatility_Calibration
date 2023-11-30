% Find the local vol using SVI framework for the SPX-fit surface
% Use interpolation to estimate thetaT for maturities not in the data set
% localVol = SSVIlocalVol(discountData_df, T, k, rho, eps, gamma1, gamma2, beta1, beta2)
function localVol = SSVIlocalVol(discountData_df, T, k, rho, eps, ...
    gamma1, gamma2, beta1, beta2)
    % Define eta
    eta = 2.016048*exp(eps);

    % Find thetaT -- extrapolate from TotImplVar = 0 at T=0
    T0 = [0; discountData_df.T]; TotImplVar0 = [0; discountData_df.TotImplVar];
    thetaT = interp1(T0,TotImplVar0,T, 'linear', 'extrap');

    % Set up total implied variance w as a function of T and k
    w = @(k,T) SSVItotalImpliedVariance(discountData_df, T, k, ...
            rho, eps, gamma1, gamma2, beta1, beta2);

    % Define phi
    phi = @(eta, theta) eta./(theta.^gamma1.*(1+beta1.*theta).^gamma2.* ...
        (1+beta2.*theta).^(1-gamma1-gamma2));

    % Central finite difference estimate for the partial derivative wrt T
    dT=0.0001;
    dwdT = (w(k, T+dT) - w(k, T-dT))/(2*dT);

    % Compute partial derivatives wrt to k (analytical solutions)
    dwdk = thetaT/2.*(rho.*phi(eta,thetaT) + ...
        (phi(eta,thetaT).*(phi(eta,thetaT).*k + rho))./...
        sqrt((phi(eta,thetaT).*k + rho).^2 + (1-rho.^2)));
    d2wdk2 = thetaT/2.*(phi(eta,thetaT).^2.*(1-rho.^2) ./...
        ((phi(eta,thetaT).*k + rho).^2 + (1-rho.^2)).^(3/2));
    
    % Local volatility function (implemented as vol^2)
    localVar = dwdT./(1 - k./w(k,T) .* dwdk + ...
        1/4 * (-1/4 - 1./w(k,T) + k.^2./w(k,T).^2).*(dwdk).^2 ...
        + 1/2* d2wdk2);

    localVol =sqrt(localVar);

    % Conditions for absence of butterfly arbitrage:
    % i)  θϕ(θ) (1 + |ρ|)  <  4,  for  all  θ  >  0
    % ii) θϕ(θ)^2  (1 + |ρ|)  ≤  4,  for  all  θ  >  0
    if (thetaT.*phi(eta,thetaT)*(1+abs(rho)) >= 4)
        disp("Butterfly arbitrage: condition 1 violated")
    end
    if (thetaT.*phi(eta,thetaT).^2*(1+abs(rho)) > 4)
        disp("Butterfly arbitrage: condition 2 violated")
    end

end