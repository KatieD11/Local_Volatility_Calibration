% Find the total implied variance using SVI framework for the SPX-fit surface
% Use interpolation to estimate thetaT for maturities not in the data set
% w = SSVItotalImpliedVariance(thetaT, T, k, rho, eps, gamma1, gamma2, beta1, beta2)
function w = SSVItotalImpliedVariance(discountData_df, T, k, rho, eps, ...
    gamma1, gamma2, beta1, beta2)
    % Define eta
    eta = 2.016048*exp(eps);

    % Extrapolate from TotImplVar = 0 at T=0
    T0 = [0; discountData_df.T]; TotImplVar0 = [0; discountData_df.TotImplVar];
    thetaT = interp1(T0,TotImplVar0,T, 'linear', 'extrap');

    % Define phi and w (for given thetaT, T, k)
    phi = @(eta, theta) eta./(theta.^gamma1.*(1+beta1.*theta).^gamma2.* ...
        (1+beta2.*theta).^(1-gamma1-gamma2));
    w = thetaT/2.*(1+rho.*phi(eta,thetaT).*k + ...
        sqrt((phi(eta,thetaT).*k + rho).^2 + (1-rho.^2)));

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