% Find implied vols using SVI framework
% impliedVol = SSVIimpliedVolatility(thetaT, T, k, rho, eps)
function impliedVol = SSVIimpliedVolatility(thetaT, T, k, rho, eps)
    % Define constants (S3)
    gamma1 = 0.238; gamma2 = 0.253; 
    beta1 = exp(5.18); beta2 = exp(-3);
    eta = 2.016048*exp(eps);
    
    % Define phi and w (for given thetaT, T, k)
    phi = @(eta, theta) eta/(theta^gamma1*(1+beta1*theta)^gamma2* ...
        (1+beta2*theta)^(1-gamma1-gamma2));
    w = thetaT/2*(1+rho*phi(eta,thetaT)*k + ...
        sqrt((phi(eta,thetaT)*k + rho)^2 + (1-rho^2)));

    % Implied volatility
    impliedVol = sqrt(w/T);
end