% Find implied vols using SVI framework (with Heston-like phi fnc)
% impliedVol = SSVIimpliedVolatility_PowerLaw(thetaT, T, k, rho, eta, gamma)
function impliedVol = SSVIimpliedVolatility_PowerLaw(thetaT, T, k, ...
    rho, eta, gamma)
    % Define functions phi (Heston) and w
    phi = @(theta) eta/(theta^gamma*(1+theta)^(1-gamma));
    w = thetaT/2*(1+rho*phi(thetaT)*k + ...
        sqrt((phi(thetaT)*k + rho).^2 + (1-rho^2)));

    % Implied volatility
    impliedVol = sqrt(w/T);

    % Conditions for absence of butterfly arbitrage:
    % i)  θϕ(θ) (1 + |ρ|)  <  4,  for  all  θ  >  0
    % ii) θϕ(θ)^2  (1 + |ρ|)  ≤  4,  for  all  θ  >  0
    if (thetaT.*phi(thetaT)*(1+abs(rho)) >= 4)
        disp("Butterfly arbitrage: condition 1 violated")
    end
    if (thetaT.*phi(thetaT).^2*(1+abs(rho)) > 4)
        disp("Butterfly arbitrage: condition 2 violated")
    end

end