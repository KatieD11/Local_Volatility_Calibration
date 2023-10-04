% Find implied vols using SVI framework (with Heston-like phi fnc)
% impliedVol = SSVIimpliedVolatility_Heston(thetaT, T, k, rho, lambda)
function impliedVol = SSVIimpliedVolatility_Heston(thetaT, T, k, rho, lambda)
    % Define functions phi (Heston) and w
    phi = @(lambda, theta) 1/(lambda*theta)* ... 
        (1-(1-exp(-lambda*theta))/(lambda*theta));
    w = thetaT/2*(1+rho*phi(lambda, thetaT)*k + ...
        sqrt((phi(lambda,thetaT)*k + rho).^2 + (1-rho^2)));

    % Implied volatility
    impliedVol = sqrt(w/T);

    % Conditions for absence of butterfly arbitrage:
    % i)  θϕ(θ) (1 + |ρ|)  <  4,  for  all  θ  >  0
    % ii) θϕ(θ)^2  (1 + |ρ|)  ≤  4,  for  all  θ  >  0
    if (thetaT.*phi(lambda,thetaT)*(1+abs(rho)) >= 4)
        disp("Butterfly arbitrage: condition 1 violated")
    end
    if (thetaT.*phi(lambda,thetaT).^2*(1+abs(rho)) > 4)
        disp("Butterfly arbitrage: condition 2 violated")
    end

end