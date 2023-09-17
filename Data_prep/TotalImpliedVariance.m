% At-the-money total implied variance, for a given strike (Ki) and maturity (Ti)
% optionData columns: TimeToExpiration, Strike, CallMktPrice, PutMktPrice 
% S0 (Current price of the underlying asset), r (risk-free rate)
% [var, BSvol] = TotalImpliedVariance(Ki, Ti, r, S0, optionData)
function [var, BSvol] = TotalImpliedVariance(Ki, Ti, QT, BT, S0, optionData)
    % Find the strikes just below and just above the target strike
    strikes = unique(optionData.Strike(optionData.TimeToExpiration == Ti)); % Strikes from option data
    strike_below = max(strikes(strikes <= Ki));
    strike_above = min(strikes(strikes >= Ki));
    
    % Calculate interpolation weights
    weight_below = (strike_above - Ki) / (strike_above - strike_below);
    weight_above = 1 - weight_below;
    
    % Get the call and put prices corresponding to strike_above & strike_below
    filter_below = (optionData.TimeToExpiration == Ti) & ...
        (optionData.Strike == strike_below);
    filter_above = (optionData.TimeToExpiration == Ti) & ...
        (optionData.Strike == strike_above);
    call_below = optionData.CallMktPrice(filter_below);
    call_above = optionData.CallMktPrice(filter_above);
    put_below = optionData.PutMktPrice(filter_below);
    put_above = optionData.PutMktPrice(filter_above);
    % If more than one call / put price was returned (ie. 2 prices for same T and K)
    % average the prices
    if (length(call_below) > 1)
        call_below = mean(call_below);
        call_above = mean(call_above);
        put_below = mean(put_below);
        put_above = mean(put_above);
    end 

    % Interpolate call and put prices for the target strike and maturity
    call_i = weight_below * call_below + weight_above * call_above;
    put_i = weight_below * put_below + weight_above * put_above;
    
    % Find Black-Scholes implied volatilities (using interpolated prices)
    call_BSvol = fzero(@(BSvol) BScall(Ti,Ki,S0,BSvol,QT, BT) - call_i,0.1);
    put_BSvol = fzero(@(BSvol) BSput(Ti,Ki,S0,BSvol,QT, BT) - put_i,0.1); 

    % Average of the put and call Black-Scholes implied volatilities
    BSvol = (call_BSvol + put_BSvol)/2;
    
    % Forward at-the-money Black-Scholes total implied variance for maturity T 
    % θTi =σ^2(Ti,FTi)* Ti
    var = BSvol^2*Ti;
end
% Notes:
% Since the value of FTi may not lie directly on a strike in the data,
% compute the Black-Scholes implied volatility σ2 (Ti,FT ) by (linearly)
% interpolating the option prices between the strikes either side of the ATM strike. 
% Take the average of the put and call Black-Scholes implied volatilities computed using the interpolated prices.

