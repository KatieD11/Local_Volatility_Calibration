% Finite diff. approx. of the local volatility function
% given the total implied variance w(k, T) parameterisation
% vol = LocalVolFD(dT, dk, w, k, T)
function vol = LocalVolFD(dT, dk, w, k, T)
    % Central finite difference estimates for the partial derivatives
    dwdT = @(k, T) (w(k, T+dT) - w(k, T-dT))/(2*dT);
    dwdk = @(k, T) (w(k+dk, T) - w(k-dk, T))/(2*dk);
    d2wdk2 = @(k, T) (w(k+2*dk, T) - 2*w(k, T) + w(k-2*dk, T))/(4*dk^2);
    
    % Local volatility function (implemented as vol^2)
    local_var = dwdT(k, T)./(1 - k./w(k,T) .* dwdk(k, T) + ...
        1/4 * (-1/4 - 1./w(k,T) + k.^2./w(k,T).^2).*(dwdk(k, T)).^2 ...
        + 1/2* d2wdk2(k,T));

    vol =sqrt(local_var);
end