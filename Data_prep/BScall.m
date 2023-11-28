% BS call price using discount factors QT, BT
% c = BScall(T,K,S0,sigma,QT, BT)
function c = BScall(T,K,S0,sigma,QT, BT)
    d1 = (log(S0/K)+log(QT/BT)+0.5*sigma^2*T)/(sigma*sqrt(T));
    d2 = d1 - sigma*sqrt(T);
    c = S0*QT*normcdf(d1) - K*BT*normcdf(d2);
end