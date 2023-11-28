% BS put price using discount factors QT, BT
% p = BSput(T,K,S0,sigma,QT, BT)
function p = BSput(T,K,S0,sigma,QT, BT)
    d1 = (log(S0/K)+log(QT/BT)+0.5*sigma^2*T)/(sigma*sqrt(T));
    d2 = d1 - sigma*sqrt(T);
    p = K*BT*normcdf(-d2) - S0*QT*normcdf(-d1);
end