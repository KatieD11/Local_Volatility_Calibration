% Given a dataset of option data info,
% calculate the corresponding discount factors (B(T), Q(T))
% optionData columns: TimeToExpiration, Strike, CallMktPrice, PutMktPrice 
% DFs = DiscountFactors(optionData, S0)
function DFs = DiscountFactors(optionData, S0)
    % Find the set of unique maturity times
    T_vals = unique(optionData.TimeToExpiration);
    
    % Put-Call parity
    % KjBTi − S0QTi = P(Ti,Kj)−C(Ti,Kj)
    BT = zeros(length(T_vals),1);
    QT = zeros(length(T_vals),1);
    for i=1:length(T_vals)
        C = optionData.CallMktPrice( ... 
            optionData.TimeToExpiration==T_vals(i));
        P = optionData.PutMktPrice( ... 
            optionData.TimeToExpiration==T_vals(i));
        % B -> [P(Ti,Kj)−C(Ti,Kj)] vector 
        B = P-C;
        % A -> [Kj, −S0] matrix
        K_vals = optionData.Strike( ... 
            optionData.TimeToExpiration==T_vals(i)); 
        A = [K_vals, - S0*ones(length(K_vals),1)];
    
        % Ax = B with x=[BTi; QTi]
        % Solve for x using Moore-Penrose pseudoinverse (pinv)
        Ainv = pinv(A);
        x = Ainv*B;
    
        % Store BTi and QTi
        BT(i) = x(1); QT(i) = x(2);
    end
    % Return discount factors 
    DFs = [BT, QT];
end