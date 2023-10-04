% Find bid-ask implied vol spread (to use for calibration)
clear; clc;
% Filtered dataset with BS implied volatilities
spx_df=readtable("Data/spx_20220401_filtered_optionDataWithImplVol.csv");
% Data set with discount factors (BT, QT) and ATM total implied variance
discountData_df=readtable("Data/spx_20220401_discountData.csv");

S0 = 4545.86;
total_open_interest_call = sum(spx_df.CallOpenInterest);
total_open_interest_put = sum(spx_df.PutOpenInterest);

% spx_df.Properties.VariableNames
% discountData_df.Properties.VariableNames

% Find sigma_ask and sigma_bid for the option dataset
% σask (Ti , Kj ) is the smallest Black-Scholes implied volatility for the call and put ask prices at maturity Ti and strike Kj
% σbid(Ti,Kj) is the largest Black-Scholes implied volatility for the call and put bid prices at maturity Ti and strike Kj
col_names = {'TimeToExpiration', 'Strike', 'logStrike', 'sigma_ask', ...
    'sigma_bid', 'open_interest_ask', 'open_interest_bid'};
col_types = {'double','double','double','double','double','double','double'}; 

% Create an empty table with specified column names and data types
option_df = table('Size', [0, length(col_names)], 'VariableNames', ...
    col_names, 'VariableTypes', col_types);

% [Start with nested loops, then try to make more efficient]
for i = 1:length(discountData_df.T)
    Ti = discountData_df.T(i);
    Ks = spx_df.Strike(spx_df.TimeToExpiration == Ti);
    %for Kj = spx_df.Strike(spx_df.TimeToExpiration == Ti)
    for j = 1:length(Ks)
        % Skip sets that are already recorded
        if (j>1 && Ks(j)==Ks(j-1))
            continue
        end
        Kj = Ks(j);
        % Filter by Ti and Kj
        filter = spx_df.TimeToExpiration == Ti & ...
            spx_df.Strike == Kj;
        % Find sigma_ask from call and put ask vols
        callAsk_BSvol = spx_df.callAsk_BSvol(filter);
        putAsk_BSvol = spx_df.putAsk_BSvol(filter);
        sigma_ask = min(callAsk_BSvol, putAsk_BSvol);
        % Find sigma_bid from call and put bid vols
        callBid_BSvol = spx_df.callBid_BSvol(filter);
        putBid_BSvol = spx_df.putBid_BSvol(filter);
        sigma_bid = max(callBid_BSvol, putBid_BSvol);

        % Store log strike
        kj = spx_df.logStrike(filter);

        if length(sigma_ask)>1
            sigma_ask=min(sigma_ask);
            sigma_bid=max(sigma_bid);
            kj = kj(1);
        end
        % Add open interest to table
        % Ask:
        open_int_call_ask = spx_df.CallOpenInterest(spx_df.callAsk_BSvol == sigma_ask);
        open_int_put_ask = spx_df.PutOpenInterest(spx_df.putAsk_BSvol == sigma_ask);
        if ~isempty(open_int_call_ask)
            open_int_ask = open_int_call_ask; 
        elseif ~isempty(open_int_put_ask)
            open_int_ask = open_int_put_ask; 
        else
            disp("Issue: no open interest for ask")
            open_int_ask=0;
        end
        % Bid:
        open_int_call_bid = spx_df.CallOpenInterest(spx_df.callBid_BSvol == sigma_bid);
        open_int_put_bid = spx_df.PutOpenInterest(spx_df.putBid_BSvol == sigma_bid);
        if ~isempty(open_int_call_bid)
            open_int_bid = open_int_call_bid; 
        elseif ~isempty(open_int_put_bid)
            open_int_bid = open_int_put_bid; 
        else
            disp("Issue: no open interest for bid")
            open_int_bid=0;
        end
        % Check for more than one value in open_int_ask or open_int_bid
        if length(open_int_ask) > 1
            open_int_ask = sum(open_int_ask);
        end
        if length(open_int_bid) > 1
            open_int_bid = sum(open_int_bid);
        end

        % Add the values to the table
        new_row = {Ti, Kj, kj, sigma_ask,sigma_bid, ...
            open_int_ask, open_int_bid};
        option_df = [option_df; new_row];
    end
end

option_df.open_interest_weight = (option_df.open_interest_ask + option_df.open_interest_bid)/ ...
    sum(option_df.open_interest_ask + option_df.open_interest_bid);
%% Save bid-ask spread data 
% Find target vol, computed as the mid point between the smallest bid-ask spread
option_df.sigma_target = 0.5*(option_df.sigma_ask + option_df.sigma_bid);
writetable(option_df, "Data/spx_20220401_option_bid_ask_spread.csv")
