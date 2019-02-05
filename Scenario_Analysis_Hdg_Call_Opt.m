%% Scenario Analysis for a Hedged Call Option in the Black-Scholes Model

% Given parameters
r   = 0.0132; 
mu1 = 0.15475; 
vol = 0.2214;
S0  = 158.12; 
Stk = 170;
t   = 0;
T   = 0.25; 
del = 1/252; 

% Functions to calculate call price and delta
bs_call_pr = @(t,y,x) blsprice(y, Stk, r, T-t, x);
bs_del_hdg = @(t,y,x) blsdelta(y, Stk, r, T-t, x);

% Given scenarios and corresponding weights
y     = [-0.6 0.6 -0.4 0.4 -0.2 0.2]; % 5 days Log return 
wgt_y = [0.5 0.5 0.75 0.75 1 1]; % 5 days Log return weights

beta  = [0.5 0.75 1.25 1.5 1.75 2]; % 5 days Beta
x     = beta*vol; % 5 days volatility
wgt_x = [0.5; 0.75; 1; 1; 0.75; 0.5];
tau   = 5*del; % 5 days out of a year

scnro_num  = length(y);
scnro_num1 = length(x);

losses = zeros(6,6);

for i = 1:scnro_num
    for j = 1:scnro_num1
        s_t         = S0*(exp(y(i)));
        hdg         = bs_del_hdg(t,S0,vol);
        price_t     = bs_call_pr(t,S0,vol);
        price_5t    = bs_call_pr(t+tau,S0*exp(y(i)), x(j));
        losses(i,j) = -hdg*s_t + price_5t - price_t + hdg*S0;
        %losses(i,j) = -1*bs_del_hdg(t,S0,x(j))*S0*(exp(y(i))-1)+bs_call_pr(t+tau,S0*exp(y(i)),x(j))-bs_call_pr(t,S0,x(j));
    end
end

% worst case scenario with equal weights
worst_scenario_rm = max(losses(:));
[row,col]         = find(losses == worst_scenario_rm);
log_ret_vol_combo = [y(row) beta(col)];


% weighted worst case scenario
wgt_mat    = wgt_x*wgt_y;
wgt_losses = wgt_mat.*losses;

worst_scenario_wgt_rm = max(wgt_losses(:));
[wgt_row,wgt_col]     = find(wgt_losses == worst_scenario_wgt_rm);
log_ret_vol_wgt_combo = [y(wgt_row) beta(wgt_col)];
disp("The worst case scenario risk measure:");
disp(worst_scenario_rm);
disp("The log return/volatlity combination which achieves this measure:");
disp(log_ret_vol_combo);

disp("The weighted worst case scenario risk measure:");
disp(worst_scenario_wgt_rm);
disp("The log return/volatlity combination which achieves this measure:");
disp(log_ret_vol_wgt_combo);


