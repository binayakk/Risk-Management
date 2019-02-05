data_file = 'Prices.csv'; % Prices file from Data folder
data      = csvread(file_path,1,0);

% Reading data without date column from data
prices    = data(:,2:5);

% Calculating log return
log_ret   = diff(log(prices));

ewma_mean = zeros(1,4);
ewma_covr = zeros(4,4);

% Usinglambda and theta for Exponential Weighted Moving Average (EWMA)
% to update mean and covariance matrix for the given period.
lambda    = 0.94;
theta     = 0.97;

% Using gamma for Spectral Risk Measure
gamma     = 30;


% Using EWMA to update mean and covariance matrix for the given period.
N = length(log_ret);
for i = 1:N
    ewma_covr = theta*ewma_covr + (1-theta)*transpose(log_ret(N+1-i,:)-ewma_mean)*(log_ret(N+1-i,:)-ewma_mean);
    
    ewma_mean = lambda*ewma_mean + (1-lambda)*log_ret(N+1-i,:);   
end

% Given market caps of 4 companies, value of portfolio and Alpha for VaR, ES

market_comp = [196.94 125.86 131.57 282.87];
wgt_per_stk = market_comp/sum(market_comp);
val_pf      = 1000000;
num_simul   = 50000;
alpha       = 0.99;

% Calculating One day VaR
loss_1_day      = -val_pf * (wgt_per_stk * transpose(exp(mvnrnd(ewma_mean, ewma_covr, num_simul))-1));
sort_loss_1_day = transpose (sort(loss_1_day, 'ascend'));
VaR_1_day       = sort_loss_1_day(ceil(num_simul * alpha),1);

% Calculating One day Expected Shortfall
exp_srt_f_1_day = 0;

for k = (num_simul*alpha+1) : num_simul
    exp_srt_f_1_day = exp_srt_f_1_day + sort_loss_1_day(ceil(k),1);
end

exp_srt_f_1_day = 1/(num_simul-num_simul*alpha) * (VaR_1_day * (ceil(num_simul*alpha) - num_simul*alpha) + exp_srt_f_1_day);
    
% Calculating One day Spectral Risk Measure - exponential weighting function
spec_1_day = 0;

for k = 1:num_simul
    spec_1_day = spec_1_day + sort_loss_1_day(ceil(k),1) * ( exp(-gamma) / (1 - exp(-gamma)) )*(exp(gamma^(k/num_simul))-exp(gamma^((k-1)/num_simul)));
end

%% Simulation to calculate 10 days risk meaasures (VaR, ES, Spectral)
num_days       = 10;
simul_ret      = zeros(10,4);
simul_loss     = zeros(10,1);
simul_day_loss = zeros(num_simul,1);

% Simulating 10 days losses
for i = 1:num_simul
    samp_mean = ewma_mean;
    covr = ewma_covr;
    for k = 1:num_days
        simul_ret(k,:) = mvnrnd(samp_mean,covr);
        samp_mean = samp_mean*lambda + (1-lambda)*simul_ret(k,:);
        covr = theta*covr + (1-theta)*(simul_ret(k,:)-samp_mean)*transpose(simul_ret(k,:)-samp_mean);
        simul_loss(k) = wgt_per_stk * transpose(exp(simul_ret(k,:)));
    end
    simul_day_loss(i) = -val_pf * (prod(simul_loss)-1); 
end

% Calculating 10 days VaR
sort_loss_K_day = sort(simul_day_loss,'ascend');
VaR_K_day = sort_loss_K_day(ceil(num_simul*alpha),1);

% Calculating 10 days ES
exp_srt_f_K_day=0;
for k =(num_simul*alpha+1) : num_simul
    exp_srt_f_K_day = exp_srt_f_K_day + sort_loss_K_day(ceil(k),1);
end
exp_srt_f_K_day = 1/(num_simul-num_simul*alpha) * (VaR_K_day*(ceil(num_simul*alpha) - num_simul*alpha) + exp_srt_f_K_day);

% Calculating 10 days Spectral Risk Measure
spec_K_day = 0;
for k = 1:num_simul
    spec_K_day = spec_K_day + sort_loss_K_day(ceil(k),1)* ( exp(-gamma) / (1 - exp(-gamma)) )*(exp(gamma^(k/num_simul))-exp(gamma^((k-1)/num_simul)));
end


%%
% Plotting Component Risk Percentages - Variance, VaR, ES for all 4 Stk
% returns
alpha22      = 0.97; 
num_trial    = 50; 
dollar_alloc = ones(1,4)*250000; 

ewma_mean22 = mean(log_ret(1:num_trial,:));
ewma_covr22 = cov(log_ret(1:num_trial,:));

VaR_c        = zeros(N-num_trial+1,4); 
exp_srt_f_c  = zeros(N-num_trial+1,4); 
variance_c   = zeros(N-num_trial+1,4); 
start_VaR       = norminv(alpha22); 
start_exp_srt_f = 1/(1-alpha22) * dollar_alloc * norminv(alpha22); 

for i = num_trial:N
    ewma_covr22 = theta*ewma_covr22 + (1-theta)*transpose(log_ret(i,:) - ewma_mean22)*(log_ret(i,:)-ewma_mean22);
    ewma_mean22 = lambda*ewma_mean22 + (1-lambda)*log_ret(i,:);
    sigma_theta = ewma_covr22*transpose(dollar_alloc);
    for j = 1:4
        VaR_c(i-49,j) = 100*(( -dollar_alloc(j)*ewma_mean22(j) ) + (dollar_alloc(j) *sigma_theta(j)*start_VaR) / sqrt( dollar_alloc*ewma_covr22*transpose(dollar_alloc) ))...
                                /(-dollar_alloc*transpose(ewma_mean22) + sqrt(dollar_alloc*ewma_covr22*transpose(dollar_alloc))*start_VaR);
                            
        exp_srt_f_c(i-49,j) = 100*(( -dollar_alloc(j)*ewma_mean22(j) ) + (dollar_alloc(j) *sigma_theta(j)*start_exp_srt_f(j)) / sqrt( dollar_alloc*ewma_covr22*transpose(dollar_alloc) ))...
                               /(-dollar_alloc*transpose(ewma_mean22) + sqrt(dollar_alloc*ewma_covr22*transpose(dollar_alloc))*start_exp_srt_f(j)); 
                            
        variance_c(i-49,j) = 100*(dollar_alloc(j)*sigma_theta(j) / (dollar_alloc*ewma_covr22*transpose(dollar_alloc)));
    end
end


date = data(50:N,1);
date = x2mdate(date);
subplot(3,1,1)
plot(date, variance_c);datetick('x', 12); xlabel('Date'); ylabel('Component Risk Percentage');title('Variance Component Risk Percentage of 4 Stocks Portfolio')
subplot(3,1,2)
plot(date, exp_srt_f_c);datetick('x', 12); xlabel('Date'); ylabel('Component Risk Percentage');title('Expected Shortfall Component Risk Percentage of 4 Stocks Portfolio')
subplot(3,1,3)
plot(date, VaR_c);datetick('x', 12); xlabel('Date'); ylabel('Component Risk Percentage');title('Value-at-Risk Component Risk Percentage of 4 Stocks Portfolio')
