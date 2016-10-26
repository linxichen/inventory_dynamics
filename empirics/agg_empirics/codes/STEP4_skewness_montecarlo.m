% This script perform Monte Carlo exerises to see if Bai and Ng skewness
% test survives detrending well # Shan Table. Follows the Experiment set up
% in Psaradakis Sola 2003 JAE paper "detrending"
%% housekeeping
clearvars;
close all;
clc;
clc;
addpath(genpath('helpers'));

%% specify parameters of experiment
nsize = [124 224];
nseries = 2500;
small_size = 1;
large_size = 2;
normal = 1;
lognormal_low = 2; % 
lognormal_high = 3;
weibull_low = 4;
weibull_high = 5;
gamma_low = 6;
gamma_high = 7;
sim_obs = cell(7,2); % for innovations, two sample sizes
sim_cyc = cell(7,2);

%% log-normal play ground

%% generate data for one test
i_sample_length = 1;
i_innov = 1;

sample_length = nsize(i_sample_length);
switch i_innov % formulas are from Wiki
	case normal
		cyc_u = randn(sample_length,nseries);
	case lognormal_low
		target_skewness = 0.95;
		myfun = @(s) target_skewness - (exp(s)+2)*sqrt(exp(s)-1);
		ssigmasq = fsolve(myfun,1);
		sample = lognrnd(0,sqrt(ssigmasq),sample_length,nseries);
		transformed_sample = (sample-mean(sample))/(std(sample)); % force 0 mean, 1 std
		cyc_u = transformed_sample;
	case lognormal_high
		target_skewness = 6.2;
		myfun = @(s) target_skewness - (exp(s)+2)*sqrt(exp(s)-1);
		ssigmasq = fsolve(myfun,1);
		sample = lognrnd(0,sqrt(ssigmasq),sample_length,nseries);
		transformed_sample = (sample-mean(sample))/(std(sample)); % force 0 mean, 1 std
		cyc_u = transformed_sample;
	case weibull_low
		% #stopped here
		target_skewness = 0.95;
		k = 0.5;
		mmufn = @(l) l*gamma(1+1/k);
		ssigmafn = @(l) sqrt( l^2*( gamma(1+2/k)-(gamma(1+1/k))^2 ));
		ttaufn = @(l) (gamma(1+3/k)*l^3-3*mmufn(l)*(ssigmafn(l)^2)-mmufn(l)^3)...
			/ (ssigmafn(l)^3) - target_skewness;
		llambda = fsolve(ttaufn,10);
		cyc_u = wblrnd(llambda,k,1e7,1);
end

trend_u = randn(sample_length,nseries);
trend = zeros(sample_length,1);
cyc = zeros(sample_length,1);
obs_table = zeros(sample_length,nseries);
cyc_table = obs_table;
for i_series = 1:nseries
	trend(1) = 0 + 0.01 + 0.01*trend_u(1);
	cyc(1) = 1.5*0 - 0.58*0 + 0.01*cyc_u(1);
	cyc(2) = 1.5*cyc(1) - 0.58*0 + 0.01*cyc_u(2);
	for t = 2:length(trend)
		trend(t) = trend(t-1) + 0.01 + 0.01*trend_u(t);
		if t >= 3
			cyc(t) = 1.5*cyc(t-1) - 0.58*cyc(t-2) + 0.01*cyc_u(t);
		end
	end
	
	% store the observable and cyclical part
	obs_table(:,i_series) = trend+cyc;
	cyc_table(:,i_series) = cyc;
end
sim_obs{i_innov,i_sample_length} = obs_table;
sim_cyc{i_innov,i_sample_length} = cyc_table;


%% computes skewness for each variable
%  for ISratio level
format bank
disp('========== IS ratio =============');
raw = log(rISratio);
[~,HP1] = hpfilter(raw,1600);
[~,HP2] = hpfilter(raw,1e5);
BP1 = bandpass(raw,6,32);
BP2 = bandpass(raw,2,80);
FD = diff(raw);

disp('---------------------------')
disp('HP(1600)');
x = HP1;
[s,t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('HP(1e5)');
x = HP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(6,32)');
x = BP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(2,80)');
x = BP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%  for ISratio GROWTH
disp('---------------------------')
disp('FD');
x = FD;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%% for inventory stock
format bank
disp('========== Inventory Stock =============');
raw = log(rInventory);
[~,HP1] = hpfilter(raw,1600);
[~,HP2] = hpfilter(raw,1e5);
BP1 = bandpass(raw,6,32);
BP2 = bandpass(raw,2,80);
FD = diff(raw);

disp('---------------------------')
disp('HP(1600)');
x = HP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('HP(1e5)');
x = HP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(6,32)');
x = BP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(2,80)');
x = BP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%  for ISratio GROWTH
disp('---------------------------')
disp('FD');
x = FD;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%% for CIPI/GDP
format bank
disp('========== CIPI/GDP =============');
raw = share_CIPI;
[~,HP1] = hpfilter(raw,1600);
[~,HP2] = hpfilter(raw,1e5);
BP1 = bandpass(raw,6,32);
BP2 = bandpass(raw,2,80);
FD = diff(raw);

disp('---------------------------')
disp('HP(1600)');
x = HP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('HP(1e5)');
x = HP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(6,32)');
x = BP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(2,80)');
x = BP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%  for steepness
disp('---------------------------')
disp('FD');
x = FD;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%% for sales
disp('========== Sales =============');
raw = log(rSalesGoods);
[~,HP1] = hpfilter(raw,1600);
[~,HP2] = hpfilter(raw,1e5);
BP1 = bandpass(raw,6,32);
BP2 = bandpass(raw,2,80);
FD = diff(raw);

disp('---------------------------')
disp('HP(1600)');
x = HP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('HP(1e5)');
x = HP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(6,32)');
x = BP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(2,80)');
x = BP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%  for deepness
disp('---------------------------')
disp('FD');
x = FD;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

