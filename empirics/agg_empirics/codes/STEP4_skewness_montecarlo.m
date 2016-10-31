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
nseries = 2500;

% list of filters
filters = {'hp1600','hp1e5','bp6_32','bp2_80','fd'};
for i_filter = 1:length(filters)
	expr = strcat(filters{i_filter},'=',num2str(i_filter),';');
	eval(expr);
end
nfilters = length(filters);

% sample sizes
sample_sizes =  [124,224];
nsizes = length(sample_sizes);

% set of innovations
innovations = {'normal',...
	           'lognormal_low',...
			   'lognormal_high',...
			   'weibull_low',...
			   'weibull_high',...
			   'gamma_low',...
			   'gamma_high',...
			   };
for i_innov = 1:length(innovations)
	expr = strcat(innovations{i_innov},'=',num2str(i_innov),';');
	eval(expr);
end
ninnov = length(innovations);


%
sim_obs = cell(7,2); % for innovations, two sample sizes
sim_cyc = cell(7,2);

%% generate data for all innovations and sample lengths
for i_innov = 1:7
	for i_sample_length = 1:2
		
		sample_length = sample_sizes(i_sample_length);
		switch i_innov % formulas are from Wiki
			case normal
				cyc_u = randn(sample_length,nseries);
			case lognormal_low
				target_skewness = 0.95;
				myfun = @(s) target_skewness - (exp(s)+2)*sqrt(exp(s)-1);
				ssigmasq = fsolve(myfun,1);
				sample = lognrnd(0,sqrt(ssigmasq),sample_length,nseries);
				transformed_sample = zscore(sample); % force 0 mean, 1 std
				cyc_u = transformed_sample;
			case lognormal_high
				target_skewness = 6.2;
				myfun = @(s) target_skewness - (exp(s)+2)*sqrt(exp(s)-1);
				ssigmasq = fsolve(myfun,1);
				sample = lognrnd(0,sqrt(ssigmasq),sample_length,nseries);
				transformed_sample = zscore(sample); % force 0 mean, 1 std
				cyc_u = transformed_sample;
			case weibull_low
				% Note that skewness depend on k only!!
				target_skewness = 0.95;
				part1 = @(k) (gamma(1+2/k)-gamma(1+1/k)^2)^(-3/2);
				part2 = @(k) gamma(1+3/k)-3*gamma(1+1/k)*(gamma(1+2/k)-gamma(1+1/k)^2)...
					-gamma(1+1/k)^3;
				ttaufn = @(k) part1(k)*part2(k) - target_skewness;
				k = fsolve(ttaufn,1.5);
				llambda = sqrt( 1/( gamma(1+2/k)-gamma(1+1/k)^2 ) );
				sample = wblrnd(llambda,k,sample_length,nseries);
				cyc_u =  zscore(sample); % force 0mean 1variance
			case weibull_high
				% Note that skewness depend on k only!!
				target_skewness = 6.2;
				part1 = @(k) (gamma(1+2/k)-gamma(1+1/k)^2)^(-3/2);
				part2 = @(k) gamma(1+3/k)-3*gamma(1+1/k)*(gamma(1+2/k)-gamma(1+1/k)^2)...
					-gamma(1+1/k)^3;
				ttaufn = @(k) part1(k)*part2(k) - target_skewness;
				k = fsolve(ttaufn,1.5);
				llambda = sqrt( 1/( gamma(1+2/k)-gamma(1+1/k)^2 ) );
				sample = wblrnd(llambda,k,sample_length,nseries);
				cyc_u =  zscore(sample); % force 0mean 1variance
			case gamma_low
				target_skewness = 0.95;
				k = (2/target_skewness)^2;
				ttheta = sqrt(1/k);
				sample = gamrnd(k,ttheta,sample_length,nseries);
				cyc_u =  zscore(sample); % force 0mean 1variance
			case gamma_high
				target_skewness = 6.2;
				k = (2/target_skewness)^2;
				ttheta = sqrt(1/k);
				sample = gamrnd(k,ttheta,sample_length,nseries);
				cyc_u =  zscore(sample); % force 0mean 1variance				
		end
		
		trend_u = randn(sample_length,nseries);
		trend = zeros(sample_length,1);
		cyc = zeros(sample_length,1);
		obs_table = zeros(sample_length,nseries);
		cyc_table = obs_table;
		for i_series = 1:nseries
			trend(1) = 0 + 0.01 + 0.01*trend_u(1,i_series);
			cyc(1) = 1.5*0 - 0.58*0 + 0.01*cyc_u(1,i_series);
			cyc(2) = 1.5*cyc(1) - 0.58*0 + 0.01*cyc_u(2,i_series);
			for t = 2:length(trend)
				trend(t) = trend(t-1) + 0.01 + 0.01*trend_u(t,i_series);
				if t >= 3
					cyc(t) = 1.5*cyc(t-1) - 0.58*cyc(t-2) + 0.01*cyc_u(t,i_series);
				end
			end
			
			% store the observable and cyclical part
			obs_table(:,i_series) = trend+cyc;
			cyc_table(:,i_series) = cyc;
		end
		sim_obs{i_innov,i_sample_length} = obs_table;
		sim_cyc{i_innov,i_sample_length} = cyc_table;
	end
end

%% For each data series
rej_obs_table = cell(nfilters,ninnov,length(sample_sizes));
rej_cyc_table = cell(nfilters,ninnov,length(sample_sizes));

for i_innov = 1:ninnov
	for i_sample_length = 1:length(sample_sizes)
		raw_obs = sim_obs{i_innov,i_sample_length};
		raw_cyc = sim_cyc{i_innov,i_sample_length};
		for i_filter = 1:nfilters
			rej_obs_count = 0;
			rej_cyc_count = 0;
			for i_series = 1:nseries
				% apply filter by case
				switch i_filter
					case hp1600
						[~,detrended_obs] = hpfilter(raw_obs(:,i_series),1600);
						[~,detrended_cyc] = hpfilter(raw_cyc(:,i_series),1600);
					case hp1e5
						[~,detrended_obs] = hpfilter(raw_obs(:,i_series),1e5);
						[~,detrended_cyc] = hpfilter(raw_cyc(:,i_series),1e5);
					case bp6_32
						detrended_obs = bandpass(raw_obs(:,i_series),6,32);
						detrended_cyc = bandpass(raw_cyc(:,i_series),6,32);
					case bp2_80
						detrended_obs = bandpass(raw_obs(:,i_series),2,80);
						detrended_cyc = bandpass(raw_cyc(:,i_series),2,80);	
					case fd
						detrended_obs = diff(raw_obs(:,i_series));
						detrended_cyc = diff(raw_cyc(:,i_series));
				end
				
				% compute skewness and test
				[s_obs,t_obs] = skewTestAC(detrended_obs);
				[s_cyc,t_cyc] = skewTestAC(detrended_cyc);
				% skewCoeffs = s;
				if pnorm(abs(t_obs)) > 0.95
					rej_obs_count = rej_obs_count + 1;
				end
				if pnorm(abs(t_cyc)) > 0.95
					rej_cyc_count = rej_cyc_count + 1;
				end
				% # stopped here
			end
			rej_obs_table{i_filter,i_innov,i_sample_length}=rej_obs_count/nseries;
			rej_cyc_table{i_filter,i_innov,i_sample_length}=rej_cyc_count/nseries;
		end
	end
end
