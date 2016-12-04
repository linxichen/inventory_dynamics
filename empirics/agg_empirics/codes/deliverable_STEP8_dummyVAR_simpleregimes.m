%% housekeeping
clear;
clc;
close all
addpath(genpath('~/Dropbox/matlabtools/BVARtoolbox'))

load('../data/manual_select_dummy_FRED')
y_table = [diff(ln_rSales) diff(ln_rGDP) share_rCIPI_potential(2:end,:)];
regimes = simple_regimes(2:end,:);

%% specify the empirical model
p = 2; N = size(y_table,2);
M = max(regimes);
model.N = N;
model.p = p;
model.K = model.N*model.p+1;
model.T = length(regimes);

%% specify priors
priors.mmu_mean = zeros(N*M,1);
priors.mmu_cov = 1e7*eye(N*M);
priors.pphi_mean = zeros(N*N*p,1);
priors.pphi_cov = 1e7*eye(N*N*p);
priors.Ssigma_SSR = 0;
priors.Ssigma_nnu = 0;

%% A quick OLS on RFVAR for initial conditions
xxx = lagmatrix(y_table,1:p);
xxx = [ones(size(y_table,1),1) xxx];
X_OLS = xxx(p+1:end,:);
Y_OLS = y_table(p+1:end,:);
est_OLS = RFVAR_OLS(Y_OLS,X_OLS,N);
Ssigma_OLS = est_OLS.Ssigma;
options.mmu_init = repmat(est_OLS.coeff(1:N),M,1);
for m = 1:M
	options.Ssigma_array_init(:,:,m) = est_OLS.Ssigma;
end
pphi_OLS = est_OLS.coeff(N+1:end);
options.pphi_init = pphi_OLS;

%% specify options for the Gibbs sampler
options.burnin = 2e4;
options.R = 5e4;
tic
draws = dummyVAR_Gibbs(y_table,regimes,model,priors,options);
toc

save dummyVAR_Gibbs_simpleregimes
%% look at result
median(draws.Ssigma_array,4)
mean(draws.pphi,3)
mean(draws.mmu,3)


sig_level = 10;
mmu_table = [mean(draws.mmu,3),median(draws.mmu,3),prctile(draws.mmu,sig_level/2,3),prctile(draws.mmu,100-sig_level/2,3)];
Ssigma_table = [mean(draws.Ssigma_array,4),median(draws.Ssigma_array,4),prctile(draws.Ssigma_array,sig_level/2,4),prctile(draws.Ssigma_array,100-sig_level/2,4)];

%% Some OIRF
% specification
sales = 1; gdp = 2; cipi = 3;
irf_horizon = 10;

% intialize
OIRF_regime1_salesshock = zeros(model.N,irf_horizon,options.R);
OIRF_regime2_salesshock = zeros(model.N,irf_horizon,options.R);
OIRF_regime1_gdpshock = zeros(model.N,irf_horizon,options.R);
OIRF_regime2_gdpshock = zeros(model.N,irf_horizon,options.R);
OIRF_regime1_cipishock = zeros(model.N,irf_horizon,options.R);
OIRF_regime2_cipishock = zeros(model.N,irf_horizon,options.R);


for i_draw = 1:options.R
	pphi = draws.pphi(:,i_draw);
	Ssigma_array = draws.Ssigma_array(:,:,:,i_draw);
	OIRF_regime1_salesshock(:,:,i_draw) = OIRF_RFVAR(sales,irf_horizon,pphi,Ssigma_array(:,:,1),model.N,model.p);
	OIRF_regime2_salesshock(:,:,i_draw) = OIRF_RFVAR(sales,irf_horizon,pphi,Ssigma_array(:,:,2),model.N,model.p);
	OIRF_regime1_gdpshock(:,:,i_draw) = OIRF_RFVAR(gdp,irf_horizon,pphi,Ssigma_array(:,:,1),model.N,model.p);
	OIRF_regime2_gdpshock(:,:,i_draw) = OIRF_RFVAR(gdp,irf_horizon,pphi,Ssigma_array(:,:,2),model.N,model.p);
	OIRF_regime1_cipishock(:,:,i_draw) = OIRF_RFVAR(cipi,irf_horizon,pphi,Ssigma_array(:,:,1),model.N,model.p);
	OIRF_regime2_cipishock(:,:,i_draw) = OIRF_RFVAR(cipi,irf_horizon,pphi,Ssigma_array(:,:,2),model.N,model.p);
end
%% Plotting
% Defaults for this section
width = 4;     % Width in inches
height = 2;    % Height in inches
alw = 1.2;    % AxesLineWidth
fsz = 13;      % Fontsize
lw = 2.0;      % LineWidth
msz = 8;       % MarkerSize

%% plot for RFVAR_OLS
OIRF_OLS_salesshock = OIRF_RFVAR(sales,irf_horizon,pphi_OLS,Ssigma_OLS,model.N,model.p);
% sales shock
figure('Units','inches',...
	'PaperPositionMode','auto');
subplot(3,1,1)
plot(0:irf_horizon-1,OIRF_OLS_salesshock(sales,:)','b-', ...
	0:irf_horizon-1,OIRF_OLS_salesshock(sales,:)','r--', ...
	0:irf_horizon-1,zeros(1,irf_horizon),'k:', ...
	'LineWidth',lw,'MarkerSize',msz);
title('Sales Response')
ylabel('Percent')
xlabel('Periods from Shock')
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

subplot(3,1,2)
plot(0:irf_horizon-1,OIRF_OLS_salesshock(gdp,:)','b-', ...
	0:irf_horizon-1,OIRF_OLS_salesshock(gdp,:)','r--', ...
	0:irf_horizon-1,zeros(1,irf_horizon),'k:', ...
	'LineWidth',lw,'MarkerSize',msz);
title('GDP Response')
ylabel('Percent')
xlabel('Periods from Shock')
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

subplot(3,1,3)
plot(0:irf_horizon-1,OIRF_OLS_salesshock(cipi,:)','b-', ...
	0:irf_horizon-1,OIRF_OLS_salesshock(cipi,:)','r--', ...
	0:irf_horizon-1,zeros(1,irf_horizon),'k:', ...
	'LineWidth',lw,'MarkerSize',msz);
title('Inventory Investment Response')
ylabel('Percent')
xlabel('Periods from Shock')
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

%% sales shock comparison with unconditional VAR
figure('Units','inches',...
	'PaperPositionMode','auto');
set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
subplot(3,1,1)
plot(0:irf_horizon-1,median(OIRF_regime1_salesshock(sales,:,:),3)','b-', ...
	0:irf_horizon-1,median(OIRF_regime2_salesshock(sales,:,:),3)','r--', ...
	0:irf_horizon-1,OIRF_OLS_salesshock(sales,:)','k-.', ...
	0:irf_horizon-1,zeros(1,irf_horizon),'k:', ...
	'LineWidth',lw,'MarkerSize',msz);
title('Sales Response')
ylabel('Percent')
xlabel('Periods from Shock')
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

subplot(3,1,2)
plot(0:irf_horizon-1,median(OIRF_regime1_salesshock(gdp,:,:),3)','b-', ...
	0:irf_horizon-1,median(OIRF_regime2_salesshock(gdp,:,:),3)','r--', ...
	0:irf_horizon-1,OIRF_OLS_salesshock(gdp,:)','k-.', ...
	0:irf_horizon-1,zeros(1,irf_horizon),'k:', ...
	'LineWidth',lw,'MarkerSize',msz);
title('GDP Response')
ylabel('Percent')
xlabel('Periods from Shock')
legend('Expansion','Recession','Unconditional')
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

subplot(3,1,3)
plot(0:irf_horizon-1,median(OIRF_regime1_salesshock(cipi,:,:),3)','b-', ...
	0:irf_horizon-1,median(OIRF_regime2_salesshock(cipi,:,:),3)','r--', ...
	0:irf_horizon-1,OIRF_OLS_salesshock(cipi,:)','k-.', ...
	0:irf_horizon-1,zeros(1,irf_horizon),'k:', ...
	'LineWidth',lw,'MarkerSize',msz);
title('Inventory Investment Response')
ylabel('Percent')
xlabel('Periods from Shock')
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

print('../figures/salesshock_median','-depsc2','-r300');

%% sales shock, cipi response with CIPI across 2 regimes
figure('Units','inches',...
	'PaperPositionMode','auto');
set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
subplot(2,1,1)
plot(0:irf_horizon-1,median(OIRF_regime1_salesshock(cipi,:,:),3)','b-', ...
	0:irf_horizon-1,prctile(OIRF_regime1_salesshock(cipi,:,:),20,3)','r-.', ...
	0:irf_horizon-1,prctile(OIRF_regime1_salesshock(cipi,:,:),80,3)','r.-.', ...
	0:irf_horizon-1,zeros(1,irf_horizon),'k:', ...
	'LineWidth',lw,'MarkerSize',msz);
title('Inventory Investment Response in Expansions')
legend('median','20 Percentile','80 Percentile')
ylabel('Percent')
xlabel('Periods from Shock')
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

subplot(2,1,2)
plot(0:irf_horizon-1,median(OIRF_regime2_salesshock(cipi,:,:),3)','b-', ...
	0:irf_horizon-1,prctile(OIRF_regime2_salesshock(cipi,:,:),20,3)','r-.', ...
	0:irf_horizon-1,prctile(OIRF_regime2_salesshock(cipi,:,:),80,3)','r.-.', ...
	0:irf_horizon-1,zeros(1,irf_horizon),'k:', ...
	'LineWidth',lw,'MarkerSize',msz);
title('Inventory Investment Response in Recessions')
legend('median','20 Percentile','80 Percentile')
ylabel('Percent')
xlabel('Periods from Shock')
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

print('../figures/salesshock_withCI','-depsc2','-r300');

%% Get the RF covariance at all lags
Ssiga_array = median(draws.Ssigma_array,4);
sim_data_regime1 = simulate_RFVAR(1e5,pphi,Ssigma_array(:,:,1),model.N,model.p);
sim_data_regime2 = simulate_RFVAR(1e5,pphi,Ssigma_array(:,:,2),model.N,model.p);

[sales_CIPI_corr_regime,pval1] = corr(sim_data_regime1(sales,:)',sim_data_regime1(cipi,:)');
[sales_CIPI_corr_regime2,pval2] = corr(sim_data_regime2(sales,:)',sim_data_regime2(cipi,:)');

save results_step8