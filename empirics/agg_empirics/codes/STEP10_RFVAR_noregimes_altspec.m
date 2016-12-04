%% housekeeping
clear;
clc;
close all
addpath(genpath('~/Dropbox/matlabtools/BVARtoolbox'))

load('../data/manual_select_dummy_FRED')
y_table = [diff(ln_rSales) diff(ln_rGDP) diff(ln_rISratio)];
regimes = ones(size(y_table,1),1);

%% specify the empirical model
p = 3; N = size(y_table,2);
M = max(regimes);
model.N = N;
model.p = p;
model.K = model.N*model.p+1;
model.T = length(regimes);

%% A quick OLS on RFVAR for initial conditions
xxx = lagmatrix(y_table,1:p);
xxx = [ones(size(y_table,1),1) xxx];
X_OLS = xxx(p+1:end,:);
Y_OLS = y_table(p+1:end,:);
est_OLS = RFVAR_OLS(Y_OLS,X_OLS,N);
options.mmu_init = repmat(est_OLS.coeff(1:N),M,1);
for m = 1:M
	options.Ssigma_array_init(:,:,m) = est_OLS.Ssigma;
end
options.pphi_init = est_OLS.coeff(N+1:end);

%% specify priors
priors.mmu_mean = zeros(N*M,1);
priors.mmu_cov = 1e7*eye(N*M);
priors.pphi_mean = zeros(N*N*p,1);
priors.pphi_cov = 1e7*eye(N*N*p);
priors.Ssigma_SSR = 0;
priors.Ssigma_nnu = 0;

%% specify options for the Gibbs sampler
options.burnin = 2e4;
options.R = 5e4;
tic
draws = dummyVAR_Gibbs(y_table,regimes,model,priors,options);
toc

save dummyVAR_Gibbs_simpleregimes_altspec
%% look at result
median(draws.Ssigma_array,4)
mean(draws.pphi,3)

sig_level = 10;
mmu_table = [mean(draws.mmu,3),median(draws.mmu,3),prctile(draws.mmu,sig_level/2,3),prctile(draws.mmu,100-sig_level/2,3)];
Ssigma_table = [mean(draws.Ssigma_array,4),median(draws.Ssigma_array,4),prctile(draws.Ssigma_array,sig_level/2,4),prctile(draws.Ssigma_array,100-sig_level/2,4)];

%% Some OIRF
sales = 1; gdp = 2; cipi = 3;
irf_horizon = 10;
OIRF_regime1_salesshock = OIRF_RFVAR(sales,irf_horizon,mean(draws.pphi,3),mean(draws.Ssigma_array(:,:,1),3),model.N,model.p);
OIRF_regime2_salesshock = OIRF_RFVAR(sales,irf_horizon,mean(draws.pphi,3),mean(draws.Ssigma_array(:,:,2),3),model.N,model.p);
OIRF_regime1_gdpshock = OIRF_RFVAR(gdp,irf_horizon,mean(draws.pphi,3),mean(draws.Ssigma_array(:,:,1),3),model.N,model.p);
OIRF_regime2_gdpshock = OIRF_RFVAR(gdp,irf_horizon,mean(draws.pphi,3),mean(draws.Ssigma_array(:,:,2),3),model.N,model.p);
OIRF_regime1_cipishock = OIRF_RFVAR(cipi,irf_horizon,mean(draws.pphi,3),mean(draws.Ssigma_array(:,:,1),3),model.N,model.p);
OIRF_regime2_cipishock = OIRF_RFVAR(cipi,irf_horizon,mean(draws.pphi,3),mean(draws.Ssigma_array(:,:,2),3),model.N,model.p);

%% Plotting
% Defaults for this section
width = 4;     % Width in inches
height = 2;    % Height in inches
alw = 1.2;    % AxesLineWidth
fsz = 13;      % Fontsize
lw = 2.0;      % LineWidth
msz = 8;       % MarkerSize

% sales shock
figure('Units','inches',...
	'PaperPositionMode','auto');
subplot(3,1,1)
plot(0:irf_horizon-1,OIRF_regime1_salesshock(sales,:)','b-', ...
	0:irf_horizon-1,OIRF_regime2_salesshock(sales,:)','r--', ...
	0:irf_horizon-1,zeros(1,irf_horizon),'k:', ...
	'LineWidth',lw,'MarkerSize',msz);
title('Sales Response')
ylabel('Percent')
xlabel('Periods from Shock')
legend('Recession','Expansion')
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

subplot(3,1,2)
plot(0:irf_horizon-1,OIRF_regime1_salesshock(gdp,:)','b-', ...
	0:irf_horizon-1,OIRF_regime2_salesshock(gdp,:)','r--', ...
	0:irf_horizon-1,zeros(1,irf_horizon),'k:', ...
	'LineWidth',lw,'MarkerSize',msz);
title('GDP Response')
ylabel('Percent')
xlabel('Periods from Shock')
legend('Recession','Expansion')
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

subplot(3,1,3)
plot(0:irf_horizon-1,OIRF_regime1_salesshock(cipi,:)','b-', ...
	0:irf_horizon-1,OIRF_regime2_salesshock(cipi,:)','r--', ...
	0:irf_horizon-1,zeros(1,irf_horizon),'k:', ...
	'LineWidth',lw,'MarkerSize',msz);
title('Inventory Investment Response')
ylabel('Percent')
xlabel('Periods from Shock')
legend('Recession','Expansion')
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

print('../figures/salesshock','-depsc2','-r300');
if ispc % Use Windows ghostscript call
  system('gswin64c -o -q -sDEVICE=png256 -dEPSCrop -r300 -oimprovedExample_eps.png improvedExample.eps');
else % Use Unix/OSX ghostscript call
  system('gs -o -q -sDEVICE=png256 -dEPSCrop -r300 -oimprovedExample_eps.png improvedExample.eps');
end


%% plot gdp shock
figure
subplot(3,1,1)
plot(1:irf_horizon,OIRF_regime1_gdpshock(sales,:)',1:irf_horizon,OIRF_regime2_salesshock(sales,:)')
title('Sales Response')
legend('Recession','Expansion')

subplot(3,1,2)
plot(1:irf_horizon,OIRF_regime1_gdpshock(gdp,:)',1:irf_horizon,OIRF_regime2_salesshock(gdp,:)')
title('GDP Response')
legend('Recession','Expansion')

subplot(3,1,3)
plot(1:irf_horizon,OIRF_regime1_gdpshock(cipi,:)',1:irf_horizon,OIRF_regime2_salesshock(cipi,:)')
title('Inventory Investment Response')
legend('Recession','Expansion')

% CIPI shock
figure
subplot(3,1,1)
plot(1:irf_horizon,OIRF_regime1_cipishock(sales,:)',1:irf_horizon,OIRF_regime2_salesshock(sales,:)')
title('Sales Response')
legend('Recession','Expansion')

subplot(3,1,2)
plot(1:irf_horizon,OIRF_regime1_cipishock(gdp,:)',1:irf_horizon,OIRF_regime2_salesshock(gdp,:)')
title('GDP Response')
legend('Recession','Expansion')

subplot(3,1,3)
plot(1:irf_horizon,OIRF_regime1_cipishock(cipi,:)',1:irf_horizon,OIRF_regime2_salesshock(cipi,:)')
title('Inventory Investment Response')
legend('Recession','Expansion')
