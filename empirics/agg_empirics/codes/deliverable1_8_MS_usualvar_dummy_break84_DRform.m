% Example Script MS_Regress_Fit.m

clear;
close all;
clc;

diary('deliverable_1_7.txt')
diary on

addpath('./MS_Regress_FEX_1.09/m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');
addpath(genpath('~/dropbox/matlabtools/BVARtoolbox'));

% load data
load manual_select_dummy_FRED
%data = 100*diff(data(:,2:3));
start = 0; % first start date
shift = 0; % shift some of the variable forward or backward
% forward_resid_rCIPI = data(start+2+shift:end,2);
% share_rCIPI_pot = (rCIPI(2:end) - rCIPI(1:end-1))./GDPPOT(1:end-1);
share_rCIPI_pot = share_rCIPI_potential(2:end);
% forward_stock = data(start+1+shift:end,3);
% match_GDP = ln_rGDP;
match_GDP = (rGDP(2:end) - rGDP(1:end-1))./rGDP(1:end-1);
match_sales = (rSalesGoods(2:end) - rSalesGoods(1:end-1))./rSalesGoods(1:end-1);
data = 100*[match_sales match_GDP share_rCIPI_pot];

% data label
datelabel = (1949.25:0.25:2016.50)'; % because first differenced
yearnum = floor(datelabel);
monthnum = 12*(datelabel - yearnum)+2;
date_serial = datenum(yearnum,monthnum,ones(size(yearnum)));
break_dummy = (datelabel >= 1984.00);

% crated dependent vars and regressors
lags = 5;                            % decide on lags
Y = data;
nvar = size(Y,2);
YLAG = lagmatrix(Y,1:lags);
dep=Y(lags+1:end,:);                  % Defining dependent variable from .mat file
constVec=ones(length(dep),1);       % Defining a constant vector in mean equation (just an example of how to do it)
safe_dum = break_dummy(lags+1:end,:);
const_dum = constVec.*safe_dum;
YLAG = YLAG(lags+1:end,:); % discard first p lags
YLAG_dum = YLAG.*safe_dum;
indep{1}=[constVec const_dum YLAG ];                  % Defining some explanatory variables
indep{2}=[constVec const_dum YLAG ];                  % Defining some explanatory variables
indep{3}=[constVec const_dum YLAG ];                  % Defining some explanatory variables

k=2;                                % Number of States
S{1}=ones(1,2+lags*nvar+1);S{1}(1,end)=0;        % Defining which parts of the equation will switch states (column 1 and variance only)
S{2}=ones(1,2+lags*nvar+1);S{2}(1,end)=0;    % Defining which parts of the equation will switch states (column 1 and variance only)
S{3}=ones(1,2+lags*nvar+1);S{3}(1,end)=0;     % Defining which parts of the equation will switch states (column 1 and variance only)

advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details
advOpt.diagCovMat=1;
advOpt.useMEX=0;
advOpt.optimizer='fmincon';

[step1_Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model
save deliverable1_8.mat

%% plot regimes 2
safe_dates = date_serial(start+1+lags:end-shift,:);
figure
plot(safe_dates,step1_Spec_Out.smoothProb(:,2));
h = gca;
datetick('x','yyyy','keepticks')
xlabel('Time');
% ylabel('Smoothed States Probabilities');
legend('Recession Probability');
axis tight
recessband = recessionplot;

%% collect residual as z
rec_flag = step1_Spec_Out.smoothProb(:,2) > 0.7;
z = step1_Spec_Out.resid;
lags = 2;
zlag = lagmatrix(z,1:lags);
zlag_safe = zlag(lags+1:end,:);
z_safe = z(lags+1:end,:);
rec_flag_safe = rec_flag(lags+1:end,:);

% select those that are in recessions
z_rec = z_safe(rec_flag_safe,:);
zlag_rec = zlag_safe(rec_flag_safe,:);
temp = z_rec';
Y = temp(:);
X = kron(zlag_rec,eye(nvar));
pphi_rec = (X'*X)\(X'*Y);
resid_rec_vec = Y-X*pphi_rec;
resid_rec = reshape(resid_rec_vec,nvar,length(z_rec)); resid_rec = resid_rec';

% select those that are in expansions
z_exp = z_safe(~rec_flag_safe,:);
zlag_exp = zlag_safe(~rec_flag_safe,:);
temp = z_exp';
Y = temp(:);
X = kron(zlag_exp,eye(nvar));
pphi_exp = (X'*X)\(X'*Y); 
resid_exp_vec = Y-X*pphi_exp;
resid_exp = reshape(resid_exp_vec,nvar,length(z_exp)); resid_exp = resid_exp';

%% collect the second step residual as epsilon and regress again
E = zeros(length(rec_flag_safe),nvar);
E(rec_flag_safe,:) = resid_rec;
E(~rec_flag_safe,:) = resid_exp;

dep=E;                  % Defining dependent variable from .mat file
indep{1}=[ones(length(E),1)];                  % Defining some explanatory variables
indep{2}=[ones(length(E),1)];                  % Defining some explanatory variables
indep{3}=[ones(length(E),1)];                  % Defining some explanatory variables

k=2;                                % Number of States
S{1}=[0 1];        % Defining which parts of the equation will switch states (column 1 and variance only)
S{2}=[0 1];      % Defining which parts of the equation will switch states (column 1 and variance only)
S{3}=[0 1];       % Defining which parts of the equation will switch states (column 1 and variance only)

advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details
advOpt.diagCovMat=0;
advOpt.useMEX=1;
advOpt.optimizer='fminsearch';


[step2_Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model
% constraint the intercept to be zero
for i_var = 1:nvar
	step2_Spec_Out.advOpt.constCoeff.nS_Param{i_var} = {0};
end
[step2_Spec_Out]=MS_Regress_Fit(dep,indep,k,S,step2_Spec_Out.advOpt); % Estimating the model
save deliverable1_7.mat

%% plot regimes
% Plotting
% Defaults for this section
width = 4;     % Width in inches
height = 2;    % Height in inches
alw = 1.2;    % AxesLineWidth
fsz = 13;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

safe_dates = date_serial(start+1+lags:end-shift,:);
figure('Units','inches',...
	'PaperPositionMode','auto');
plot(safe_dates,step1_Spec_Out.smoothProb(lags+1:end,2),...
	'LineWidth',lw,'MarkerSize',msz);
% set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
h = gca;
datetick('x','yyyy','keepticks')
xlabel('Time');
ylabel('Probability');
title('Autoregressive Regime, Smoothed')
axis tight
recessband1 = recessionplot;
% legend('Probability');
print('../../figures/AR_regime','-depsc2','-r300');


safe_dates = date_serial(start+1+lags:end-shift,:);
figure
plot(safe_dates,step2_Spec_Out.smoothProb(:,2),...
	'LineWidth',lw,'MarkerSize',msz);
% set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
h = gca;
datetick('x','yyyy','keepticks')
xlabel('Time');
ylabel('Probability');
% legend('Probability');
title('High Volatility Regime, Smoothed')
axis tight
recessband2 = recessionplot;
print('../../figures/COV_regime','-depsc2','-r300');

%% conditonal variance
low_vola = step2_Spec_Out.Coeff.covMat{1};
high_vola = step2_Spec_Out.Coeff.covMat{2};
rec_high_sim = simulate_RFVAR(1e5,pphi_rec,high_vola,nvar,lags);
exp_low_sim = simulate_RFVAR(1e5,pphi_exp,low_vola,nvar,lags);
rec_high_corr = corr(rec_high_sim');
exp_low_corr = corr(exp_low_sim');
rec_high_cov = cov(rec_high_sim');
exp_low_cov = cov(exp_low_sim');


%% Discrete simulation for 1981 recession
Pphi_exp = reshape(pphi_exp,nvar,nvar,lags);
Pphi_rec = reshape(pphi_rec,nvar,nvar,lags);
peakdate = 130; % 1981:III
troughdate = 135; % 1982:IV
ylags_demean = zeros(nvar,lags);
ysim_demean = zeros(nvar,troughdate-peakdate);
exp_mean(1,1) = step1_Spec_Out.Coeff.S_Param{1,1}(1,1);
rec_mean(1,1) = step1_Spec_Out.Coeff.S_Param{1,1}(2,1);
exp_mean(2,1) = step1_Spec_Out.Coeff.S_Param{1,2}(1,1);
rec_mean(2,1) = step1_Spec_Out.Coeff.S_Param{1,2}(2,1);
exp_mean(3,1) = step1_Spec_Out.Coeff.S_Param{1,3}(1,1);
rec_mean(3,1) = step1_Spec_Out.Coeff.S_Param{1,3}(2,1);

for i_lag = 1:lags
	ylags_demean(:,i_lag) = data(peakdate-i_lag+1,:)'-exp_mean;
end
ysim_demean(:,1) =  Pphi_rec(:,:,1)*ylags_demean(:,1) ...
	+ Pphi_rec(:,:,2)*ylags_demean(:,2) ...
	%+ Pphi_rec(:,:,3)*ylags_demean(:,3) ...
	%+ Pphi_rec(:,:,4)*ylags_demean(:,4) ...
	;
for date = 2:troughdate-peakdate
	%ylags_demean(:,4) = ylags_demean(:,3);
	%ylags_demean(:,3) = ylags_demean(:,2);
	ylags_demean(:,2) = ylags_demean(:,1);
	ylags_demean(:,1) = ysim_demean(:,date-1); 
	ysim_demean(:,date) = Pphi_rec(:,:,1)*ylags_demean(:,1) ...
	+ Pphi_rec(:,:,2)*ylags_demean(:,2) ...
	%+ Pphi_rec(:,:,3)*ylags_demean(:,3) ...
	%+ Pphi_rec(:,:,4)*ylags_demean(:,4) ...
	;
end
ysim = ysim_demean + repmat(rec_mean,1,troughdate-peakdate);

data_seq = [data(peakdate:troughdate,:)]';
sim_seq = [data(peakdate,:)',ysim];

figure
subplot(3,1,1)
plot(0:troughdate-peakdate,[data(peakdate:troughdate,1)'],...
	0:troughdate-peakdate,[data(peakdate,1)' ysim(1,:)],'r--','LineWidth',lw,'MarkerSize',msz)
title('Sales Growth')
subplot(3,1,2)
plot(0:troughdate-peakdate,[data(peakdate:troughdate,2)'],...
	0:troughdate-peakdate,[data(peakdate,2)' ysim(2,:)],'r--','LineWidth',lw,'MarkerSize',msz)
ylabel('Percent')
title('GDP Growth')
subplot(3,1,3)
plot(0:troughdate-peakdate,[data(peakdate:troughdate,3)'],...
	0:troughdate-peakdate,[data(peakdate,3)' ysim(3,:)],'r--','LineWidth',lw,'MarkerSize',msz)
xlabel('Periods From Peak')
title('CIPI/GDP^{Pot}')
legend('Data','Discrete Shock')
print('../../figures/discrete_demeaned','-depsc2','-r300');

% revert to back to level
data_level = [rSalesGoods(peakdate+1:troughdate+1) rGDP(peakdate+1:troughdate+1) rCIPI(peakdate+1:troughdate+1)]';
sim_level = zeros(3,6);
sim_level(:,1) = data_level(:,1);
for date = 2:6
	sim_level(1:2,date) = (1+sim_seq(1:2,date)/100).*sim_level(1:2,date-1);
	sim_level(3,date) = sim_seq(3,date)/100.*GDPPOT(peakdate+1+date-1);
end

figure
subplot(3,1,1)
plot(0:troughdate-peakdate,data_level(1,:),...
	0:troughdate-peakdate,sim_level(1,:),'r--','LineWidth',lw,'MarkerSize',msz)
xticks(0:troughdate-peakdate)
title('Sales Level')

subplot(3,1,2)
plot(0:troughdate-peakdate,data_level(2,:),...
	0:troughdate-peakdate,sim_level(2,:),'r--','LineWidth',lw,'MarkerSize',msz)
xticks(0:troughdate-peakdate)
ylabel('Billions of 2009 Dollar')
title('GDP Level')

subplot(3,1,3)
plot(0:troughdate-peakdate,data_level(3,:),...
	0:troughdate-peakdate,sim_level(3,:),'r--','LineWidth',lw,'MarkerSize',msz)
xticks(0:troughdate-peakdate)
xlabel('Periods From Peak')
title('CIPI Level')
legend('Data','Discrete Shock')
print('../../figures/discrete_level','-depsc2','-r300');

% compute share
gdp_change = sim_level(2,end)-sim_level(2,1);
cipi_change = sim_level(3,end)-sim_level(3,1);
disp(cipi_change/gdp_change)

%% Discrete simulation for 1981 recession
Pphi_exp = reshape(pphi_exp,nvar,nvar,lags);
Pphi_rec = reshape(pphi_rec,nvar,nvar,lags);
peakdate = 135; % 1981:III
troughdate = 166; % 1982:IV
ylags_demean = zeros(nvar,lags);
ysim_demean = zeros(nvar,troughdate-peakdate);
exp_mean(1,1) = step1_Spec_Out.Coeff.S_Param{1,1}(1,1)+step1_Spec_Out.Coeff.S_Param{1,1}(1,2);
rec_mean(1,1) = step1_Spec_Out.Coeff.S_Param{1,1}(2,1)+step1_Spec_Out.Coeff.S_Param{1,1}(2,2);
exp_mean(2,1) = step1_Spec_Out.Coeff.S_Param{1,2}(1,1)+step1_Spec_Out.Coeff.S_Param{1,2}(1,2);
rec_mean(2,1) = step1_Spec_Out.Coeff.S_Param{1,2}(2,1)+step1_Spec_Out.Coeff.S_Param{1,2}(2,2);
exp_mean(3,1) = step1_Spec_Out.Coeff.S_Param{1,3}(1,1)+step1_Spec_Out.Coeff.S_Param{1,3}(1,2);
rec_mean(3,1) = step1_Spec_Out.Coeff.S_Param{1,3}(2,1)+step1_Spec_Out.Coeff.S_Param{1,3}(2,2);

for i_lag = 1:lags
	ylags_demean(:,i_lag) = data(peakdate-i_lag+1,:)'-rec_mean;
end
ysim_demean(:,1) =  Pphi_exp(:,:,1)*ylags_demean(:,1) ...
	+ Pphi_exp(:,:,2)*ylags_demean(:,2) ...
	%+ Pphi_rec(:,:,3)*ylags_demean(:,3) ...
	%+ Pphi_rec(:,:,4)*ylags_demean(:,4) ...
	;
for date = 2:troughdate-peakdate
	%ylags_demean(:,4) = ylags_demean(:,3);
	%ylags_demean(:,3) = ylags_demean(:,2);
	ylags_demean(:,2) = ylags_demean(:,1);
	ylags_demean(:,1) = ysim_demean(:,date-1); 
	ysim_demean(:,date) = Pphi_exp(:,:,1)*ylags_demean(:,1) ...
	+ Pphi_exp(:,:,2)*ylags_demean(:,2) ...
	%+ Pphi_rec(:,:,3)*ylags_demean(:,3) ...
	%+ Pphi_rec(:,:,4)*ylags_demean(:,4) ...
	;
end
ysim = ysim_demean + repmat(exp_mean,1,troughdate-peakdate);

data_seq = [data(peakdate:troughdate,:)]';
sim_seq = [data(peakdate,:)',ysim];

figure
subplot(3,1,1)
plot(0:troughdate-peakdate,[data(peakdate:troughdate,1)'],...
	0:troughdate-peakdate,[data(peakdate,1)' ysim(1,:)],'r--','LineWidth',lw,'MarkerSize',msz)
title('Sales Growth')
subplot(3,1,2)
plot(0:troughdate-peakdate,[data(peakdate:troughdate,2)'],...
	0:troughdate-peakdate,[data(peakdate,2)' ysim(2,:)],'r--','LineWidth',lw,'MarkerSize',msz)
ylabel('Percent')
title('GDP Growth')
subplot(3,1,3)
plot(0:troughdate-peakdate,[data(peakdate:troughdate,3)'],...
	0:troughdate-peakdate,[data(peakdate,3)' ysim(3,:)],'r--','LineWidth',lw,'MarkerSize',msz)
xlabel('Periods From Peak')
title('CIPI/GDP^{Pot}')
legend('Data','Discrete Shock')
print('../../figures/discrete_demeaned_t2p','-depsc2','-r300');

% revert to back to level
data_level = [rSalesGoods(peakdate+1:troughdate+1) rGDP(peakdate+1:troughdate+1) rCIPI(peakdate+1:troughdate+1)]';
sim_level = zeros(3,length(sim_seq));
sim_level(:,1) = data_level(:,1);
for date = 2:length(sim_seq)
	sim_level(1:2,date) = (1+sim_seq(1:2,date)/100).*sim_level(1:2,date-1);
	sim_level(3,date) = sim_seq(3,date)/100.*GDPPOT(peakdate+1+date-1);
end

figure
subplot(3,1,1)
plot(0:troughdate-peakdate,data_level(1,:),...
	0:troughdate-peakdate,sim_level(1,:),'r--','LineWidth',lw,'MarkerSize',msz)
xticks(0:troughdate-peakdate)
title('Sales Level')

subplot(3,1,2)
plot(0:troughdate-peakdate,data_level(2,:),...
	0:troughdate-peakdate,sim_level(2,:),'r--','LineWidth',lw,'MarkerSize',msz)
xticks(0:troughdate-peakdate)
ylabel('Billions of 2009 Dollar')
title('GDP Level')

subplot(3,1,3)
plot(0:troughdate-peakdate,data_level(3,:),...
	0:troughdate-peakdate,sim_level(3,:),'r--','LineWidth',lw,'MarkerSize',msz)
xticks(0:troughdate-peakdate)
xlabel('Periods From Peak')
title('CIPI Level')
%legend('Data','Discrete Shock')
print('../../figures/discrete_level_t2p','-depsc2','-r300');

% compute share
gdp_change = sim_level(2,end)-sim_level(2,1);
cipi_change = sim_level(3,end)-sim_level(3,1);
disp(cipi_change/gdp_change)
data_level(3,end)-data_level(3,1)

%% OIRF for all, fixed at low volatility regime
irf_periods = 10;
sales_shock = 1;
gdp_shock = 2;
cipi_shock = 3;

turbulent = zeros(nvar,nvar,irf_periods);
quite = zeros(nvar,nvar,irf_periods);

turbulent(:,sales_shock,:) = OIRF_RFVAR(sales_shock,irf_periods,pphi_rec,low_vola,nvar,lags);
quite(:,sales_shock,:) = OIRF_RFVAR(sales_shock,irf_periods,pphi_exp,low_vola,nvar,lags);

turbulent(:,gdp_shock,:) = OIRF_RFVAR(gdp_shock,irf_periods,pphi_rec,low_vola,nvar,lags);
quite(:,gdp_shock,:) = OIRF_RFVAR(gdp_shock,irf_periods,pphi_exp,low_vola,nvar,lags);

turbulent(:,cipi_shock,:) = OIRF_RFVAR(cipi_shock,irf_periods,pphi_rec,low_vola,nvar,lags);
quite(:,cipi_shock,:) = OIRF_RFVAR(cipi_shock,irf_periods,pphi_exp,low_vola,nvar,lags);

figure
subplot(3,3,1)
plot(0:irf_periods-1,squeeze(quite(1,1,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(1,1,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
title('Sales Shock')
ylabel('Sales Response')

subplot(3,3,2)
plot(0:irf_periods-1,squeeze(quite(1,2,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(1,2,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
title('GDP Shock')

subplot(3,3,3)
plot(0:irf_periods-1,squeeze(quite(1,3,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(1,3,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
title('CIPI Shock')

subplot(3,3,4)
plot(0:irf_periods-1,squeeze(quite(2,1,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(2,1,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
ylabel('GDP Response')


subplot(3,3,5)
plot(0:irf_periods-1,squeeze(quite(2,2,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(2,2,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)

subplot(3,3,6)
plot(0:irf_periods-1,squeeze(quite(2,3,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(2,3,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)

subplot(3,3,7)
plot(0:irf_periods-1,squeeze(quite(3,1,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(3,1,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
ylabel('CIPI Response')

subplot(3,3,8)
plot(0:irf_periods-1,squeeze(quite(3,2,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(3,2,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
xlabel('Periods from Shock')

subplot(3,3,9)
plot(0:irf_periods-1,squeeze(quite(3,3,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(3,3,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
legend('Expansion','Recession')

print('../../figures/9_OIRF_lowvola','-depsc2','-r300');

%% OIRF for all, fixed at high volatility regime
irf_periods = 10;
sales_shock = 1;
gdp_shock = 2;
cipi_shock = 3;

turbulent = zeros(nvar,nvar,irf_periods);
quite = zeros(nvar,nvar,irf_periods);

turbulent(:,sales_shock,:) = OIRF_RFVAR(sales_shock,irf_periods,pphi_rec,high_vola,nvar,lags);
quite(:,sales_shock,:) = OIRF_RFVAR(sales_shock,irf_periods,pphi_exp,high_vola,nvar,lags);

turbulent(:,gdp_shock,:) = OIRF_RFVAR(gdp_shock,irf_periods,pphi_rec,high_vola,nvar,lags);
quite(:,gdp_shock,:) = OIRF_RFVAR(gdp_shock,irf_periods,pphi_exp,high_vola,nvar,lags);

turbulent(:,cipi_shock,:) = OIRF_RFVAR(cipi_shock,irf_periods,pphi_rec,high_vola,nvar,lags);
quite(:,cipi_shock,:) = OIRF_RFVAR(cipi_shock,irf_periods,pphi_exp,high_vola,nvar,lags);

figure
subplot(3,3,1)
plot(0:irf_periods-1,squeeze(quite(1,1,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(1,1,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
title('Sales Shock')
ylabel('Sales Response')

subplot(3,3,2)
plot(0:irf_periods-1,squeeze(quite(1,2,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(1,2,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
title('GDP Shock')

subplot(3,3,3)
plot(0:irf_periods-1,squeeze(quite(1,3,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(1,3,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
title('CIPI/POTGDP Shock')

subplot(3,3,4)
plot(0:irf_periods-1,squeeze(quite(2,1,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(2,1,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
ylabel('GDP Response')


subplot(3,3,5)
plot(0:irf_periods-1,squeeze(quite(2,2,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(2,2,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)

subplot(3,3,6)
plot(0:irf_periods-1,squeeze(quite(2,3,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(2,3,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)

subplot(3,3,7)
plot(0:irf_periods-1,squeeze(quite(3,1,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(3,1,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
ylabel('CIPI/POTGDP Response')

subplot(3,3,8)
plot(0:irf_periods-1,squeeze(quite(3,2,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(3,2,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
xlabel('Periods from Shock')

subplot(3,3,9)
plot(0:irf_periods-1,squeeze(quite(3,3,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(3,3,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
legend('Expansion','Recession')

print('../../figures/9_OIRF_highvola','-depsc2','-r300');

%% RF-IRF for all,
irf_periods = 10;
sales_shock = 1;
gdp_shock = 2;
cipi_shock = 3;

turbulent = zeros(nvar,nvar,irf_periods);
quite = zeros(nvar,nvar,irf_periods);

turbulent(:,sales_shock,:) = OIRF_RFVAR(sales_shock,irf_periods,pphi_rec,eye(nvar),nvar,lags);
quite(:,sales_shock,:) = OIRF_RFVAR(sales_shock,irf_periods,pphi_exp,eye(nvar),nvar,lags);

turbulent(:,gdp_shock,:) = OIRF_RFVAR(gdp_shock,irf_periods,pphi_rec,eye(nvar),nvar,lags);
quite(:,gdp_shock,:) = OIRF_RFVAR(gdp_shock,irf_periods,pphi_exp,eye(nvar),nvar,lags);

turbulent(:,cipi_shock,:) = OIRF_RFVAR(cipi_shock,irf_periods,pphi_rec,eye(nvar),nvar,lags);
quite(:,cipi_shock,:) = OIRF_RFVAR(cipi_shock,irf_periods,pphi_exp,eye(nvar),nvar,lags);

figure
subplot(3,3,1)
plot(0:irf_periods-1,squeeze(quite(1,1,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(1,1,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
title('Sales Shock')
ylabel('Sales Response')

subplot(3,3,2)
plot(0:irf_periods-1,squeeze(quite(1,2,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(1,2,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
title('GDP Shock')

subplot(3,3,3)
plot(0:irf_periods-1,squeeze(quite(1,3,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(1,3,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
title('CIPI/POTGDP Shock')

subplot(3,3,4)
plot(0:irf_periods-1,squeeze(quite(2,1,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(2,1,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
ylabel('GDP Response')


subplot(3,3,5)
plot(0:irf_periods-1,squeeze(quite(2,2,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(2,2,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)

subplot(3,3,6)
plot(0:irf_periods-1,squeeze(quite(2,3,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(2,3,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)

subplot(3,3,7)
plot(0:irf_periods-1,squeeze(quite(3,1,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(3,1,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
ylabel('CIPI/POTGDP Response')

subplot(3,3,8)
plot(0:irf_periods-1,squeeze(quite(3,2,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(3,2,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
xlabel('Periods from Shock')

subplot(3,3,9)
plot(0:irf_periods-1,squeeze(quite(3,3,:)),'b-',...
	0:irf_periods-1,squeeze(turbulent(3,3,:)),'r-.',...
	'LineWidth',lw,'MarkerSize',msz)
legend('Expansion','Recession')

print('../../figures/9_IRF_eye','-depsc2','-r300');

%% OIRF for all, fixed at high volatility regime
irf_periods = 10;
sales_shock = 1;
gdp_shock = 2;
cipi_shock = 3;

turbulent = zeros(nvar,nvar,irf_periods);
quite = zeros(nvar,nvar,irf_periods);

turbulent(:,sales_shock,:) = OIRF_RFVAR(sales_shock,irf_periods,pphi_rec,high_vola,nvar,lags);
quite(:,sales_shock,:) = OIRF_RFVAR(sales_shock,irf_periods,pphi_exp,high_vola,nvar,lags);

turbulent(:,gdp_shock,:) = OIRF_RFVAR(gdp_shock,irf_periods,pphi_rec,high_vola,nvar,lags);
quite(:,gdp_shock,:) = OIRF_RFVAR(gdp_shock,irf_periods,pphi_exp,high_vola,nvar,lags);

turbulent(:,cipi_shock,:) = OIRF_RFVAR(cipi_shock,irf_periods,pphi_rec,high_vola,nvar,lags);
quite(:,cipi_shock,:) = OIRF_RFVAR(cipi_shock,irf_periods,pphi_exp,high_vola,nvar,lags);

figure
subplot(1,2,1)
plot(0:irf_periods-1,-squeeze(quite(2,1,:)),'b-',...
	0:irf_periods-1,-squeeze(quite(3,1,:)),'b-.',...
	0:irf_periods-1,zeros(1,irf_periods),'k.-',...
	'LineWidth',lw,'MarkerSize',msz)
title('Expansion')
legend('GDP','CIPI/POTGDP')
ylabel('Percent')
xlabel('Periods From Shock')

subplot(1,2,2)
plot(0:irf_periods-1,-squeeze(turbulent(2,1,:)),'r-',...
	0:irf_periods-1,-squeeze(turbulent(3,1,:)),'r-.',...
	0:irf_periods-1,zeros(1,irf_periods),'k.-',...
	'LineWidth',lw,'MarkerSize',msz)
legend('GDP','CIPI/POTGDP')
xlabel('Periods From Shock')
title('Recession')

print('../../figures/focus_salesshock','-depsc2','-r300');

%%
rmpath('m_Files');
rmpath('data_Files');
diary off

