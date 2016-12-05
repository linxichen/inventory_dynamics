% Example Script MS_Regress_Fit.m

clear;
close all;
clc;

diary('deliverable_1_7.txt')
diary on

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');

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
lags = 0;                            % decide on lags
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
advOpt.std_method=2;                % Defining the method for calculation of standard errors. See pdf file for more details
advOpt.diagCovMat=1;
advOpt.useMEX=1;
advOpt.optimizer='fminsearch';


[step1_Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model
save deliverable1_7.mat

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
lags = 3;
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

%% IRF
irf_periods = 24;
conditonal_mean_tur(1,1) = step1_Spec_Out.Coeff.S_Param{1}(1,2);
conditonal_mean_tur(2,1) = step1_Spec_Out.Coeff.S_Param{2}(1,2);
conditonal_mean_tur(3,1) = step1_Spec_Out.Coeff.S_Param{3}(1,2);

conditonal_mean_qui(1,1) = step1_Spec_Out.Coeff.S_Param{1}(1,1);
conditonal_mean_qui(2,1) = step1_Spec_Out.Coeff.S_Param{2}(1,1);
conditonal_mean_qui(3,1) = step1_Spec_Out.Coeff.S_Param{3}(1,1);


turbulent = OIRF_RFVAR(1,irf_periods,pphi_rec,step2_Spec_Out.Coeff.covMat{2},nvar,lags);
% turbulent = turbulent + repmat(conditonal_mean_tur,1,size(turbulent,2));
quite = OIRF_RFVAR(1,irf_periods,pphi_exp,step2_Spec_Out.Coeff.covMat{2},nvar,lags);
% quite = quite + repmat(conditonal_mean_qui,1,size(quite,2));

turbulent_cumsum = cumsum(turbulent,2);
quite_cumsum = cumsum(quite,2);

quite_cumsum(2,6)
quite(3,6)-quite(3,1)

turbulent_cumsum(2,6)
turbulent(3,6)-turbulent(3,1)

figure
subplot(1,2,1)
plot(quite')
ylim([-0.2 1.2])

subplot(1,2,2)
plot(turbulent')
ylim([-0.2 1.2])

figure
subplot(1,2,1)
plot(quite(3,:)')
ylim([-0.05 0.15])
subplot(1,2,2)
plot(turbulent(3,:)')
ylim([-0.05 0.15])

%% Discrete simulation



%%
rmpath('m_Files');
rmpath('data_Files');
diary off

