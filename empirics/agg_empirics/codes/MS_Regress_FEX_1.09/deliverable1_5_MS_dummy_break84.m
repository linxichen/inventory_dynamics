% Example Script MS_Regress_Fit.m

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');

% load data
load manual_select_dummy_FRED
%data = 100*diff(data(:,2:3));
start = 0; % first start date
shift = 0; % shift some of the variable forward or backward
% forward_resid_rCIPI = data(start+2+shift:end,2);
share_rCIPI_pot = share_rCIPI_potential;
% forward_stock = data(start+1+shift:end,3);
% match_GDP = ln_rGDP;
match_GDP = (rGDP(2:end) - rGDP(1:end-1))./GDPPOT(1:end-1);
match_sales = (rSalesGoods(2:end) - rSalesGoods(1:end-1))./GDPPOT(1:end-1);
data = 100*[match_sales match_GDP share_rCIPI_pot(2:end)];

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
S{2}=ones(1,2+lags*nvar+1);S{2}(1,end)=0;     % Defining which parts of the equation will switch states (column 1 and variance only)
S{3}=ones(1,2+lags*nvar+1);S{3}(1,end)=0;      % Defining which parts of the equation will switch states (column 1 and variance only)

advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=2;                % Defining the method for calculation of standard errors. See pdf file for more details
advOpt.diagCovMat=1;
advOpt.useMEX=1;
advOpt.optimizer='fminsearch';


[step1_Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model
save deliverable1_5.mat

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
save deliverable1_5.mat

%% plot regimes
safe_dates = date_serial(start+1+lags:end-shift,:);
figure
plot(safe_dates,step1_Spec_Out.smoothProb(lags+1:end,2));
h = gca;
datetick('x','yyyy','keepticks')
xlabel('Time');
ylabel('Smoothed States Probabilities');
legend('Recession Regimes');
title('Autoregressive Regime')
axis tight
recessband = recessionplot;

safe_dates = date_serial(start+1+lags:end-shift,:);
figure
plot(safe_dates,step2_Spec_Out.smoothProb(:,2));
h = gca;
datetick('x','yyyy','keepticks')
xlabel('Time');
ylabel('Probability');
legend('Moderation');
title('Smoothed States Probabilities , Moderation Covariance Regime')
axis tight
recessband = recessionplot;


%% IRF
regime = 1;
irfperiods = 20;
impulsevar = 2;

% initialzation
impulsevec = zeros(3,irfperiods);
Yimpulse = zeros(3,irfperiods);

B = [step1_Spec_Out.Coeff.S_Param{1,1}(:,regime)'; ...
	        step1_Spec_Out.Coeff.S_Param{1,2}(:,regime)'; ...
			step1_Spec_Out.Coeff.S_Param{1,3}(:,regime)'; ...
			];
B0 = B(:,1);
B1 = B(:,2:end);
Ymean = (eye(3) - B1)\B0;
L = chol(step1_Spec_Out.Coeff.covMat{regime})';
impulsevec(impulsevar,1) = 1;
Yimpulse(:,1) = B0 + B1*Ymean + L*impulsevec(:,1);
for i_period = 2:irfperiods
	Yimpulse(:,i_period) = B0 + B1*Yimpulse(:,i_period-1) + L*impulsevec(:,i_period-1);
end
Yirf_regime1 = Yimpulse - Ymean;

regime = 2;
% initialzation
impulsevec = zeros(3,irfperiods);
Yimpulse = zeros(3,irfperiods);

B = [step1_Spec_Out.Coeff.S_Param{1,1}(:,regime)'; ...
	        step1_Spec_Out.Coeff.S_Param{1,2}(:,regime)'; ...
			step1_Spec_Out.Coeff.S_Param{1,3}(:,regime)'; ...
			];
B0 = B(:,1);
B1 = B(:,2:end);
Ymean = (eye(3) - B1)\B0;
L = chol(step1_Spec_Out.Coeff.covMat{regime})';
impulsevec(impulsevar,1) = 1;
Yimpulse(:,1) = B0 + B1*Ymean + L*impulsevec(:,1);
for i_period = 2:irfperiods
	Yimpulse(:,i_period) = B0 + B1*Yimpulse(:,i_period-1) + L*impulsevec(:,i_period-1);
end
Yirf_regime2 = Yimpulse - Ymean;

figure 
subplot(2,1,1)
plot(Yirf_regime1')
legend('resid_CIPI','GDP','ISratio')
subplot(2,1,2)
plot(Yirf_regime2')
legend('resid_CIPI','GDP','ISratio')

rmpath('m_Files');
rmpath('data_Files'); 

