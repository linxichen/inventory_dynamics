% Example Script MS_Regress_Fit.m - MS-VAR estimation

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');

% load data
imported=importdata('./data_Files/trial1.txt');  % load some Data.
data = imported.data;
%data = 100*diff(data(:,2:3));
start = 0;
shift = 0;
forward_resid_rCIPI = data(start+2+shift:end,2);
forward_IS = data(start+1+shift:end,4);
match_sales = data(start+1:end-shift,1);
data = 100*[forward_resid_rCIPI diff(match_sales) diff(forward_IS)];

% data label
datelabel = (1947.25:0.25:2016.25)'; % because first differenced
yearnum = floor(datelabel);
monthnum = 12*(datelabel - yearnum)+2;
date_serial = datenum(yearnum,monthnum,ones(size(yearnum)));

dep=data;                  % Defining dependent variables in system
nLag=2;                             % Number of lags in system
k=4;                                % Number of States
doIntercept=1;                      % add intercept to equations?
advOpt.distrib='Normal';            % The Distribution assumption (only 'Normal' for MS VAR models)
advOpt.std_method=2;                % Defining the method for calculation of standard errors. See pdf file for more details
advOpt.diagCovMat=0;                % since it reduced form, diagonal is stupid
advOpt.useMex=1;                % uses mex version of hamilton filter
advOpt.optimizer='fminsearch';     % use fmincon instead

[Spec_Out]=MS_VAR_Fit(dep,nLag,k,doIntercept,advOpt);
save trials1.mat

%% plot
safe_dates = date_serial(start+1:end-shift,:);
figure
plot(safe_dates,Spec_Out.smoothProb(:,2));
h = gca;
datetick('x','yyyy','keepticks')
xlabel('Time');
ylabel('Smoothed States Probabilities');
legend('Regime 1');
axis tight
recessband = recessionplot;

%% IRF
regime = 1;
irfperiods = 6;
impulsevar = 1;

% initialzation
impulsevec = zeros(3,irfperiods);
Yimpulse = zeros(3,irfperiods);

B = [Spec_Out.Coeff.S_Param{1,1}(:,regime)'; ...
	        Spec_Out.Coeff.S_Param{1,2}(:,regime)'; ...
			Spec_Out.Coeff.S_Param{1,3}(:,regime)'; ...
			];
B0 = B(:,regime);
B1 = B(:,2:end);
Ymean = (eye(3) - B1)\B0;
L = chol(Spec_Out.Coeff.covMat{regime})';
impulsevec(impulsevar,1) = 1;
Yimpulse(:,1) = B0 + B1*Ymean + L*impulsevec(:,1);
for i_period = 2:irfperiods
	Yimpulse(:,i_period) = B0 + B1*Yimpulse(:,i_period-1) + L*impulsevec(:,i_period);
end
Yirf_regime1 = Yimpulse - Ymean;

regime = 2;
% initialzation
impulsevec = zeros(3,irfperiods);
Yimpulse = zeros(3,irfperiods);

B = [Spec_Out.Coeff.S_Param{1,1}(:,regime)'; ...
	        Spec_Out.Coeff.S_Param{1,2}(:,regime)'; ...
			Spec_Out.Coeff.S_Param{1,3}(:,regime)'; ...
			];
B0 = B(:,regime);
B1 = B(:,2:end);
Ymean = (eye(3) - B1)\B0;
L = chol(Spec_Out.Coeff.covMat{regime})';
impulsevec(impulsevar,1) = 1;
Yimpulse(:,1) = B0 + B1*Ymean + L*impulsevec(:,1);
for i_period = 2:irfperiods
	Yimpulse(:,i_period) = B0 + B1*Yimpulse(:,i_period-1) + L*impulsevec(:,i_period);
end
Yirf_regime2 = Yimpulse - Ymean;

figure 
subplot(2,1,1)
plot(Yirf_regime1')
legend('resid_CIPI','sales','ISratio')
subplot(2,1,2)
plot(Yirf_regime2')
legend('resid_CIPI','sales','ISratio')

%% aftermath

rmpath('m_Files');
rmpath('data_Files'); 