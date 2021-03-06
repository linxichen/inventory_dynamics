% Example Script MS_Regress_Fit.m

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');

% load data
imported=importdata('./data_Files/trial1.txt');  % load some Data.
data = imported.data;
%data = 100*diff(data(:,2:3));
start = 0; % first start date
shift = 0; % shift some of the variable forward or backward
forward_resid_rCIPI = data(start+2+shift:end,2);
forward_IS = data(start+1+shift:end,4);
match_GDP = data(start+1:end-shift,3);
match_sales = data(start+1:end-shift,1);
data = 100*[diff(match_sales) forward_resid_rCIPI diff(match_GDP) diff(forward_IS)];

% data label
datelabel = (1947.25:0.25:2016.25)'; % because first differenced
yearnum = floor(datelabel);
monthnum = 12*(datelabel - yearnum)+2;
date_serial = datenum(yearnum,monthnum,ones(size(yearnum)));
break_dummy = (datelabel >= 1984.00);

% crated dependent vars and regressors
lags = 2;                            % decide on lags
Y = data;
nvar = size(Y,2);
YLAG = lagmatrix(Y,1:lags);
dep=Y;                  % Defining dependent variable from .mat file
constVec=ones(length(dep),1);       % Defining a constant vector in mean equation (just an example of how to do it)
indep{1}=[constVec constVec.*break_dummy];                  % Defining some explanatory variables
indep{2}=[constVec constVec.*break_dummy];                  % Defining some explanatory variables
indep{3}=[constVec constVec.*break_dummy];                  % Defining some explanatory variables
indep{4}=[constVec constVec.*break_dummy];                  % Defining some explanatory variables

k=9;                                % Number of States
S{1}=[1 1 1];                        % Defining which parts of the equation will switch states (column 1 and variance only)
S{2}=[1 1 1];                        % Defining which parts of the equation will switch states (column 1 and variance only)
S{3}=[1 1 1];                        % Defining which parts of the equation will switch states (column 1 and variance only)
S{4}=[1 1 1];                        % Defining which parts of the equation will switch states (column 1 and variance only)

advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details
advOpt.diagCovMat=0;
advOpt.useMex=1;


[Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model
save trial4.mat

%% plot regimes
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
irfperiods = 20;
impulsevar = 2;
after = 1;

% initialzation
impulsevec = zeros(4,irfperiods);
Yimpulse = zeros(4,irfperiods);

B = [Spec_Out.Coeff.S_Param{1,1}(:,regime)'; ...
	        Spec_Out.Coeff.S_Param{1,2}(:,regime)'; ...
			Spec_Out.Coeff.S_Param{1,3}(:,regime)'; ...
			Spec_Out.Coeff.S_Param{1,4}(:,regime)'; ...
			];
B0 = B(:,1);
B1 = B(:,2);
Ymean = B0+after*B1;
L = chol(Spec_Out.Coeff.covMat{regime})';
impulsevec(impulsevar,1) = 1;
Yimpulse(:,1) = B0 + B1*after + L*impulsevec(:,1);
for i_period = 2:irfperiods
	Yimpulse(:,i_period) = B0 + B1*after + L*impulsevec(:,i_period);
end
Yirf_regime1 = Yimpulse - repmat(Ymean,1,irfperiods);

regime = 2;
% initialzation
impulsevec = zeros(4,irfperiods);
Yimpulse = zeros(4,irfperiods);

B = [Spec_Out.Coeff.S_Param{1,1}(:,regime)'; ...
	        Spec_Out.Coeff.S_Param{1,2}(:,regime)'; ...
			Spec_Out.Coeff.S_Param{1,3}(:,regime)'; ...
			Spec_Out.Coeff.S_Param{1,4}(:,regime)'; ...
			];
B0 = B(:,1);
B1 = B(:,2);
Ymean = B0+after*B1;
L = chol(Spec_Out.Coeff.covMat{regime})';
impulsevec(impulsevar,1) = 1;
Yimpulse(:,1) = B0 + B1*after + L*impulsevec(:,1);
for i_period = 2:irfperiods
	Yimpulse(:,i_period) = B0 + B1*after + L*impulsevec(:,i_period);
end
Yirf_regime2 = Yimpulse - repmat(Ymean,1,irfperiods);

figure 
subplot(2,1,1)
plot(Yirf_regime1')
legend('resid_CIPI','GDP','ISratio')
subplot(2,1,2)
plot(Yirf_regime2')
legend('resid_CIPI','GDP','ISratio')

rmpath('m_Files');
rmpath('data_Files'); 

