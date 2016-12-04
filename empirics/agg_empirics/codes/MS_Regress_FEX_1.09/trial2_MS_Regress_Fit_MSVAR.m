% Example Script MS_Regress_Fit.m - MS-VAR estimation

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');

imported=importdata('./data_Files/trial2.txt');  % load some Data.
data = imported.data;
start = 0;
shift = 3;
forward_stock = data(1+shift:end,3);
GDP_match = data(1:end-shift,2);
sales_match = data(1:end-shift,1);

% data label
datelabel = (1947.25:0.25:2016.25)'; % because first differenced
yearnum = floor(datelabel);
monthnum = 12*(datelabel - yearnum)+2;
date_serial = datenum(yearnum,monthnum,ones(size(yearnum)));
break_dummy = (datelabel >= 1984.00);

data = 100*diff([GDP_match forward_stock]);

dep=data;                  % Defining dependent variables in system
nLag=2;                             % Number of lags in system
k=2;                                % Number of States
doIntercept=1;                      % add intercept to equations?
advOpt.distrib='Normal';            % The Distribution assumption (only 'Normal' for MS VAR models)
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details
advOpt.diagCovMat=0;                % since it reduced form, diagonal is stupid
advOpt.useMex=1;                % uses mex version of hamilton filter
advOpt.optimizer='fminsearch';     % use fmincon instead

[Spec_Out]=MS_VAR_Fit(dep,nLag,k,doIntercept,advOpt);

save trials2.mat

rmpath('m_Files');
rmpath('data_Files'); 

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

%% look at OIRF
Ssigma_regime1 = Spec_Out.Coeff.covMat{1};
Ssigma_regime2 = Spec_Out.Coeff.covMat{2};
L_regime1 = chol(Ssigma_regime1,'lower');
L_regime2 = chol(Ssigma_regime2,'lower');