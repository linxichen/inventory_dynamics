% Example Script MS_Regress_Fit.m - MS-VAR estimation

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');

imported=importdata('./data_Files/trial2.txt');  % load some Data.
data = imported.data;
shift = 3;
forward_stock = data(1+shift:end,3);
GDP_match = data(1:end-shift,2);
sales_match = data(1:end-shift,1);

data = 100*diff([GDP_match forward_stock]);

dep=data;                  % Defining dependent variables in system
nLag=2;                             % Number of lags in system
k=2;                                % Number of States
doIntercept=1;                      % add intercept to equations?
advOpt.distrib='Normal';            % The Distribution assumption (only 'Normal' for MS VAR models)
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details
advOpt.diagCovMat=0;                % since it reduced form, diagonal is stupid
advOpt.useMex=1;                % uses mex version of hamilton filter
advOpt.optimizer='fmincon';     % use fmincon instead

[Spec_Out]=MS_VAR_Fit(dep,nLag,k,doIntercept,advOpt);

save trials2.mat

rmpath('m_Files');
rmpath('data_Files'); 