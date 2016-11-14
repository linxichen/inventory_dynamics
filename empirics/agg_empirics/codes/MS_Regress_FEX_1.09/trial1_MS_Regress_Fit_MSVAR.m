% Example Script MS_Regress_Fit.m - MS-VAR estimation

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path
addpath('data_Files');

imported=importdata('./data_Files/trial1.txt');  % load some Data.
data = imported.data;
%data = 100*diff(data(:,2:3));
shift = 3;
forward_resid_rCIPI = data(130+2+shift:end,2);
match_GDP = data(130+1:end-shift,3);
data = 100*[forward_resid_rCIPI diff(match_GDP)];

dep=data;                  % Defining dependent variables in system
nLag=1;                             % Number of lags in system
k=2;                                % Number of States
doIntercept=1;                      % add intercept to equations?
advOpt.distrib='Normal';            % The Distribution assumption (only 'Normal' for MS VAR models)
advOpt.std_method=2;                % Defining the method for calculation of standard errors. See pdf file for more details
advOpt.diagCovMat=0;                % since it reduced form, diagonal is stupid
advOpt.useMex=1;                % uses mex version of hamilton filter
advOpt.optimizer='fmincon';     % use fmincon instead

[Spec_Out]=MS_VAR_Fit(dep,nLag,k,doIntercept,advOpt);

save trials1.mat

rmpath('m_Files');
rmpath('data_Files'); 