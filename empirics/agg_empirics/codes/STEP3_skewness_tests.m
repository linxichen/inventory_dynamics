% This scripts computes skewness for test of deepness and steepness (Xia
% Table)
%% housekeeping
clear all;
close all;
clc;
global time;
clc;
load ../data/Prepared_FRED.mat;
addpath(genpath('helpers'));

%% computes skewness for each variable
%  for ISratio level
format bank
disp('========== IS ratio =============');
raw = log(rISratio);
[~,HP1] = hpfilter(raw,1600);
[~,HP2] = hpfilter(raw,1e5);
BP1 = bandpass(raw,6,32);
BP2 = bandpass(raw,2,80);
FD = diff(raw);

disp('---------------------------')
disp('HP(1600)');
x = HP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('HP(1e5)');
x = HP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(6,32)');
x = BP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(2,80)');
x = BP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%  for ISratio GROWTH
disp('---------------------------')
disp('FD');
x = FD;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%% for inventory stock
format bank
disp('========== Inventory Stock =============');
raw = log(rInventory);
[~,HP1] = hpfilter(raw,1600);
[~,HP2] = hpfilter(raw,1e5);
BP1 = bandpass(raw,6,32);
BP2 = bandpass(raw,2,80);
FD = diff(raw);

disp('---------------------------')
disp('HP(1600)');
x = HP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('HP(1e5)');
x = HP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(6,32)');
x = BP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(2,80)');
x = BP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%  for ISratio GROWTH
disp('---------------------------')
disp('FD');
x = FD;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%% for CIPI/GDP
format bank
disp('========== CIPI/GDP =============');
raw = share_CIPI;
[~,HP1] = hpfilter(raw,1600);
[~,HP2] = hpfilter(raw,1e5);
BP1 = bandpass(raw,6,32);
BP2 = bandpass(raw,2,80);
FD = diff(raw);

disp('---------------------------')
disp('HP(1600)');
x = HP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('HP(1e5)');
x = HP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(6,32)');
x = BP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(2,80)');
x = BP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%  for steepness
disp('---------------------------')
disp('FD');
x = FD;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%% for sales
disp('========== Sales =============');
raw = log(rSalesGoods);
[~,HP1] = hpfilter(raw,1600);
[~,HP2] = hpfilter(raw,1e5);
BP1 = bandpass(raw,6,32);
BP2 = bandpass(raw,2,80);
FD = diff(raw);

disp('---------------------------')
disp('HP(1600)');
x = HP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('HP(1e5)');
x = HP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(6,32)');
x = BP1;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

disp('---------------------------')
disp('BP(2,80)');
x = BP2;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

%  for deepness
disp('---------------------------')
disp('FD');
x = FD;
[s t] = skewTestAC(x);
skewCoeffs = s;
testStats = t;
pVals = pnorm(t);
fprintf('%.3f [%.3f]\n',skewCoeffs,pVals);

