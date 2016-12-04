%% housekeeping
clear;
clc;
close all
addpath(genpath('~/Dropbox/matlabtools/BVARtoolbox'))

load('../data/manual_select_dummy_FRED')
y_table = [diff(ln_rSales) diff(ln_rGDP) share_rCIPI_potential(2:end,:)];
regimes = regimes(2:end,:);

%% specify the empirical model
N = 3; p = 3; M = max(regimes);
model.N = N;
model.p = p;
model.K = model.N*model.p+1;
model.T = length(regimes);

%% specify priors
priors.mmu_mean = zeros(N*M,1);
priors.mmu_cov = 1e5*eye(N*M);
priors.pphi_mean = zeros(N*N*p,1);
priors.pphi_cov = 1e5*eye(N*N*p);

%% specify options for the Gibbs sampler
options.burnin = 2e4;
options.R = 5e4;
tic
draws = dummyVAR_Gibbs(y_table,regimes,model,priors,options);
toc

save dummyVAR_Gibbs
%% look at result
median(draws.Ssigma_array,4)
mean(draws.pphi,3)
mean(draws.mmu,3)