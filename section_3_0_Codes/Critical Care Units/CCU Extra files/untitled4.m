%% fit_mrsa_dram_full.m
% Full DRAM MCMC workflow for MRSA transmission model

clear; clc;

%% Load weekly data
drugFree = readtable('4_Gold_Observation_DrugFree_weekly_summary.csv');
drugUsers = readtable('4_Gold_Observation_DrugUsers_weekly_summary.csv');

weeks = (0:height(drugFree)-1)';

data.time = weeks;
data.S_F = drugFree.weekly_susceptible;
data.I_F = drugFree.weekly_infected;
data.R_F = drugFree.weekly_recovered;
data.S_D = drugUsers.weekly_susceptible;
data.I_D = drugUsers.weekly_infected;
data.R_D = drugUsers.weekly_recovered;

%% Initial conditions vector for ODE solver
Y0 = [data.S_D(1), data.I_D(1), data.R_D(1), ...
      data.S_F(1), data.I_F(1), data.R_F(1), ...
      50, 5];  % [S_D, I_D, R_D, S_F, I_F, R_F, H_S, H_C]

%% Define parameters as struct array
params = struct( ...
  'name', {'lambda', 'gamma', 'mu_D', 'mu_F', 'delta', 'sigma', ...
           'Lambda_S_D', 'Lambda_I_D', 'Lambda_S_F', 'Lambda_I_F'}, ...
  'init', [0.01, 0.1, 0.01, 0.01, 0.1, 0.1, 10, 1, 10, 1], ...
  'min',  [0,    0,   0,    0,    0,   0,   0,  0,  0,  0], ...
  'max',  [1,    1,   1,    1,    1,   1, 100, 10, 100, 10]);

%% Model struct for mcmcrun
model.ssfun = @(theta, data) ssfun_mrsa(theta, data, Y0);
model.sigma2 = 1;              % initial guess of observation noise variance
model.N = length(params);      % number of parameters to estimate

%% MCMC options
options = mcmcset(...
    'nsimu', 5000, ...
    'burnintime', 1000, ...
    'adaptint', 100, ...
    'drscale', 2, ...
    'method', 'dram', ...
    'updatesigma', 1, ...
    'verbosity', 1);

disp(params);
disp(class(params));
disp(isstruct(params));

%% Run DRAM MCMC
[results, chain, s2chain] = mcmcrun(model, data, params, options);

%% Add parameter names and class for plotting
results.names = {params.name};
results.class = 'MCMC';

%% Save results
save('mcmc_mrsa_results.mat', 'results', 'chain', 's2chain');

%% Plot MCMC diagnostics (using custom function)
mcmc_diagnostics_plot(chain, results);

disp('MCMC run complete.');
