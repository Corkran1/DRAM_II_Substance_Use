

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

%% Initial Conditions
Y0 = [data.S_D(1), data.I_D(1), data.R_D(1), ...
      data.S_F(1), data.I_F(1), data.R_F(1), ...
      50, 5];  % [S_D, I_D, R_D, S_F, I_F, R_F, H_S, H_C]

%% Define parameters as cell array of cells (name, init, lower, upper)
params = {
    {'lambda',     0.01, 0, 1}
    {'gamma',      0.1,  0, 1}
    {'mu_D',       0.01, 0, 1}
    {'mu_F',       0.01, 0, 1}
    {'delta',      0.1,  0, 1}
    {'sigma',      0.1,  0, 1}
    {'Lambda_S_D', 10,   0, 100}
    {'Lambda_I_D', 1,    0, 10}
    {'Lambda_S_F', 10,   0, 100}
    {'Lambda_I_F', 1,    0, 10}
};

%% Model struct for mcmcrun
model.ssfun = @(theta, data) ssfun_mrsa(theta, data, Y0);
model.sigma2 = 1;      % initial guess for observation noise
model.N = length(params);

%% MCMC Options
options = mcmcset(...
    'nsimu', 5000, ...
    'burnintime', 1000, ...
    'adaptint', 100, ...
    'drscale', 2, ...
    'method', 'dram', ...
    'updatesigma', 1, ...
    'verbosity', 1);

%% Run DRAM MCMC
[results, chain, s2chain] = mcmcrun(model, data, params, options);

%% Add parameter names for plotting

results.names = {params.name};
results.class = 'MCMC';
%% Save results
save('mcmc_mrsa_results.mat', 'results', 'chain', 's2chain');

%% Plot results
mcmc_diagnostics_plot(chain, results);
disp('MCMC run complete.');
