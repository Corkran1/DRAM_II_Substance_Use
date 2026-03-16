%% Load weekly data
drugFree = readtable('Inpatient_4_Blue_DrugUsers_weekly_summary.csv');
drugUsers = readtable('Inpatient_4_Blue_DrugUsers_weekly_summary.csv');

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

%% Define parameters: {name, init, min, max}
params = {
    {'beta_D',     0.01, 0, 0.30}
    {'beta_F',     0.01, 0, 0.30}
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
model.sigma2 = 1;
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

%% Run DRAM
[results, chain, s2chain] = mcmcrun(model, data, params, options);

%% Extract parameter names
results.names = cellfun(@(c) c{1}, params, 'UniformOutput', false);
results.class = 'MCMC';

%% Save .mat results
save('Inpatient_4_Blue.mat', 'results', 'chain', 's2chain');

%% Save chain as CSV
chain_table = array2table(chain, 'VariableNames', results.names);
chain_table.sigma2 = s2chain;
writetable(chain_table, 'Inpatient_4_Blue_chain.csv');

%% Plot diagnostics
mcmc_diagnostics_plot(chain, results);
