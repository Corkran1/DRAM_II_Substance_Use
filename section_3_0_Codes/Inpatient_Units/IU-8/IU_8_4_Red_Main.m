%% IU-8 (4 Red): extended MRSA SIR/HCW DRAM model
% Run from this folder with the two weekly-summary CSVs available on path.
clearvars;
rng(20261118, 'twister'); % Reproducible starting seed; not a convergence guarantee.

%% Load and validate weekly data
drugFree = readtable('Inpatient_4_Red_DrugFree_weekly_summary.csv');
drugUsers = readtable('Inpatient_4_Red_DrugUsers_weekly_summary.csv');
assert(height(drugFree) == height(drugUsers), 'SU and NS series must have the same number of weeks.');
required = {'weekly_susceptible','weekly_infected','weekly_recovered'};
assert(all(ismember(required, drugFree.Properties.VariableNames)) && ...
       all(ismember(required, drugUsers.Properties.VariableNames)), 'Missing weekly summary columns.');
assert(height(drugFree) >= 2, 'At least two time points are required.');

data.time = (0:height(drugFree)-1)';
data.S_F = drugFree.weekly_susceptible;
data.I_F = drugFree.weekly_infected;
data.R_F = drugFree.weekly_recovered;
data.S_D = drugUsers.weekly_susceptible;
data.I_D = drugUsers.weekly_infected;
data.R_D = drugUsers.weekly_recovered;

%% Initial state: SU (D), NS (F), susceptible/contaminated HCWs
Y0 = [data.S_D(1), data.I_D(1), data.R_D(1), ...
      data.S_F(1), data.I_F(1), data.R_F(1), 50, 5];

%% Parameter order MUST match mrsa_ode.m
% Retain original unit-specific data-informed upper bounds (0.30)
% for both transmission coefficients. Prior justification belongs to E2.
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
parameter_names = cellfun(@(p) p{1}, params, 'UniformOutput', false);

%% DRAM settings (kept from supplied script)
model.ssfun = @(theta, observations) ssfun_mrsa(theta, observations, Y0);
model.sigma2 = 1;
model.N = length(params);
options = mcmcset('nsimu', 5000, 'burnintime', 1000, ...
    'adaptint', 100, 'drscale', 2, 'method', 'dram', ...
    'updatesigma', 1, 'verbosity', 1);

%% Fit
[results, chain, s2chain] = mcmcrun(model, data, params, options);
results.names = parameter_names;
results.class = 'MCMC';
save('IU_8_4_Red_results.mat', 'results', 'chain', 's2chain', ...
    'params', 'parameter_names', 'data', 'Y0', 'options');

%% Export raw chain (do not assume burn-in has or has not been removed)
chainTable = array2table(chain, 'VariableNames', parameter_names);
writetable(chainTable, 'IU_8_4_Red_chain.csv');
if exist('mcmc_diagnostics_plot', 'file') == 2
    mcmc_diagnostics_plot(chain, results);
else
    warning('mcmc_diagnostics_plot not found; chain and results were saved.');
end
disp('IU-8 (4 Red) DRAM run complete.');
