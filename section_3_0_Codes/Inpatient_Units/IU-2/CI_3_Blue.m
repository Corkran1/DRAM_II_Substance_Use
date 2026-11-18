% compute_credible_intervals.m
% Assumes 'chain_table' and 'results.names' are in your workspace

% If not already loaded, load them:
load('Inpatient_3_Blue.mat', 'results', 'chain');

% Rebuild chain_table
param_names = results.names;
chain_table = array2table(chain, 'VariableNames', param_names);

% Initialize summary table
summary = table('Size', [length(param_names), 5], ...
                'VariableTypes', {'string', 'double', 'double', 'double', 'double'}, ...
                'VariableNames', {'Parameter', 'Mean', 'StdDev', 'CI_2_5', 'CI_97_5'});

% Compute statistics for each parameter
for i = 1:length(param_names)
    pname = param_names{i};
    samples = chain_table.(pname);
    
    summary.Parameter(i) = pname;
    summary.Mean(i) = mean(samples);
    summary.StdDev(i) = std(samples);
    ci = quantile(samples, [0.025, 0.975]);
    summary.CI_2_5(i) = ci(1);
    summary.CI_97_5(i) = ci(2);
end

% Save to CSV
writetable(summary, 'Inpatient_3_Blue_CredibleIntervals.csv');

% Display
disp('Credible intervals summary:');
disp(summary);
