%% Progressive Care Unit: Extended MRSA SIR/HCW DRAM model
% Step 6: Pilot calibration and extended independent chains
% Run from the Progressive Care Unit folder.

clearvars;

%% Load and validate weekly data
drugFree = readtable(...
    'Progressive_Care_Unit_DrugFree_weekly_summary.csv');

drugUsers = readtable(...
    'Progressive_Care_Unit_DrugUsers_weekly_summary.csv');

assert(height(drugFree) == height(drugUsers), ...
    'SU and NS series must have the same number of weeks.');

required = {'weekly_susceptible', ...
            'weekly_infected', ...
            'weekly_recovered'};

assert(all(ismember(required, drugFree.Properties.VariableNames)) && ...
       all(ismember(required, drugUsers.Properties.VariableNames)), ...
       'Missing weekly summary columns.');

assert(height(drugFree) >= 2, ...
    'At least two time points are required.');

%% Assign observed compartments
data.time = (0:height(drugFree)-1)';

data.S_F = drugFree.weekly_susceptible;
data.I_F = drugFree.weekly_infected;
data.R_F = drugFree.weekly_recovered;

data.S_D = drugUsers.weekly_susceptible;
data.I_D = drugUsers.weekly_infected;
data.R_D = drugUsers.weekly_recovered;

%% Validate observations
SU = [data.S_D(:), data.I_D(:), data.R_D(:)];
NS = [data.S_F(:), data.I_F(:), data.R_F(:)];

assert(~isequal(SU, NS), ...
    'SU and NS observations are identical.');

assert(all(isfinite(SU(:))) && all(isfinite(NS(:))), ...
    'Nonfinite observations detected.');

assert(all(SU(:) >= 0) && all(NS(:) >= 0), ...
    'Negative observations detected.');

%% Initial conditions
Y0 = [data.S_D(1), data.I_D(1), data.R_D(1), ...
      data.S_F(1), data.I_F(1), data.R_F(1), ...
      50, 5];

%% Parameter definitions
% Parameter order must match mrsa_ode.m.
% Original Progressive Care Unit transmission bounds retained.
% Comparable-bound sensitivity fits will be conducted separately.

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

parameter_names = cellfun(@(p) p{1}, params, ...
    'UniformOutput', false);

npar = length(params);

lower_bounds = cellfun(@(p) p{3}, params);
upper_bounds = cellfun(@(p) p{4}, params);

%% DRAM model structure
model.ssfun = @(theta, observations) ...
    ssfun_mrsa(theta, observations, Y0);

% Corrected observation count
model.N = 6 * (numel(data.time) - 1);

% Preserve original residual variance prior
model.sigma2 = 1;
model.S20 = 1;
model.N0 = 1;

fprintf('Number of weekly observations: %d\n', ...
    numel(data.time));

fprintf('Number of residual observations: %d\n', ...
    model.N);

%% Numerical verification
fprintf('ODE: %s\n', which('mrsa_ode'));
fprintf('Objective: %s\n', which('ssfun_mrsa'));

theta0 = cellfun(@(p) p{2}, params);

baselineSSE = ssfun_mrsa(theta0, data, Y0);

assert(isfinite(baselineSSE) && ...
       baselineSSE >= 0 && baselineSSE < 1e12, ...
       'Invalid baseline SSE.');

thetaD = theta0;
thetaD(1) = 0.02;

thetaF = theta0;
thetaF(2) = 0.02;

sseD = ssfun_mrsa(thetaD, data, Y0);
sseF = ssfun_mrsa(thetaF, data, Y0);

assert(isfinite(sseD) && sseD >= 0 && sseD < 1e12, ...
    'Invalid beta_D perturbation SSE.');

assert(isfinite(sseF) && sseF >= 0 && sseF < 1e12, ...
    'Invalid beta_F perturbation SSE.');

assert(abs(sseD - baselineSSE) > 1e-8, ...
    'Objective function insensitive to beta_D perturbation.');

assert(abs(sseF - baselineSSE) > 1e-8, ...
    'Objective function insensitive to beta_F perturbation.');

fprintf('Baseline SSE: %.12g\n', baselineSSE);
fprintf('beta_D perturbation SSE: %.12g\n', sseD);
fprintf('beta_F perturbation SSE: %.12g\n', sseF);

disp('Numerical verification passed.');

%% Pilot DRAM simulation
fprintf('\n====================================\n');
fprintf('Progressive Care Unit: Starting pilot DRAM simulation\n');
fprintf('====================================\n');

pilotSeed = 20262600;

rng(pilotSeed, 'twister');

pilotNsimu = 50000;
pilotBurnin = 10000;

pilotModel = model;

pilotModel.sigma2 = max(...
    baselineSSE / model.N, eps);

pilotOptions = mcmcset(...
    'nsimu', pilotNsimu, ...
    'burnintime', pilotBurnin, ...
    'adaptint', 100, ...
    'drscale', 2, ...
    'method', 'dram', ...
    'updatesigma', 1, ...
    'verbosity', 1);

[pilotResults, pilotChain, pilotSigma2, pilotSSE] = ...
    mcmcrun(pilotModel, data, params, pilotOptions);

%% Validate pilot chain
assert(size(pilotChain, 1) == pilotNsimu, ...
    'Unexpected pilot chain length.');

assert(size(pilotChain, 2) == npar, ...
    'Unexpected pilot parameter count.');

assert(all(isfinite(pilotChain(:))), ...
    'Pilot chain contains nonfinite values.');

assert(all(isfinite(pilotSigma2(:))) && ...
       all(pilotSigma2(:) > 0), ...
       'Invalid pilot residual variance.');

assert(all(isfinite(pilotSSE(:))) && ...
       all(pilotSSE(:) >= 0) && ...
       all(pilotSSE(:) < 1e12), ...
       'Invalid pilot SSE.');

%% Save pilot results
pilotResults.names = parameter_names;
pilotResults.class = 'MCMC';

save('PCC_Step6_Pilot.mat', ...
    'pilotResults', ...
    'pilotChain', ...
    'pilotSigma2', ...
    'pilotSSE', ...
    'pilotModel', ...
    'pilotOptions', ...
    'pilotSeed', ...
    'params', ...
    'parameter_names', ...
    'data', ...
    'Y0');

pilotTable = array2table(pilotChain, ...
    'VariableNames', parameter_names);

pilotTable.sigma2 = pilotSigma2(:);
pilotTable.SSE = pilotSSE(:);

pilotTable.delta_beta = ...
    pilotTable.beta_D - pilotTable.beta_F;

pilotTable.iteration = (1:pilotNsimu)';

writetable(pilotTable, ...
    'PCC_Step6_Pilot_chain.csv');

%% Pilot sampling efficiency
pilotRetained = pilotChain(pilotBurnin+1:end, :);

pilotMovement = mean(any(...
    diff(pilotRetained, 1, 1) ~= 0, 2));

fprintf('Pilot post-burn-in movement: %.4f%%\n', ...
    100 * pilotMovement);

if pilotMovement < 0.01
    warning(['Pilot movement below 1%%. ', ...
             'Proposal calibration may be unreliable.']);
end

%% Calibrate proposal covariance
% Pilot samples are used only for proposal calibration.
% They are excluded from final posterior inference.

pilotSegment = pilotChain(30001:50000, :);

pilotCov = cov(pilotSegment);

minimumSD = [
    0.0001
    0.0001
    0.0010
    0.0001
    0.0001
    0.0001
    0.0010
    0.0010
    0.0010
    0.0010
    0.0010
];

%% Regularize proposal covariance
pilotCov = 0.95 * pilotCov ...
         + 0.05 * diag(diag(pilotCov));

pilotCov = pilotCov + diag(minimumSD.^2);

%% Scale proposal covariance
qcov = (2.38^2 / npar) * pilotCov;

qcov = (qcov + qcov') / 2;

[~, covarianceFlag] = chol(qcov);

assert(covarianceFlag == 0, ...
    'Proposal covariance is not positive definite.');

disp('Pilot-calibrated proposal covariance constructed.');

%% Extended DRAM simulation settings
nsimu = 200000;
burnin = 20000;
nchains = 4;

options = mcmcset(...
    'nsimu', nsimu, ...
    'burnintime', burnin, ...
    'adaptint', 100, ...
    'drscale', 2, ...
    'method', 'dram', ...
    'updatesigma', 1, ...
    'qcov', qcov, ...
    'verbosity', 1);

%% Independent chain starting values
startIndices = [15000, 23000, 35000, 48000];

startingValues = pilotChain(startIndices, :);

if size(unique(startingValues, 'rows'), 1) < nchains
    warning('Some pilot-derived starting values are identical.');
end

%% Independent random seeds
seeds = [
    20262601
    20262602
    20262603
    20262604
];

%% Allocate chain storage
allChains = zeros(nsimu, npar, nchains);

allSigma2 = zeros(nsimu, nchains);
allSSE = zeros(nsimu, nchains);

allResults = cell(nchains, 1);

movementRates = zeros(nchains, 1);
acceptanceRates = zeros(nchains, 1);

meanBetaD = zeros(nchains, 1);
meanBetaF = zeros(nchains, 1);
meanDeltaBeta = zeros(nchains, 1);
probBetaDGreater = zeros(nchains, 1);

meanSSE = zeros(nchains, 1);
meanSigma2 = zeros(nchains, 1);

%% Run four independent DRAM chains
for c = 1:nchains

    fprintf('\n====================================\n');
    fprintf('Progressive Care Unit: Extended DRAM Chain %d of %d\n', ...
        c, nchains);
    fprintf('====================================\n');

    rng(seeds(c), 'twister');

    chainParams = params;

    %% Assign chain-specific starting values
    for j = 1:npar

        margin = 1e-8 * ...
            (upper_bounds(j) - lower_bounds(j));

        initialValue = startingValues(c, j);

        initialValue = max(initialValue, ...
            lower_bounds(j) + margin);

        initialValue = min(initialValue, ...
            upper_bounds(j) - margin);

        chainParams{j}{2} = initialValue;

    end

    %% Check starting SSE
    thetaStart = cellfun(@(p) p{2}, chainParams);

    initialSSE = ssfun_mrsa(thetaStart, data, Y0);

    assert(isfinite(initialSSE) && ...
           initialSSE >= 0 && ...
           initialSSE < 1e12, ...
           'Invalid initial SSE.');

    fprintf('Initial SSE: %.12g\n', initialSSE);

    %% Initialize residual variance
    chainModel = model;

    chainModel.sigma2 = max(...
        initialSSE / model.N, eps);

    % Preserve original residual variance prior
    chainModel.S20 = 1;
    chainModel.N0 = 1;

    fprintf('Initial residual variance: %.12g\n', ...
        chainModel.sigma2);

    %% Run DRAM
    [results, chain, s2chain, sschain] = ...
        mcmcrun(chainModel, data, chainParams, options);

    results.names = parameter_names;
    results.class = 'MCMC';

    %% Validate chain dimensions
    assert(size(chain, 1) == nsimu, ...
        'Unexpected chain length.');

    assert(size(chain, 2) == npar, ...
        'Unexpected number of parameters.');

    assert(numel(s2chain) == nsimu, ...
        'Unexpected residual variance chain length.');

    assert(numel(sschain) == nsimu, ...
        'Unexpected SSE chain length.');

    %% Validate numerical values
    assert(all(isfinite(chain(:))), ...
        'Nonfinite parameter samples.');

    assert(all(isfinite(s2chain(:))) && ...
           all(s2chain(:) > 0), ...
           'Invalid residual variance samples.');

    assert(all(isfinite(sschain(:))) && ...
           all(sschain(:) >= 0) && ...
           all(sschain(:) < 1e12), ...
           'Invalid SSE samples.');

    %% Verify parameter bounds
    assert(all(all(bsxfun(@ge, chain, ...
        lower_bounds(:)'))), ...
        'Lower parameter bound violation.');

    assert(all(all(bsxfun(@le, chain, ...
        upper_bounds(:)'))), ...
        'Upper parameter bound violation.');

    %% Store results
    allChains(:, :, c) = chain;

    allSigma2(:, c) = s2chain(:);
    allSSE(:, c) = sschain(:);

    allResults{c} = results;

    %% Remove burn-in for preliminary summaries
    retained = chain(burnin+1:end, :);

    retainedSigma2 = s2chain(burnin+1:end);
    retainedSSE = sschain(burnin+1:end);

    %% Calculate chain movement
    movementRate = mean(any(...
        diff(retained, 1, 1) ~= 0, 2));

    acceptanceRate = 1 - results.rejected;

    movementRates(c) = movementRate;
    acceptanceRates(c) = acceptanceRate;

    fprintf('Post-burn-in movement rate: %.4f%%\n', ...
        100 * movementRate);

    fprintf('Overall acceptance rate: %.4f%%\n', ...
        100 * acceptanceRate);

    %% Preliminary transmission comparison
    beta_D = retained(:,1);
    beta_F = retained(:,2);

    delta_beta = beta_D - beta_F;

    meanBetaD(c) = mean(beta_D);
    meanBetaF(c) = mean(beta_F);
    meanDeltaBeta(c) = mean(delta_beta);

    probBetaDGreater(c) = mean(delta_beta > 0);

    meanSSE(c) = mean(retainedSSE);
    meanSigma2(c) = mean(retainedSigma2);

    fprintf('Mean beta_D: %.8f\n', meanBetaD(c));
    fprintf('Mean beta_F: %.8f\n', meanBetaF(c));

    fprintf('Mean transmission difference: %.8f\n', ...
        meanDeltaBeta(c));

    fprintf('P(beta_D > beta_F): %.6f\n', ...
        probBetaDGreater(c));

    fprintf('Mean retained SSE: %.8f\n', ...
        meanSSE(c));

    fprintf('Mean retained residual variance: %.8f\n', ...
        meanSigma2(c));

    %% Export complete individual chain
    chainTable = array2table(chain, ...
        'VariableNames', parameter_names);

    chainTable.sigma2 = s2chain(:);
    chainTable.SSE = sschain(:);

    chainTable.delta_beta = ...
        chainTable.beta_D - chainTable.beta_F;

    chainTable.iteration = (1:nsimu)';
    chainTable.chain_id = repmat(c, nsimu, 1);

    outputCSV = sprintf(...
        'PCC_Step6_Extended_Chain_%d.csv', c);

    writetable(chainTable, outputCSV);

    %% Save individual MATLAB results
    outputMAT = sprintf(...
        'PCC_Step6_Extended_Chain_%d.mat', c);

    save(outputMAT, ...
        'results', ...
        'chain', ...
        's2chain', ...
        'sschain', ...
        'chainParams', ...
        'chainModel', ...
        'options', ...
        'movementRate', ...
        'acceptanceRate', ...
        'parameter_names', ...
        'burnin');

    fprintf('Chain %d complete and saved.\n', c);

end

%% Verify chains are not identical
for c1 = 1:nchains-1

    for c2 = c1+1:nchains

        assert(~isequal(allChains(:,:,c1), ...
                        allChains(:,:,c2)), ...
               'Duplicate parameter chains detected.');

    end

end

%% Save combined MATLAB results
save('PCC_Step6_Extended_All_Chains.mat', ...
    'allChains', ...
    'allSigma2', ...
    'allSSE', ...
    'allResults', ...
    'parameter_names', ...
    'params', ...
    'startingValues', ...
    'qcov', ...
    'seeds', ...
    'pilotSeed', ...
    'pilotMovement', ...
    'data', ...
    'Y0', ...
    'model', ...
    'options', ...
    'burnin', ...
    'movementRates', ...
    'acceptanceRates');

%% Export chain summary
summaryTable = table(...
    (1:nchains)', ...
    100 * movementRates, ...
    100 * acceptanceRates, ...
    meanBetaD, ...
    meanBetaF, ...
    meanDeltaBeta, ...
    probBetaDGreater, ...
    meanSSE, ...
    meanSigma2, ...
    'VariableNames', ...
    {'Chain_ID', ...
     'PostBurnin_Movement_Percent', ...
     'Overall_Acceptance_Percent', ...
     'Mean_beta_D', ...
     'Mean_beta_F', ...
     'Mean_delta_beta', ...
     'Probability_beta_D_Greater', ...
     'Mean_SSE', ...
     'Mean_sigma2'});

writetable(summaryTable, ...
    'PCC_Step6_Extended_Chain_Summary.csv');

disp(summaryTable);

%% Compare transmission trace plots
figure('Color', 'w');

subplot(3,1,1);
hold on;

for c = 1:nchains
    plot(allChains(:,1,c));
end

xline(burnin, '--k');

title('Progressive Care Unit: beta_D Independent Chains');
xlabel('Iteration');
ylabel('beta_D');

legend('Chain 1','Chain 2','Chain 3','Chain 4');

hold off;

subplot(3,1,2);
hold on;

for c = 1:nchains
    plot(allChains(:,2,c));
end

xline(burnin, '--k');

title('Progressive Care Unit: beta_F Independent Chains');
xlabel('Iteration');
ylabel('beta_F');

hold off;

subplot(3,1,3);
hold on;

for c = 1:nchains
    plot(allChains(:,1,c) - allChains(:,2,c));
end

xline(burnin, '--k');
yline(0, ':k');

title('Progressive Care Unit: Transmission Difference');
xlabel('Iteration');
ylabel('beta_D - beta_F');

hold off;

saveas(gcf, ...
    'PCC_Step6_Extended_Trace_Comparison.png');

%% Compare SSE across chains
figure('Color', 'w');

hold on;

for c = 1:nchains
    plot(allSSE(:,c));
end

xline(burnin, '--k');

title('Progressive Care Unit: Sum of Squared Errors');
xlabel('Iteration');
ylabel('SSE');

legend('Chain 1','Chain 2','Chain 3','Chain 4');

hold off;

saveas(gcf, ...
    'PCC_Step6_Extended_SSE_Comparison.png');

%% Compare residual variance across chains
figure('Color', 'w');

hold on;

for c = 1:nchains
    plot(allSigma2(:,c));
end

xline(burnin, '--k');

title('Progressive Care Unit: Residual Variance');
xlabel('Iteration');
ylabel('sigma^2');

legend('Chain 1','Chain 2','Chain 3','Chain 4');

hold off;

saveas(gcf, ...
    'PCC_Step6_Extended_Variance_Comparison.png');

%% Completion message
fprintf('\n====================================\n');
fprintf('Progressive Care Unit Extended DRAM Analysis Complete\n');

fprintf('Pilot iterations: %d\n', pilotNsimu);
fprintf('Independent chains: %d\n', nchains);
fprintf('Iterations per chain: %d\n', nsimu);
fprintf('Burn-in per chain: %d\n', burnin);

fprintf('Total retained samples: %d\n', ...
    nchains * (nsimu - burnin));

fprintf('Mean chain movement: %.4f%%\n', ...
    100 * mean(movementRates));

fprintf('====================================\n');

disp('Convergence diagnostics must be evaluated separately.');