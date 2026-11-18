function mcmc_diagnostics_plot(chain, results, saveDir)
% MCMC_DIAGNOSTICS_PLOT - Generate and save diagnostic plots for MCMC output
%
% Inputs:
%   chain   - [nSimu × nParams] MCMC samples
%   results - struct with 'names' field: {param1, param2, ...}
%   saveDir - optional: directory to save figures (default = current folder)

    if nargin < 2 || isempty(results) || ~isfield(results, 'names')
        paramNames = strcat("param", string(1:size(chain,2)));
    else
        paramNames = string(results.names);
    end

    if nargin < 3 || isempty(saveDir)
        saveDir = pwd;  % current directory
    end

    nParams = size(chain, 2);
    nSimu = size(chain, 1);

    %% === TRACE PLOTS ===
    figure('Name', 'Trace Plots', 'NumberTitle', 'off', 'Position', [100, 100, 800, 200 + 150 * nParams]);
    for i = 1:nParams
        subplot(nParams, 1, i);
        plot(chain(:, i), 'b');
        title(sprintf('Trace: %s', paramNames(i)), 'Interpreter', 'none');
        ylabel(paramNames(i));
        xlabel('Iteration');
    end
    saveas(gcf, fullfile(saveDir, 'trace_plots.png'));

    %% === HISTOGRAMS ===
    figure('Name', 'Posterior Histograms', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
    nRows = ceil(sqrt(nParams));
    nCols = ceil(nParams / nRows);
    for i = 1:nParams
        subplot(nRows, nCols, i);
        histogram(chain(:, i), 40, 'FaceColor', [0.2 0.4 0.6]);
        title(sprintf('Posterior: %s', paramNames(i)), 'Interpreter', 'none');
        xlabel(paramNames(i));
        ylabel('Frequency');
    end
    saveas(gcf, fullfile(saveDir, 'posterior_histograms.png'));

    %% === PAIRWISE SCATTER PLOTS ===
    figure('Name', 'Pairwise Scatter Plots', 'NumberTitle', 'off', 'Position', [100, 100, 800, 800]);
    for i = 1:nParams
        for j = 1:nParams
            subplot(nParams, nParams, (i-1)*nParams + j);
            if i == j
                histogram(chain(:, i), 20, 'EdgeColor', 'none', 'FaceColor', [0.6 0.6 0.6]);
            else
                plot(chain(:, j), chain(:, i), '.', 'MarkerSize', 2);
            end
            if i == nParams
                xlabel(paramNames(j), 'FontSize', 6);
            else
                xticks([]);
            end
            if j == 1
                ylabel(paramNames(i), 'FontSize', 6);
            else
                yticks([]);
            end
        end
    end
    saveas(gcf, fullfile(saveDir, 'pairwise_scatter_plots.png'));

    disp("✅ MCMC diagnostic plots saved to: " + saveDir);
end
