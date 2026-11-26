%% IU-2 (3 Blue): posterior credible intervals
% Run after IU_2_3_Blue_Final.m and independent convergence verification.
% Uses complete main-chain CSVs, not pilot samples.
clearvars;
nchains = 4;
nsimu = 200000;
burnin = 20000;
names = {'beta_D','beta_F','gamma','mu_D','mu_F','delta','sigma', ...
    'Lambda_S_D','Lambda_I_D','Lambda_S_F','Lambda_I_F','sigma2','SSE','delta_beta'};
upper = [0.30 0.30 1 1 1 1 1 100 10 100 10];
chains = cell(nchains,1);
for c = 1:nchains
    filename = sprintf('IU_2_3_Blue_Step6_Extended_Chain_%d.csv',c);
    T = readtable(filename);
    assert(height(T)==nsimu,'Unexpected chain length.');
    assert(all(ismember([names {'iteration','chain_id'}],T.Properties.VariableNames)), ...
        'Missing chain columns.');
    assert(isequal(T.iteration,(1:nsimu)') && all(T.chain_id==c), ...
        'Incorrect iteration or chain identifiers.');
    X = T{burnin+1:end,names};
    assert(all(isfinite(X(:))),'Nonfinite retained samples.');
    assert(all(all(X(:,1:11)>=0 & X(:,1:11)<=upper)), 'Parameter bound violation.');
    assert(all(X(:,12)>0) && all(X(:,13)>=0 & X(:,13)<1e12), 'Invalid variance or SSE.');
    assert(all(abs(X(:,14)-(X(:,1)-X(:,2)))<1e-10),'Incorrect transmission difference.');
    for j = 1:c-1
        assert(~isequal(X(:,1:11),chains{j}(:,1:11)),'Duplicate chains.');
    end
    chains{c} = X;
end
pooled = vertcat(chains{:});
n = size(pooled,1);
assert(n==720000,'Unexpected retained sample count.');
ci = quantile(pooled,[0.025 0.975],1);
summary = table(string(names(:)),mean(pooled,1)',std(pooled,0,1)', ...
    median(pooled,1)',ci(1,:)',ci(2,:)',repmat(n,numel(names),1), ...
    'VariableNames',{'Parameter','Mean','StdDev','Median','CI_2_5','CI_97_5','Samples'});
d = pooled(:,14);
comparison = table(mean(pooled(:,1)),mean(pooled(:,2)),mean(d),std(d), ...
    ci(1,14),ci(2,14),mean(d>0),mean(d<0),mean(d==0),n, ...
    'VariableNames',{'Beta_D','Beta_F','Difference','Difference_SD','CI_2_5', ...
    'CI_97_5','P_SU_Higher','P_NS_Higher','P_Equal','Samples'});
writetable(summary,'IU_2_3_Blue_CredibleIntervals.csv');
writetable(comparison,'IU_2_3_Blue_TransmissionComparison.csv');
disp(summary);
disp(comparison);
