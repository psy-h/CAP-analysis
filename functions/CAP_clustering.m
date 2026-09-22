function results = CAP_clustering(all_data, k_range)
%CAP_CLUSTERING Perform CAP clustering across a range of cluster numbers.
%
%   results = CAP_clustering(all_data, k_range)
%
%   Inputs:
%       all_data - [n_samples × n_ROIs] concatenated fMRI data
%       k_range  - Vector specifying the numbers of CAPs to evaluate
%
%   Outputs:
%       results  - Structure containing clustering results and evaluation
%                  metrics for each value of k
results = struct();

for i = 1:numel(k_range)

    k = k_range(i);

    fprintf('Running k = %d...\n', k);

    [idx, ~, sumd] = kmeans(all_data, k, ...
        'Distance', 'correlation', ...
        'Replicates', 100, ...
        'MaxIter', 500, ...
        'Display', 'final');

    mean_s = mean( ...
        silhouette(all_data, idx, 'correlation'));


    results(i).k = k;
    results(i).idx = idx;
    results(i).sumd = sumd;
    results(i).WCSS = sum(sumd);
    results(i).mean_silhouette = mean_s;

end

end
