%% CAP clustering based on healthy control participants

clear;
clc;

% Load data and group information

load('fmri_data.mat');
info = readtable('participant_info.xlsx');

% Check that the participant information matches the fMRI data
n_total = numel(fmri_data);

if height(info) ~= n_total
    error('The number of participants in participant_info.xlsx does not match fmri_data.');
end

% Group information
group = info.Group;

idx_healthy = find(group == 0);
idx_patient = find(group == 1);

n_H = numel(idx_healthy);
n_P = numel(idx_patient);

nTR = size(fmri_data{1}, 1);
nROI = size(fmri_data{1}, 2);

for i = 1:n_total
    if size(fmri_data{i}, 1) ~= nTR || size(fmri_data{i}, 2) ~= nROI
        error('All participants must have the same number of TRs and ROIs.');
    end
end

% Prepare concatenated healthy-control data

total_TR = n_H * nTR;
all_data = zeros(total_TR, nROI);

row_start = 1;

for i = 1:n_H

    ts = fmri_data{idx_healthy(i)};

    % Z-score each ROI within subject
    ts_z = zscore(ts, 0, 1);

    row_end = row_start + nTR - 1;

    all_data(row_start:row_end, :) = ts_z;

    row_start = row_end + 1;

end

fprintf('Total concatenated TRs: %d\n', size(all_data, 1));

% CAP clustering

k_range = 2:11;

results = CAP_clustering(all_data, k_range);

%% Plot clustering evaluation metrics

CAP_plot(results, k_range);

%% Select the number of CAPs

k_best = 4;
k_best_index = find(k_range == k_best, 1);

%% Calculate CAP centroids and Z-maps

[centroids, z_centroids] = CAP_centroids( ...
    all_data, ...
    results(k_best_index).idx, ...
    k_best);

%% Assign CAPs to all participants

cap_all = zeros(n_total, nTR);

% Healthy controls

labels_hc = results(k_best_index).idx;

for i = 1:n_H

    orig_idx = idx_healthy(i);

    range_hc = (i - 1) * nTR + 1 : i * nTR;

    cap_all(orig_idx, :) = labels_hc(range_hc)';

end


% Patients

for i = 1:n_P

    orig_idx = idx_patient(i);

    ts = fmri_data{orig_idx};

    ts_z = zscore(ts, 0, 1);

    [~, labels_patient] = pdist2( ...
        centroids, ...
        ts_z, ...
        'correlation', ...
        'Smallest', 1);

    cap_all(orig_idx, :) = labels_patient';

end

%% Calculate CAP dynamics
[fraction_all, persistence_all, counts_all, ...
    trans_mats_all, out_degree_all, in_degree_all] = ...
    CAP_dynamics(cap_all, k_best);

%% Group comparison

use_covariates = true;
covariate_names = {};   % Empty = use all available covariates

% Fraction
fraction_results = compare_CAP( ...
    fraction_all, info, use_covariates, covariate_names);

% Persistence
persistence_results = compare_CAP( ...
    persistence_all, info, use_covariates, covariate_names);

% Counts
counts_results = compare_CAP( ...
    counts_all, info, use_covariates, covariate_names);

% Transition probabilities
trans_data = zeros(n_total, k_best * k_best);

for i = 1:n_total
    tm = trans_mats_all(:, :, i);

    trans_data(i, :) = reshape(tm', 1, []);
end

transition_results = compare_CAP( ...
    trans_data, info, use_covariates, covariate_names);

%% FDR correction

all_p = [
    fraction_results.P_Value;
    persistence_results.P_Value;
    counts_results.P_Value;
    transition_results.P_Value
];

adjusted_p = mafdr(all_p, 'BHFDR', true);

n_fraction = height(fraction_results);
n_persistence = height(persistence_results);
n_counts = height(counts_results);

fraction_results.FDR_P_Value = ...
    adjusted_p(1:n_fraction);

persistence_results.FDR_P_Value = ...
    adjusted_p(n_fraction + 1 : ...
               n_fraction + n_persistence);

counts_results.FDR_P_Value = ...
    adjusted_p(n_fraction + n_persistence + 1 : ...
               n_fraction + n_persistence + n_counts);

transition_results.FDR_P_Value = ...
    adjusted_p(n_fraction + n_persistence + n_counts + 1 : end);

%% Save analysis results

save('CAP_analysis_results.mat', ...
    'fraction_results', ...
    'persistence_results', ...
    'counts_results', ...
    'transition_results', ...
    'centroids', ...
    'z_centroids', ...
    'cap_all', ...
    'fraction_all', ...
    'persistence_all', ...
    'counts_all', ...
    'trans_mats_all', ...
    'k_best', ...
    'k_range');