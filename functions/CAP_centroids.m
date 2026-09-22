function [centroids, z_centroids] = CAP_centroids(data, labels, k)
%CAP_centroids Calculate cluster centroids and z_centroids.
%
%   [centroids, z_centroids] = CAP_centroids(data, labels, k)
%
%   Inputs:
%       data   - [n_samples × n_ROIs] data matrix
%       labels - [n_samples × 1] cluster labels
%       k      - number of clusters
%
%   Outputs:
%       centroids   - [k × n_ROIs] cluster centroids
%       z_centroids - [k × n_ROIs] z_centroids

    n_ROI = size(data, 2);

    centroids = zeros(k, n_ROI);
    z_centroids = zeros(k, n_ROI);

    for i = 1:k

        cluster_data = data(labels == i, :);

        mu = mean(cluster_data, 1);

        sigma = std(cluster_data, 0, 1);

        centroids(i, :) = mu;

        z_centroids(i, :) = mu ./ sigma;

    end

end