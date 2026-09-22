function CAP_plot(results, k_range)
%CAP_PLOT Plot clustering evaluation metrics.
%
%   CAP_plot(results, k_range)
%
%   Inputs:
%       results - Output structure from CAP_clustering
%       k_range - Number of CAPs evaluated
%
%   The function plots:
%       1. WCSS
%       2. silhouette coefficient

    % Extract clustering metrics
    wcss = [results.WCSS];
    silhouette_values = [results.mean_silhouette];

    % Check input consistency
    if numel(wcss) ~= numel(k_range)
        error('The number of clustering results does not match k_range.');
    end

    if numel(silhouette_values) ~= numel(k_range)
        error('The number of silhouette values does not match k_range.');
    end

    % WCSS curve

    fig1 = figure('Color', 'w');

    plot(k_range, wcss, '-o', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6);

    xlabel('Number of CAPs (k)');
    title('WCSS');

    xticks(k_range);
    grid on;
    box off;

    exportgraphics(fig1, 'WCSS_curve.tif', ...
    'Resolution', 600);

    % Silhouette coefficient curve

    fig2 = figure('Color', 'w');

    plot(k_range, silhouette_values, '-o', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6);

    xlabel('Number of CAPs (k)');
    title('Silhouette Coefficient');

    xticks(k_range);
    grid on;
    box off;

    exportgraphics(fig2, 'silhouette_curve.tif', ...
    'Resolution', 600);

end
