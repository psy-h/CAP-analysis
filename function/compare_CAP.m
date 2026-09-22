function results = compare_CAP( ...
    data, info, use_covariates, covariate_names)
%COMPARE_CAP Compare CAP metrics between two groups.
%
%   results = compare_CAP(data, info, use_covariates, covariate_names)
%
%   Inputs:
%       data             - [n_subjects x n_metrics] data matrix
%       info             - participant information table
%       use_covariates   - true/false, whether to adjust for covariates
%       covariate_names  - cell array of covariate names
%                          {} means all variables except Group
%
%   Outputs:
%       results          - table containing group means, statistic,
%                          and p-values
%
%   Notes:
%       Group must be coded as 0 = HC and 1 = Patient.
%       The row order of info must match the participant order in data.

    % Check number of participants
    if size(data, 1) ~= height(info)
        error(['The number of participants in data does not match ' ...
               'the number of rows in info.']);
    end

    % Check Group variable
    if ~ismember('Group', info.Properties.VariableNames)
        error('The participant information table must contain "Group".');
    end

    group = info.Group;

    % Check that Group is coded as 0 = HC and 1 = Patient
    if ~all(ismember(unique(group), [0 1]))
        error('Group must be coded as 0 = HC and 1 = Patient.');
    end

    % Automatically use all available covariates
    if use_covariates && isempty(covariate_names)

        all_variables = info.Properties.VariableNames;

        covariate_names = setdiff(all_variables, {'Group'});

    end

    % Check selected covariates
    if use_covariates

        for c = 1:numel(covariate_names)

            if ~ismember(covariate_names{c}, ...
                    info.Properties.VariableNames)

                error('Covariate "%s" was not found in info.', ...
                    covariate_names{c});

            end

        end

    end

    n_metrics = size(data, 2);

    % Initialize results
    results = table( ...
        zeros(n_metrics, 1), ...
        zeros(n_metrics, 1), ...
        zeros(n_metrics, 1), ...
        zeros(n_metrics, 1), ...
        'VariableNames', ...
        {'Mean_Patient', ...
         'SD_Patient', ...
         'Mean_HC', ...
         'SD_HC'});

    % Add statistic and p-value columns
    results.Statistic = zeros(n_metrics, 1);
    results.P_Value = zeros(n_metrics, 1);

    % Group comparison
    for m = 1:n_metrics

        y = data(:, m);

        % Group means and SDs
        hc = y(group == 0);
        patient = y(group == 1);

        results.Mean_HC(m) = mean(hc, 'omitnan');
        results.SD_HC(m) = std(hc, 'omitnan');

        results.Mean_Patient(m) = mean(patient, 'omitnan');
        results.SD_Patient(m) = std(patient, 'omitnan');

        % Build regression table
        T = table( ...
            y, ...
            categorical(group), ...
            'VariableNames', {'Metric', 'Group'});

        % Add covariates if requested
        if use_covariates

            for c = 1:numel(covariate_names)

                name = covariate_names{c};

                T.(name) = info.(name);

            end

        end

        % Build regression formula
        formula = 'Metric ~ Group';

        if use_covariates

            for c = 1:numel(covariate_names)

                formula = [ ...
                    formula ' + ' covariate_names{c}];

            end

        end

        % Fit linear model
        model = fitlm(T, formula);

        % Find group effect
        group_idx = strcmp( ...
            model.CoefficientNames, 'Group_1');

        if ~any(group_idx)
            error('Could not identify the group coefficient.');
        end

        % Extract t statistic and p-value
        results.Statistic(m) = ...
            model.Coefficients.tStat(group_idx);

        results.P_Value(m) = ...
            model.Coefficients.pValue(group_idx);

    end

end
