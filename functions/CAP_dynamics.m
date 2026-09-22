function [fraction, persistence, counts, ...
          trans_mats, out_degree, in_degree] = ...
          CAP_dynamics(cap_all, k_best)
%CAP_dynamics Calculate CAP temporal dynamics for each participant.
%
%   Inputs:
%       cap_all  - [n_subjects × nTR] CAP sequence for each participant
%       k_best   - number of CAPs
%
%   Outputs:
%       fraction    - Fraction of TRs occupied by each CAP
%       persistence  - Average persistence of each CAP in TRs
%       counts   - Number of occurrences of each CAP
%       trans_mats   - Transition probability matrix for each participant
%       out_degree   - Number of non-self outgoing transitions
%       in_degree   - Number of non-self incoming transitions

    n_subjects = size(cap_all, 1);
    nTR = size(cap_all, 2);

    fraction    = zeros(n_subjects, k_best);
    persistence = zeros(n_subjects, k_best);
    counts      = zeros(n_subjects, k_best);

    trans_mats  = zeros(k_best, k_best, n_subjects);
    out_degree  = zeros(n_subjects, k_best);
    in_degree   = zeros(n_subjects, k_best);

    % Calculate CAP dynamics for each participant
    for i = 1:n_subjects
        
        sub_seq = cap_all(i, :);

        % Fraction, persistence, and counts

        for k = 1:k_best

            logical_vec = (sub_seq == k);

            % Fraction of TRs occupied by the current CAP
            
            fraction(i, k) = sum(logical_vec) / nTR;

            % Number of occurrences
            
            changes = diff([0, logical_vec, 0]);
            counts(i, k) = sum(changes == 1);

            % Average persistence in TRs
            
            if counts(i, k) > 0
                persistence(i, k) = ...
                    sum(logical_vec) / counts(i, k);
            else
                persistence(i, k) = 0;
            end

        end

        % Transition matrix

        tm = zeros(k_best, k_best);

        for t = 1:nTR-1

            curr = sub_seq(t);
            next = sub_seq(t + 1);

            tm(curr, next) = tm(curr, next) + 1;

        end

        % In- and out-transition counts

        tm_no_diag = tm;

        tm_no_diag(logical(eye(k_best))) = 0;

        out_degree(i, :) = sum(tm_no_diag, 2)';
        in_degree(i, :)  = sum(tm_no_diag, 1);

        % Normalize transition matrix

        for k = 1:k_best

            row_total = sum(sub_seq(1:end-1) == k);

            if row_total > 0
                tm(k, :) = tm(k, :) / row_total;
            end

        end

        trans_mats(:, :, i) = tm;

    end

end