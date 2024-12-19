function [Z_, Psi_, costf_] = PECOTED_multi_init_runs(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts, num_runs, flag_algebraic_or_optim)
% =========================================================================
% runs the algorihtm several times with different random initializations 
% and picks the one which gives the best reconstruction error
%
% This code implements algorithms related to the paper:
%   Personalized Coupled Tensor Decomposition for Multimodal Data Fusion: Uniqueness and Algorithms
%   R.A. Borsoi, K. Usevich, D. Brie, T. Adali
%   IEEE Transactions on Signal Processing, 2024.
%
% Ricardo Borsoi
% =========================================================================

% choose which of the algorihtms to run
if nargin < 10
    flag_algebraic_or_optim = 'optimization'; % 'optimization' or 'algebraic'
end


if strcmp(flag_algebraic_or_optim, 'optimization')
    % opts.initialCPD_initoption = 'gesvd'; % 'auto', 'random', 'gesvd' % initialization option for CPD inside the codes
    opts.initialCPD_initoption = 'auto';
    [Z_, Psi_, Bz1, Bz2, Bz3, B1, B2, B3, costf_, C1, C2, C3] = solve_PECOTED_optimization(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts);
    % disp(costf_)

    opts.initialCPD_initoption = 'random';
    for i=1:(num_runs-1)
        [Z, Psi, Bz1, Bz2, Bz3, B1, B2, B3, costf, C1, C2, C3] = solve_PECOTED_optimization(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts);
        % disp(costf)

        if costf < costf_
            Z_ = Z;
            Psi_ = Psi;
            costf_ = costf;
        end
    end

elseif strcmp(flag_algebraic_or_optim, 'algebraic')
    % opts.initialCPD_initoption = 'gesvd'; % 'auto', 'random', 'gesvd' % initialization option for CPD inside the codes
    opts.initialCPD_initoption = 'auto';
    [Z_, Psi_, C, costf_] = solve_PECOTED_semialgebraic(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts);
    % disp(costf_)

    opts.initialCPD_initoption = 'random';
    for i=1:(num_runs-1)
        [Z, Psi, C, costf] = solve_PECOTED_semialgebraic(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts);
        % disp(costf)
        
        if costf < costf_
            Z_ = Z;
            Psi_ = Psi;
            costf_ = costf;
        end
    end

else
    error('unknown algorithm option, must be "optimization" or "algebraic"')
end







