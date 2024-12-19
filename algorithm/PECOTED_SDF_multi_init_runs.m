function [Z_, Psi_, costf_] = PECOTED_SDF_multi_init_runs(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts, num_runs)
% =========================================================================
% runs the SDF framework of Sorber et al. algorihtm several
% times with different random initializations and 
% picks the one which gives the best reconstruction error
%
% Thhis code implements algorithms related to the paper:
%   Personalized Coupled Tensor Decomposition for Multimodal Data Fusion: Uniqueness and Algorithms
%   R.A. Borsoi, K. Usevich, D. Brie, T. Adali
%   IEEE Transactions on Signal Processing, 2024.
%
% Ricardo Borsoi
% =========================================================================



% compute first run
[Z_, Psi_, C_, D_, costf_] = solve_PECOTED_with_SDF(Y, rankZ, rankPsi, P1, P2, P3);

% compute remaining runs to see if we find a better one
for i=1:(num_runs-1)
    [Z, Psi, C, D, costf] = solve_PECOTED_with_SDF(Y, rankZ, rankPsi, P1, P2, P3);
    % disp(costf)

    if costf < costf_
        Z_ = Z;
        Psi_ = Psi;
        costf_ = costf;
    end
end

