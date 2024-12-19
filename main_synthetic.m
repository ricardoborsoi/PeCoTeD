% ========================================================================= 
% This code implements an example with synthetic data to compare the methods 
% presented in the paper:
%
%   Personalized Coupled Tensor Decomposition for Multimodal Data Fusion: Uniqueness and Algorithms
%   R.A. Borsoi, K. Usevich, D. Brie, T. Adali
%   IEEE Transactions on Signal Processing, 2024.
%
% Ricardo Borsoi
% =========================================================================

clear all
clc

input('WARNING: to run the codes, you need to add the freeLYAP, Tensorlab, and YALMIP packages \n in the "assets" folder! (press enter to continue...)')

addpath(genpath('assets'))
addpath(genpath('algorithm'))
% addpath(genpath('freeLYAP'))
% addpath(genpath('tensorlab_2016-03-28'))
% addpath(genpath('YALMIP'))

K       = 3; % number of observed tensors
Gamma = {[1,2,3], ... k=1
         [1,2,3], ... k=2
         [1,2,3], ... k=3
         }; % couplings

% dimensions of Z
M1 = 7;
M2 = 11;
M3 = 9;

% one mode per tensor has full rank Pj{k}
%    k=1   k=2   k=3   
N1 = {M1+3, 5,    5};
N2 = {5,    M2+1, 7};
N3 = {7,    7,    M3+1};

% ranks
rankZ   = 5;
rankPsi = 5; % use the same for all modes, for simplicity
disp('only tensor k=2 is fully identifiable, and each k is mode-k identifiable')


P1 = cell(K,1); P2 = cell(K,1); P3 = cell(K,1);
% set degradation in some modes
for k=1:K, P1{k} = rand(N1{k}, M1); end
for k=1:K, P2{k} = rand(N2{k}, M2); end
for k=1:K, P3{k} = rand(N3{k}, M3); end



num_MC_runs = 20;
SNRs_all = 20:10:70;


% average errors and times
mc_err_z_cp_opt_a = zeros(num_MC_runs, length(SNRs_all));
mc_err_z_cp_opt_r = zeros(num_MC_runs, length(SNRs_all));
mc_err_z_cp_alg   = zeros(num_MC_runs, length(SNRs_all));
mc_err_z_cp_SDF_r = zeros(num_MC_runs, length(SNRs_all));

mc_time_cp_opt_a = zeros(num_MC_runs, length(SNRs_all));
mc_time_cp_opt_r = zeros(num_MC_runs, length(SNRs_all));
mc_time_cp_alg   = zeros(num_MC_runs, length(SNRs_all));
mc_time_cp_SDF_r = zeros(num_MC_runs, length(SNRs_all));

rng("shuffle")




for snr_count=1:length(SNRs_all)
SNR = SNRs_all(snr_count);
for n_run=1:num_MC_runs

    fprintf('MC run: %d, SNR: %d \n', n_run, SNR)

    % generate Z
    Z = cpdgen({randn(M1,rankZ), randn(M2,rankZ), randn(M3,rankZ)});
    
    % generate Psis and Y
    Y = cell(K,1);
    for k=1:K
        Y{k} = tmprod(Z, {P1{k}, P2{k}, P3{k}}, [1,2,3]) ...
            + cpdgen({randn(size(P1{k},1),rankPsi), randn(size(P2{k},1),rankPsi), randn(size(P3{k},1),rankPsi)});
        % add noise
        Y{k} = Y{k} + sqrt(sum((Y{k}(:)).^2)/numel(Y{k})/10^(SNR/10)) * randn(size(Y{k}));
    end
    
    
    
    % reconstruct Z
    opts.num_iters = 1000;
    % opts.initOpt = 'random';
    % opts.initOpt = 'semialgebraic';
    % opts.initOpt = 'cpd';
    % opts.initialCPD_initoption = 'auto'; % 'auto', 'random', 'gesvd' % initialization option for CPD inside the codes
    opts.alg_use_regression = false; % use regression to estimate the factors in the algebraic algorithm (can be useful when some Pjk are identity)


    % runs different initializations and select the solution with the best data fit
    num_runs = 50;
    
    tic
    opts.initOpt = 'semialgebraic';
    flag_algebraic_or_optim = 'optimization';
    [Z_cp_opt_a, Psi_cp_opt_a, costf_cp_opt_a] = PECOTED_multi_init_runs(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts, num_runs, flag_algebraic_or_optim);
    t_opt_a = toc / num_runs;

    tic
    flag_algebraic_or_optim = 'algebraic';
    [Z_cp_alg, Psi_cp_alg, costf_cp_alg] = PECOTED_multi_init_runs(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts, num_runs, flag_algebraic_or_optim);
    t_alg = toc / num_runs;
    
    tic
    opts.initOpt = 'random';
    flag_algebraic_or_optim = 'optimization';
    [Z_cp_opt_r, Psi_cp_opt_r, costf_cp_opt_r] = PECOTED_multi_init_runs(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts, num_runs, flag_algebraic_or_optim);
    t_opt_r = toc / num_runs;

    % this one is only with random initialization
%     tic
%     [Z_SDF, Psi_SDF, costf_SDF] = fPECOTED_SDF_multi_init_runs(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts, num_runs);
%     t_SDF_r = toc / num_runs;


    % compute avg errors and times
    mc_err_z_cp_opt_a(n_run, snr_count) = comp_NRMSE(Z, Z_cp_opt_a);
    mc_err_z_cp_opt_r(n_run, snr_count) = comp_NRMSE(Z, Z_cp_opt_r);
    mc_err_z_cp_alg(n_run, snr_count)   = comp_NRMSE(Z, Z_cp_alg);
%     mc_err_z_cp_SDF_r(n_run, snr_count) = comp_NRMSE(Z, Z_SDF);

    

    % store times:
    mc_time_cp_opt_a(n_run, snr_count) = t_opt_a;
    mc_time_cp_opt_r(n_run, snr_count) = t_opt_r;
    mc_time_cp_alg(n_run, snr_count)   = t_alg;
%     mc_time_cp_SDF_r(n_run, snr_count) = t_SDF_r;
end
end

%%

% plot NRMSEs
figure, hold on
plot(SNRs_all, mean(mc_err_z_cp_alg,1), '-x')
plot(SNRs_all, mean(mc_err_z_cp_opt_a,1), '-*')
plot(SNRs_all, mean(mc_err_z_cp_opt_r,1), '-d')
% plot(SNRs_all, mean(mc_err_z_cp_SDF_r,1), '--')

set(gca, 'YScale', 'log')
legend({'Semi-algebraic','Optimization (init. 1)','Optimization (init. 2)', 'SDF'})
xlabel('SNR [dB]')
xlim([20, 70])
% ylim([0.0010, 152])
ylabel('NRMSE')
grid


% show times 
fprintf('\n\n')
fprintf('Semi-Alg Time  %.3f, std= %1.3f \n', mean(mc_time_cp_alg(:)),   std(mc_time_cp_alg(:)) )
fprintf('Opt in.1 Time  %.3f, std= %1.3f \n', mean(mc_time_cp_opt_a(:)), std(mc_time_cp_opt_a(:)) )
fprintf('Opt in.2 Time  %.3f, std= %1.3f \n', mean(mc_time_cp_opt_r(:)), std(mc_time_cp_opt_r(:)) )
% fprintf('SDF Time       %.3f, std= %1.3f \n', mean(mc_time_cp_SDF_r(:)), std(mc_time_cp_SDF_r(:)) )



%%
% =========================================================================
function [nrmsez]=comp_NRMSE(Ztrue, Zhat)
    nrmsez = norm(Zhat(:)-Ztrue(:))/norm(Ztrue(:));
end
