function [Z, Psi, C, costf, Bpsi] = solve_PECOTED_semialgebraic(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts)
% =========================================================================
% Code for personalized tensor-based coupled data fusion with variability
% Semi-algebraic version
%
% This code implements algorithms related to the paper:
%   Personalized Coupled Tensor Decomposition for Multimodal Data Fusion: Uniqueness and Algorithms
%   R.A. Borsoi, K. Usevich, D. Brie, T. Adali
%   IEEE Transactions on Signal Processing, 2024.
%
%
% Y  : cell array of order-3 tensors 
% P1 : cell array 
% Gamma : K-dim cell array whose k-th cell contains the set of indices of
%         which modes are coupled for the k-th image (deprecated!!!)
% opts : struct with lots of options (see examle file)
% 
% Ricardo Borsoi
% =========================================================================

if nargin < 8
    opts = struct();
end

% initializatio of the CPD algorithms used to initialize the method
if ~isfield(opts, 'initialCPD_initoption')
    initialCPD_initoption = 'auto'; % options: 'auto', 'random', 'gesvd'
else, initialCPD_initoption = opts.initialCPD_initoption; end

% whether to estimate the common factors using an all-algebraic approach or
% using a corase but more robust regression (treating the distinct part of
% the model as additive noise)
if ~isfield(opts, 'alg_use_regression')
    opts.alg_use_regression = false;
end

% make struct with CPD init options for tensorlab
% options.Algorithm = @cpd_nls; @cpd_als; %[@cpd_als|@cpd_minf|{@cpd_nls}|@cpd3_sd|@cpd3_sgsd|@cpd_rbs]
opts_cpd_tensorlab = struct();
if strcmp(initialCPD_initoption, 'auto')
    % pass nothing
elseif strcmp(initialCPD_initoption, 'random')
    opts_cpd_tensorlab.Initialization = @cpd_rnd;
elseif strcmp(initialCPD_initoption, 'gesvd')
    opts_cpd_tensorlab.Initialization = @cpd_gevd;
end



% get the number of measurements
K = length(Y);
R = rankZ;
if length(rankPsi) == 1, rankPsi = rankPsi*ones(K,1); end

if ~iscell(P1)
    P1_const = P1;
    P1 = cell(K,1);
    for k=1:K, P1{k} = P1_const; end
end
if ~iscell(P2)
    P2_const = P2;
    P2 = cell(K,1);
    for k=1:K, P2{k} = P2_const; end
end
if ~iscell(P3)
    P3_const = P3;
    P3 = cell(K,1);
    for k=1:K, P3{k} = P3_const; end
end

P_all = {P1; P2; P3};

% get latent tensor size
M1 = size(P1{1},2);
M2 = size(P2{1},2);
M3 = size(P3{1},2);


% find the unique tensors =================================================
eta_idx = []; % indices of unique tensors
for k=1:K
    % check if tensor k is unique
    if min(rank(P1{k},1e-5), R+rankPsi(k)) ...
       + min(rank(P2{k},1e-5), R+rankPsi(k)) ...
       + min(rank(P3{k},1e-5), R+rankPsi(k)) ...
       >= 2*(R+rankPsi(k))+2
        eta_idx = [eta_idx, k];
    end
end

% find the uni-mode unique tensors with full rank P
xi_idx = cell(3,1); 

xi_idx{1} = []; % indices of mode-1 unique tensors with full rank P1
for k=1:K
    % check if P1 is full column rank
    if rank(P1{k},1e-5) == size(P1{k},2)
        % check if mode-i identifiable
        if min(size(P1{k},1), min(size(P1{k},1),R)+rankPsi(k)) ...
           + min(rank(P2{k},1e-5), R+rankPsi(k)) ...
           + min(rank(P3{k},1e-5), R+rankPsi(k)) ...
           >= 2*(R+rankPsi(k))+2
            xi_idx{1} = [xi_idx{1}, k];
        end
    end
end

% find the uni-mode unique tensors with full rank P
xi_idx{2} = []; % indices of mode-1 unique tensors with full rank P2
for k=1:K
    % check if P2 is full column rank
    if rank(P2{k},1e-5) == size(P2{k},2)
        % check if mode-i identifiable
        if min(rank(P1{k},1e-5), R+rankPsi(k)) ...
           + min(size(P2{k},1), min(size(P2{k},1),R)+rankPsi(k)) ...
           + min(rank(P3{k},1e-5), R+rankPsi(k)) ...
           >= 2*(R+rankPsi(k))+2
            xi_idx{2} = [xi_idx{2}, k];
        end
    end
end

% find the uni-mode unique tensors with full rank P
xi_idx{3} = []; % indices of mode-1 unique tensors with full rank P3
for k=1:K
    % check if P1 is full column rank
    if rank(P3{k},1e-5) == size(P3{k},2)
        % check if mode-i identifiable
        if min(rank(P1{k},1e-5), R+rankPsi(k)) ...
           + min(rank(P2{k},1e-5), R+rankPsi(k)) ...
           + min(size(P3{k},1), min(size(P3{k},1),R)+rankPsi(k)) ...
           >= 2*(R+rankPsi(k))+2
            xi_idx{3} = [xi_idx{3}, k];
        end
    end
end


% find distinct (mode) unique tensors =====================================
% test if there is a \xi in uni-mode uniques different from \eta
j = 0; % defaults in case uni-mode unique tensor is different from etas
for i=1:length(eta_idx)
    if numel(setdiff(xi_idx{1}, eta_idx(i))) > 0  &&  ~isempty(xi_idx{1}) % true if xi_idx{1} (mode 1) has an index different from eta_idx(i)
        j = 1;
        eta = eta_idx(i);
        xi = setdiff(xi_idx{1}, eta); xi = xi(1);
        break;
    
    elseif numel(setdiff(xi_idx{2}, eta_idx(i))) > 0  &&  ~isempty(xi_idx{2})
        j = 2;
        eta = eta_idx(i);
        xi = setdiff(xi_idx{2}, eta); xi = xi(1);
        break;
    
    elseif numel(setdiff(xi_idx{3}, eta_idx(i))) > 0  &&  ~isempty(xi_idx{3})
        j = 3;
        eta = eta_idx(i);
        xi = setdiff(xi_idx{3}, eta); xi = xi(1);
        break;
    end
end
if j == 0 
    error(['no uni-mode unique and unique tensors are distinct, cannot solve ' ...
           'ambiguities to find the common part...'])
end


% compute CPDs for step 1 =================================================
% options.Initialization = @cpd_gevd; @cpd_rnd;
% options.Algorithm = @cpd_nls; @cpd_als; %[@cpd_als|@cpd_minf|{@cpd_nls}|@cpd3_sd|@cpd3_sgsd|@cpd_rbs]

[U_eta,output] = cpd(Y{eta}, R+rankPsi(eta), opts_cpd_tensorlab);
[U_xi,output]  = cpd(Y{xi},  R+rankPsi(xi),  opts_cpd_tensorlab);
V = P_all{j}{eta} * (P_all{j}{xi} \ U_xi{j});

% find common parts and permute ===========================================
[Delta,Delta_u,Delta_v] = find_partial_assignment(U_eta{j}, V, R);

Xc_eta = cell(3,1);
Xc_xi  = cell(3,1);
C      = cell(3,1);

idx_common_u = find(sum(Delta_u,2));
Xc_eta{1} = U_eta{1}(:,idx_common_u);
Xc_eta{2} = U_eta{2}(:,idx_common_u);
Xc_eta{3} = U_eta{3}(:,idx_common_u);

% Xc_xi{j}  = U_xi{j} * Delta_v;
Xc_xi{1}  = U_xi{1} * Delta_v;
Xc_xi{2}  = U_xi{2} * Delta_v;
Xc_xi{3}  = U_xi{3} * Delta_v;


% now choose how to estimate the factors

if opts.alg_use_regression == true
    % ---------------------------------------------------------------------
    % we have the common part of tensor eta, we can regress the factors
    % from the other uni-mode unique tensors for a coraser but more robust
    % estimate
    for jj=1:3
        Xc_xi{jj} = U_xi{jj} * Delta_v;
        C{jj} = P_all{jj}{xi} \ Xc_xi{jj}; % good or rough estimate of C{jj}
    end
    % C{j} = P_all{j}{xi} \ Xc_xi{j};
    
    % Take the next mode
    for i = setdiff([1,2,3], j)
        % test if mode-i of xi has full rank Pi{xi}
        if rank(P_all{i}{xi}, 1e-5) == size(P_all{i}{xi}, 2)
            % C{i} = P_all{i}{xi} \ Xc_xi{i};
            % in this case C{i} is already good from above
        else
            xi_new = xi_idx{i}(1); % index of the uni-mode unique tensor
            Y_xi_mat = tens2mat(Y{xi_new},[],i); % mode-i matricize Y{xi}
            if i==1
                % regress
                temp = khatri_rao(P_all{3}{xi_new} * C{3}, P_all{2}{xi_new} * C{2});
                C{1} = P_all{1}{xi_new} \ ((temp \ Y_xi_mat)');
            elseif i==2
                % regress
                temp = khatri_rao(P_all{3}{xi_new} * C{3}, P_all{1}{xi_new} * C{1});
                C{2} = P_all{2}{xi_new} \ ((temp \ Y_xi_mat)');
            elseif i==3
                % regress
                temp = khatri_rao(P_all{2}{xi_new} * C{2}, P_all{1}{xi_new} * C{1});
                C{3} = P_all{3}{xi_new} \ ((temp \ Y_xi_mat)');
            end
        end
    end

    

else
    % other method --------------------------------------------------------
    % find corresponding factors across modes
    % this is the method used for general P_{j,k} operators
    
    % Lambda = diag( sqrt(sum(Xc_eta{j}.^2,1)) ./ sqrt(sum((V*Delta_v).^2,1)) ) * diag(sign(?max?));
    Lambda = compute_scaling_compensator(Xc_eta{j}, V*Delta_v);
    C{j} = (P_all{j}{xi} \ Xc_xi{j}) * Lambda;
    
    % find factors for remaining modes ========================================
    for i = setdiff([1,2,3], j)
        % test if mode-i of eta has full rank Pi{eta}
        if rank(P_all{i}{eta}, 1e-5) == size(P_all{i}{eta}, 2)
            C{i} = P_all{i}{eta} \ Xc_eta{i};
    
        % test if mode-i of xi has full rank Pi{xi}
        elseif rank(P_all{i}{xi}, 1e-5) == size(P_all{i}{xi}, 2)
            C{i} = P_all{i}{xi} \ Xc_xi{i};
    
        else % if not search in the other xi's
            if isempty(xi_idx{i})
                error(['could not find uni-mode unique tensors for mode-' num2str(i)])
            end
            xi = xi_idx{i}(1);
            % compute CPD
            [U_xi,output] = cpd(Y{xi}, R+rankPsi(xi), opts_cpd_tensorlab);
            V = P_all{i}{eta} * (P_all{i}{xi} \ U_xi{i});
            [Delta,Delta_u,Delta_v] = find_partial_assignment(Xc_eta{i}, V, R);
            Xc_xi{i}  = U_xi{i} * Delta_v;
            % Lambda = diag( sqrt(sum(Xc_eta{i}.^2,1)) ./ sqrt(sum((V*Delta_v).^2,1)) ) * diag(sign(?max?));
            Lambda = compute_scaling_compensator(Xc_eta{i}, V*Delta_v);
            C{i} = (P_all{i}{xi} \ Xc_xi{i}) * Lambda;
        end
    end
end



% reconstruct tensor
Z = cpdgen({C{1}, C{2}, C{3}});
Psi  = cell(K,1);
Bpsi = cell(K,1);
for k=1:K
    temp = Y{k} - tmprod(Z, {P1{k}, P2{k}, P3{k}}, [1,2,3]);
    U = cpd(temp, rankPsi(k), opts_cpd_tensorlab);
    Bpsi{k} = {U{1}, U{2}, U{3}};
    Psi{k} = cpdgen(U);
end

% compute reconstruction error --------------------------
costf = 0;
for k=1:K
    err = Y{k} - tmprod(Z, {P1{k}, P2{k}, P3{k}}, [1,2,3]) - Psi{k};
    costf = costf + norm(err(:))^2;
end

end







function [Lambda]=compute_scaling_compensator(U,V)
% compute diagonal scaling compensator matrix such that 
% V * Lambda =(approx) U
assert(size(U,2) == size(V,2))
R = size(U,2);

Lambda = zeros(R,R);
for r=1:R
    % Solve linear system 
    Lambda(r,r) = (V(:,r)'*V(:,r)) \ V(:,r)'*U(:,r);
end
% Lambda = diag( sqrt(sum(Xc_eta{i}.^2,1)) ./ sqrt(sum((V*Delta_v).^2,1)) ) + sign ambiguity;
end



% simple sanity test
% addpath(genpath('YALMIP'))
% R = 3;
% U = rand(10,11);
% V = [U(:,1:R)+0.1*randn(10,R), rand(10,7-R)];
% [Delta,Delta_u,Delta_v]=find_partial_assignment(U, V, R);


function [Delta,Delta_u,Delta_v]=find_partial_assignment(U, V, R)
% finds partial assignment of R elements between columns of U and V,
% according to the angle between them

% compute matrix with angles between columns
Z = zeros(size(U,2), size(V,2));
for r = 1:size(U,2)
    for s = 1:size(V,2)
        Z(r,s) = abs(U(:,r)' * V(:,s))/ (norm(U(:,r),2) * norm(V(:,s),2));
    end
end

onesu = ones(size(U,2),1);
onesv = ones(size(V,2),1);
Delta = binvar(size(U,2), size(V,2), 'full');

obj = sum(sum(Delta .* Z));
constrs = [onesu' * Delta * onesv <= R, ...
           Delta * onesv <= onesu, ...
           Delta' * onesu <= onesv];

% ops = sdpsettings('solver','gurobi');
% ops = sdpsettings('solver','bnb','bnb.solver','fmincon');
ops = sdpsettings('solver','intlinprog','verbose',0);
% ops = sdpsettings('solver','bmibnb');

optimize(constrs,-obj,ops);
Delta = value(Delta);
% disp(['value ' num2str(R) '... ' num2str(value(obj))]) % show cost of assignement

% create permutation-selection matrices
Delta_u = Delta(:, sum(Delta,1)>0);
Delta_v = Delta(sum(Delta,2)>0, :)';
end




function [ kR ] = khatri_rao( A,B )
%khatri-rao using the matlab function for kronecker
[~, F] = size(A);

kR = [];
for f=1:F
 kR = [kR kron(A(:,f),B(:,f))];
end
end






