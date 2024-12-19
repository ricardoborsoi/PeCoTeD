function [Z, Psi, Bz1, Bz2, Bz3, B1, B2, B3, costf, C1, C2, C3]=solve_PECOTED_optimization(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts)
% =========================================================================
% Code for personalized tensor-based coupled data fusion with variability
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
%         which modes are coupled for the k-th image
% opts : struct with a lot of possible options (see example file)
% 
% Ricardo Borsoi
% =========================================================================

if nargin < 8
    opts = struct();
end

% option to initialize the algorithm
if ~isfield(opts, 'initOpt')
    initOpt = 'cpd';
else, initOpt = opts.initOpt; end

% in case of random intialization, whether to initialize with rand or randn
if ~isfield(opts, 'flag_init_nonnegative')
    flag_init_nonnegative = false;
else, flag_init_nonnegative = opts.flag_init_nonnegative; end

% in case of initialization of factors of Z given a priori
if ~isfield(opts, 'U_init')
    U_init = {[],[],[]};
else, U_init = opts.U_init; end

% tolerance for stopping criterion (relative loss decrease)
if ~isfield(opts, 'mytol')
    mytol = 1e-5;
else, mytol = opts.mytol; end

% whether to initialize with complex values on a random initialization
if ~isfield(opts, 'initcomp')
    initcomp = 0;
else, initcomp = opts.initcomp; end

% number of ALS iterations
if ~isfield(opts, 'num_iters')
    num_iters = 100;
else, num_iters = opts.num_iters; end

% whether to use Bartels and Stewart method to solve the generalized
% sylvester equations (the code can be numerically unstable!)
if ~isfield(opts, 'flag_use_bartelsStewart')
    flag_use_bartelsStewart = true;
else, flag_use_bartelsStewart = opts.flag_use_bartelsStewart; end

% initializatio of the CPD algorithms used to initialize the method
if ~isfield(opts, 'initialCPD_initoption')
    initialCPD_initoption = 'auto'; % options: 'auto', 'random', 'gesvd'
else, initialCPD_initoption = opts.initialCPD_initoption; end


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


% if the matrices P1, P2, P3 are not cell arrays, they will be traated as
% constant in k, i.e., P1{1}=...=P1{K}=P1
flag_P1_constant = false;
flag_P2_constant = false;
flag_P3_constant = false;
if ~iscell(P1)
    flag_P1_constant = true; 
    P1_const = P1;
    P1 = cell(K,1);
    for k=1:K, P1{k} = P1_const; end
end
if ~iscell(P2)
    flag_P2_constant = true; 
    P2_const = P2;
    P2 = cell(K,1);
    for k=1:K, P2{k} = P2_const; end
end
if ~iscell(P3)
    flag_P3_constant = true; 
    P3_const = P3;
    P3 = cell(K,1);
    for k=1:K, P3{k} = P3_const; end
end

% get latent tensor size
M1 = size(P1{1},2);
M2 = size(P2{1},2);
M3 = size(P3{1},2);


% Initialize Bz randomly or with a CPD of interpolated and averaged data
if strcmp(initOpt, 'random')
    if flag_init_nonnegative
        Bz1 = rand(M1, rankZ) + initcomp*1i * randn(M1, rankZ);
        Bz2 = rand(M2, rankZ) + initcomp*1i * randn(M2, rankZ);
        Bz3 = rand(M3, rankZ) + initcomp*1i * randn(M3, rankZ);
    else
        Bz1 = randn(M1, rankZ) + initcomp*1i * randn(M1, rankZ);
        Bz2 = randn(M2, rankZ) + initcomp*1i * randn(M2, rankZ);
        Bz3 = randn(M3, rankZ) + initcomp*1i * randn(M3, rankZ);
    end

elseif strcmp(initOpt, 'cpd')
    X = zeros(M1,M2,M3);
    for k=1:K
        % this is a cheap and dirty initialization: try to regress C from
        % each Yk, take the mean, and compute the CPD
        temp = pinv_per_mode(Y{k}, P1{k}, P2{k}, P3{k});
        X    = X + (1/K) * temp;
    end
    U = cpd(X, rankZ, opts_cpd_tensorlab);
    Bz1 = U{1};
    Bz2 = U{2};
    Bz3 = U{3};
    Z_init = cpdgen({Bz1, Bz2, Bz3});
    

elseif strcmp(initOpt, 'apriori')
    Bz1 = U_init{1};
    Bz2 = U_init{2};
    Bz3 = U_init{3};
    Z_init = cpdgen({Bz1, Bz2, Bz3});


elseif strcmp(initOpt, 'semialgebraic')
    opts_semialgebraic.initialCPD_initoption = initialCPD_initoption; % pass same option for the CPD initialization
    opts_semialgebraic.alg_use_regression    = opts.alg_use_regression; % use a regression approach to estimate some factors
    try 
        [Z_init, Psi_init, zfacsinit] = solve_PECOTED_semialgebraic(Y, P1, P2, P3, rankZ, rankPsi, Gamma, opts_semialgebraic);
    catch me
        error('Error found during the semialgebraic initialization, try "random" or "cpd" instead.')
    end
    Bz1 = zfacsinit{1};
    Bz2 = zfacsinit{2};
    Bz3 = zfacsinit{3};
else
    error(['Unknown initialization type:' initOpt])
end


% Initialize Bk randomly or with a CPD of the residual
B1 = cell(K,1);
B2 = cell(K,1);
B3 = cell(K,1);
for k=1:K
    if strcmp(initOpt, 'random')
        if flag_init_nonnegative
            B1{k} = rand(size(P1{k},1), rankPsi(k)) + initcomp*1i * randn(size(P1{k},1), rankPsi(k));
            B2{k} = rand(size(P2{k},1), rankPsi(k)) + initcomp*1i * randn(size(P2{k},1), rankPsi(k));
            B3{k} = rand(size(P3{k},1), rankPsi(k)) + initcomp*1i * randn(size(P3{k},1), rankPsi(k));
        else
            B1{k} = randn(size(P1{k},1), rankPsi(k)) + initcomp*1i * randn(size(P1{k},1), rankPsi(k));
            B2{k} = randn(size(P2{k},1), rankPsi(k)) + initcomp*1i * randn(size(P2{k},1), rankPsi(k));
            B3{k} = randn(size(P3{k},1), rankPsi(k)) + initcomp*1i * randn(size(P3{k},1), rankPsi(k));
        end
    elseif strcmp(initOpt, 'cpd') || strcmp(initOpt, 'semialgebraic') || strcmp(initOpt, 'apriori')
        U = cpd(Y{k} - mlprod(Z_init, P1{k}, P2{k}, P3{k}), rankPsi(k), opts_cpd_tensorlab);
        B1{k} = U{1};
        B2{k} = U{2};
        B3{k} = U{3};
    end
end


% initialize Cj{k}
C1 = cell(K,1);
C2 = cell(K,1);
C3 = cell(K,1);
for k=1:K
    C1{k} = [P1{k}*Bz1, B1{k}];
    C2{k} = [P2{k}*Bz2, B2{k}];
    C3{k} = [P3{k}*Bz3, B3{k}];
end



% run ALS -----------------------------------------------------------------
costf_old = 1e10;
costf_vec = zeros(num_iters,1);
for iter=1:num_iters
    
    % optimize w.r.t. B_Z_1 -------------------
    F = zeros(rankZ, M1);
    i = 0; AA = cell(1); BB = cell(1);
    for k=1:K
        if ismember(1, Gamma{k})
            i = i+1;
            % temp = khatri_rao(P3{k}*Bz3, P2{k}*Bz2);
            temp = khatri_rao(C3{k}(:,1:R), C2{k}(:,1:R));
            AA{i} = temp'*temp;
            if ~flag_P1_constant || i == 1
                BB{i} = P1{k}'*P1{k};
            end
            F = F + temp'*(tens2mat(Y{k},1)' - khatri_rao(B3{k},B2{k})*B1{k}')*P1{k};
        end
    end
    % solve matrix equation
    if ~flag_P1_constant
        Bz1 = matrixEquationSolver(AA, BB, F, Bz1', flag_use_bartelsStewart)';
    else
        temp = AA{1};
        for i=2:length(AA), temp = temp + AA{i}; end
        Bz1 = (temp \ (F / BB{1}))'; 
    end
    
    
    
    
    % optimize w.r.t. C1_k(:,1:R) -------------------
    for k=1:K
        if ismember(1, Gamma{k})
            C1{k}(:,1:R) = P1{k}*Bz1;
        else
            temp = tens2mat(Y{k},1)' - khatri_rao(B3{k},B2{k})*B1{k}';
            C1{k}(:,1:R) = (khatri_rao(C3{k}(:,1:R), C2{k}(:,1:R)) \ temp)';
        end
    end
    
    
    
    
    % optimize w.r.t. B_Z_2 -------------------
    F = zeros(rankZ, M2);
    i = 0; AA = cell(1); BB = cell(1);
    for k=1:K
        if ismember(2, Gamma{k})
            i = i+1;
            % temp = khatri_rao(P3{k}*Bz3, P1{k}*Bz1);
            temp = khatri_rao(C3{k}(:,1:R), C1{k}(:,1:R));
            AA{i} = temp'*temp;
            if ~flag_P2_constant || i == 1
                BB{i} = P2{k}'*P2{k};
            end
            F = F + temp'*(tens2mat(Y{k},2)' - khatri_rao(B3{k},B1{k})*B2{k}')*P2{k};
        end
    end
    % solve matrix equation
    if ~flag_P2_constant
        Bz2 = matrixEquationSolver(AA, BB, F, Bz2', flag_use_bartelsStewart)';
    else
        temp = AA{1};
        for i=2:length(AA), temp = temp + AA{i}; end
        Bz2 = (temp \ (F / BB{1}))'; 
    end
    
    
    
    
    % optimize w.r.t. C2_k(:,1:R) -------------------
    for k=1:K
        if ismember(2, Gamma{k})
            C2{k}(:,1:R) = P2{k}*Bz2;
        else
            temp = tens2mat(Y{k},2)' - khatri_rao(B3{k},B1{k})*B2{k}';
            C2{k}(:,1:R) = (khatri_rao(C3{k}(:,1:R), C1{k}(:,1:R)) \ temp)';
        end
    end
    
    
    
    
    % optimize w.r.t. B_Z_3 -------------------
    F = zeros(rankZ, M3);
    i = 0; AA = cell(1); BB = cell(1);
    for k=1:K
        if ismember(3, Gamma{k})
            i = i+1;
            % temp = khatri_rao(P2{k}*Bz2, P1{k}*Bz1);
            temp = khatri_rao(C2{k}(:,1:R), C1{k}(:,1:R));
            AA{i} = temp'*temp;
            if ~flag_P3_constant || i == 1
                BB{i} = P3{k}'*P3{k};
            end
            F = F + temp'*(tens2mat(Y{k},3)' - khatri_rao(B2{k},B1{k})*B3{k}')*P3{k};
        end
    end
    % solve matrix equation
    if ~flag_P3_constant
        Bz3 = matrixEquationSolver(AA, BB, F, Bz3', flag_use_bartelsStewart)';
    else
        temp = AA{1};
        for i=2:length(AA), temp = temp + AA{i}; end
        Bz3 = (temp \ (F / BB{1}))'; 
    end
    
    
    
    
    
    % optimize w.r.t. C3_k(:,1:R) -------------------
    for k=1:K
        if ismember(3, Gamma{k})
            C3{k}(:,1:R) = P3{k}*Bz3;
        else
            temp = tens2mat(Y{k},3)' - khatri_rao(B2{k},B1{k})*B3{k}';
            C3{k}(:,1:R) = (khatri_rao(C2{k}(:,1:R), C1{k}(:,1:R)) \ temp)';
        end
    end
    
    
    
    
    
    % optimize w.r.t. B_k_1, B_Z_2, B_Z_3 -------------------
    for k=1:K
        temp  = khatri_rao(B3{k}, B2{k});
        X     = tens2mat(Y{k},1)' - khatri_rao(P3{k}*Bz3, P2{k}*Bz2)*Bz1'*P1{k}';
        B1{k} = (temp \ X)';
        
        temp  = khatri_rao(B3{k}, B1{k});
        X     = tens2mat(Y{k},2)' - khatri_rao(P3{k}*Bz3, P1{k}*Bz1)*Bz2'*P2{k}';
        B2{k} = (temp \ X)';
        
        temp  = khatri_rao(B2{k}, B1{k});
        X     = tens2mat(Y{k},3)' - khatri_rao(P2{k}*Bz2, P1{k}*Bz1)*Bz3'*P3{k}';
        B3{k} = (temp \ X)';
    end
    
    
    
    % update w.r.t. Cj_k(:,R+1:R+Lk) -------------------
    for k=1:K
        C1{k}(:,R+1:end) = B1{k};
        C2{k}(:,R+1:end) = B2{k};
        C3{k}(:,R+1:end) = B3{k};
    end
    
    
    
    
    % check if terminaton criterion is satisfied --------------------------
    costf = 0;
    for k=1:K
        err = Y{k} - cpdgen({P1{k}*Bz1, P2{k}*Bz2, P3{k}*Bz3}) ...
                   - cpdgen({B1{k}, B2{k}, B3{k}});
        costf = costf + norm(err(:))^2;
    end
    
    if abs(costf_old - costf)/costf_old < mytol
        break
    end
    costf_old = costf;
    
%     disp(costf)
%     costf_vec(iter) = costf;
%     plot(costf_vec), pause(0.1)
    
end



% normalize the mode-1 and mode-2 factors to have unit norm
[Bz1, Bz2, Bz3] = normalize_factors(Bz1, Bz2, Bz3);
for k=1:K
    [B1{k}, B2{k}, B3{k}] = normalize_factors(B1{k}, B2{k}, B3{k});
    [C1{k}, C2{k}, C3{k}] = normalize_factors(C1{k}, C2{k}, C3{k});
end


% reconstruct Z
Z = cpdgen({Bz1, Bz2, Bz3});

% reconstruct Psis
Psi = cell(K,1);
for k=1:K
    Psi{k} = cpdgen({B1{k}, B2{k}, B3{k}});
end

end



%*******************************************************************************
function [A,B,C] = normalize_factors(A,B,C)
% normalize the mode-1 and mode-2 factors to have unit norm
rankZ         = size(A,2);
temp_scalings = zeros(rankZ,1);
for i=1:rankZ
    temp_scalings(i) = norm(A(:,i));
    temp_scalings(i) = temp_scalings(i) * norm(B(:,i));
    A(:,i) = A(:,i) / norm(A(:,i)); % unit norm
    B(:,i) = B(:,i) / norm(B(:,i)); % unit norm
    % make C (third mode) absorb scaling factors
    C(:,i) = C(:,i) * temp_scalings(i);
end
end



%*******************************************************************************
function X = pinv_per_mode(Y,P1,P2,P3)
M1 = size(P1,2);
M2 = size(P2,2);
M3 = size(P3,2);
temp = Y;
temp = mldivide(P1'*P1 + 1e-6*speye(M1), P1' * tens2mat(temp,1));
temp = mat2tens(temp, [M1,size(Y,2),size(Y,3)], 1);
temp = mldivide(P2'*P2 + 1e-6*speye(M2), P2' * tens2mat(temp,2));
temp = mat2tens(temp, [M1,M2,size(Y,3)], 2);
temp = mldivide(P3'*P3 + 1e-6*speye(M3), P3' * tens2mat(temp,3));
temp = mat2tens(temp, [M1,M2,M3], 3);
X    = temp;
end



%*******************************************************************************
function X_out = mlprod(X,B1,B2,B3)
% MLPROD full multilinear product
X_out = tmprod(X,B1,1);
X_out = tmprod(X_out,B2,2);
X_out = tmprod(X_out,B3,3);
end

%*******************************************************************************
function [U1,U2,U3,S,S1,S2,S3] = mlsvd3(X,size_core)
%MLSVD3 Multilinear singular value decomposition of a third-order tensor.
[I1,I2,I3]=size(X);
[U1,S1,temp]=svd(reshape(X,I1,I3*I2),'econ'); S1=diag(S1);
[U2,S2,temp]=svd(reshape(permute(X,[2 3 1]),I2,I1*I3),'econ'); S2=diag(S2);
[U3,S3,temp]=svd(reshape(permute(X,[3 1 2]),I3,I2*I1),'econ'); S3=diag(S3);
if nargin==2
    U1=U1(:,1:min(size_core(1),I2*I3));
    U2=U2(:,1:min(size_core(2),I1*I3));
    U3=U3(:,1:min(size_core(3),I1*I2));
end
S=tmprod(tmprod(tmprod(X,U1',1),U2',2),U3',3);
end

%*******************************************************************************
function X_out = tmprod(X,U,mode)
%TMPROD mode-n tensor-matrix product.
[I,J,K]=size(X);
[M,N]=size(U);
if (mode~=1) && (mode~=2) && (mode~=3)
    error('The input variable mode should be 1, 2 or 3')
end
if N~=size(X,mode) 
    error(['The number of columns of the input matrix should be equal to dimension ',int2str(mode),' of the input tensor'])
end    
if mode==1
    X_out = reshape(U*reshape(X,I,J*K) ,M,J,K);      
elseif mode==2
    X_out = permute(reshape (U*reshape(permute(X,[2 1 3]),J,I*K), M,I,K),[2 1 3]);        
elseif mode==3
    X_out = permute(reshape (U*reshape(permute(X,[3 1 2]),K,I*J), M,I,J),[2 3 1]);
end
end

%*******************************************************************************       
function C = khatri_rao(A,B)
%KR Khatri-Rao product.
[I R1]=size(A); J=size(B,1); 
C=zeros(I*J,R1);
for j=1:R1
    C(:,j)=reshape(B(:,j)*A(:,j).',I*J,1);
end
end

%*******************************************************************************
function Mat = khatri_rao_part(B,C,partB,partC)
%KR_PART Partition-Wise Kronecker product      
[J M]=size(B);
[K N]=size(C);
if (sum(partB)~=M) 
    error(['Error: a matrix with ',int2str(M),' columns can not be partitioned in such a way'])
end
if (sum(partC)~=N) 
    error(['Error: a matrix with ',int2str(N),' columns can not be partitioned in such a way'])
end
if length(partB)~=length(partC)
     error('Error: the 2 input matrices do not have the same number of blocks')
end

indB=[0 cumsum(partB)];
indC=[0 cumsum(partC)];
indMat=[0 cumsum(partB.*partC)];

 Mat=zeros(J*K,sum(partB.*partC));
 for i=1:length(partC)  
     Mat(:,indMat(i)+1:indMat(i+1))=fast_kron( B(:,indB(i)+1:indB(i+1)) , C(:,indC(i)+1:indC(i+1)));
 end
end
         
%*******************************************************************************         
function C = fast_kron (A,B)
%FAST_KRON fast kronecker product
[I,L1]=size(A);
[J,L2]=size(B);

if (L1==1) && (L2==1)
     C=reshape(B*A.',I*J,1);
elseif (L1==1) && (L2>1)
    Bt=B.'; 
    C=reshape(Bt(:)*A.',L2,I*J).';
elseif (L2==1) && (L1>1)
     C=reshape(B*A(:).',I*J,L1);
else
     C=reshape(permute(reshape(B(:)*A(:).',[J,L2,I,L1]),[1 3 2 4]),[I*J,L1*L2]);
end
end






