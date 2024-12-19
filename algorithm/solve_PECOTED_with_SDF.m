function [Z, Psi, C, D, costf] = solve_PECOTED_with_SDF(Y, rankZ, rankPsi, P1, P2, P3, initC, initD)
% =========================================================================
% Code for Personalized tensor-based coupled data fusion with variability using 
% the structured data fusion (SDF) framework or Sorber et al.
%
% Thhis code implements algorithms related to the paper:
%   Personalized Coupled Tensor Decomposition for Multimodal Data Fusion: Uniqueness and Algorithms
%   R.A. Borsoi, K. Usevich, D. Brie, T. Adali
%   IEEE Transactions on Signal Processing, 2024.
% 
% 
% Y  : cell array of order-3 tensors 
% P1, P2, P3 : cell array 
% initC, initD : initializations for the factors (optional)
% 
% Ricardo Borsoi
% =========================================================================

% see if initialization is provided
flag_initC_random = false;
flag_initD_random = false;
if nargin < 8
    flag_initD_random = true;
    if nargin < 7
        flag_initC_random = true;
    end
end


% get latent tensor size
M1 = size(P1{1},2);
M2 = size(P2{1},2);
M3 = size(P3{1},2);

K = length(Y); % number of tensors to fuse

assert(K>1 && K<4);
if length(rankPsi == 1)
    rankPsi = rankPsi * ones(K,1);
end



model = struct;

% initialize common components
if flag_initC_random == true
    model.variables.C1 = randn(M1, rankZ);
    model.variables.C2 = randn(M2, rankZ);
    model.variables.C3 = randn(M3, rankZ);
else
    model.variables.C1 = initC{1};
    model.variables.C2 = initC{2};
    model.variables.C3 = initC{3};
end


% model tensor 1 ----------------------------------

if flag_initD_random == true
    model.variables.D11 = randn(size(Y{1},1), rankPsi(1));
    model.variables.D12 = randn(size(Y{1},2), rankPsi(1));
    model.variables.D13 = randn(size(Y{1},3), rankPsi(1));
else
    model.variables.D11 = initD{1}{1};
    model.variables.D12 = initD{1}{2};
    model.variables.D13 = initD{1}{3};
end

matprodP11 = @(z,task) struct_matvec(z,task,P1{1},[]);
matprodP12 = @(z,task) struct_matvec(z,task,P2{1},[]);
matprodP13 = @(z,task) struct_matvec(z,task,P3{1},[]);

model.factors.X11 = {{'C1', matprodP11}, 'D11'};
model.factors.X12 = {{'C2', matprodP12}, 'D12'};
model.factors.X13 = {{'C3', matprodP13}, 'D13'};


% model tensor 2 ----------------------------------
if flag_initD_random == true
    model.variables.D21 = randn(size(Y{2},1), rankPsi(2));
    model.variables.D22 = randn(size(Y{2},2), rankPsi(2));
    model.variables.D23 = randn(size(Y{2},3), rankPsi(2));
else
    model.variables.D21 = initD{2}{1};
    model.variables.D22 = initD{2}{2};
    model.variables.D23 = initD{2}{3};
end

matprodP21 = @(z,task) struct_matvec(z,task,P1{2},[]);
matprodP22 = @(z,task) struct_matvec(z,task,P2{2},[]);
matprodP23 = @(z,task) struct_matvec(z,task,P3{2},[]);

model.factors.X21 = {{'C1', matprodP21}, 'D21'};
model.factors.X22 = {{'C2', matprodP22}, 'D22'};
model.factors.X23 = {{'C3', matprodP23}, 'D23'};


% model tensor 3 ----------------------------------

if K == 3
    if flag_initD_random == true
        model.variables.D31 = randn(size(Y{3},1), rankPsi(3));
        model.variables.D32 = randn(size(Y{3},2), rankPsi(3));
        model.variables.D33 = randn(size(Y{3},3), rankPsi(3));
    else
        model.variables.D31 = initD{3}{1};
        model.variables.D32 = initD{3}{2};
        model.variables.D33 = initD{3}{3};
    end

    matprodP31 = @(z,task) struct_matvec(z,task,P1{3},[]);
    matprodP32 = @(z,task) struct_matvec(z,task,P2{3},[]);
    matprodP33 = @(z,task) struct_matvec(z,task,P3{3},[]);

    model.factors.X31 = {{'C1', matprodP31}, 'D31'};
    model.factors.X32 = {{'C2', matprodP32}, 'D32'};
    model.factors.X33 = {{'C3', matprodP33}, 'D33'};
end

model.factorizations.tensor1.data = Y{1};
model.factorizations.tensor1.cpd  = {'X11', 'X12', 'X13'};

model.factorizations.tensor2.data = Y{2};
model.factorizations.tensor2.cpd  = {'X21', 'X22', 'X23'};

if K == 3
    model.factorizations.tensor3.data = Y{3};
    model.factorizations.tensor3.cpd  = {'X31', 'X32', 'X33'};
end


% sdf_check(model, 'print'); % recommended
sdf_check(model); % recommended
[sol,ooutput] = sdf_nls(model);

% sol.variables
% sol.factors
% costf = ooutput.fval(end);

C{1} = sol.variables.C1;
C{2} = sol.variables.C2;
C{3} = sol.variables.C3;

% reconstruct tensor
Z   = cpdgen({C{1}, C{2}, C{3}});
Psi = cell(K,1);
D   = cell(K,1);

D{1}   = {sol.variables.D11, sol.variables.D12, sol.variables.D13};
D{2}   = {sol.variables.D21, sol.variables.D22, sol.variables.D23};
Psi{1} = cpdgen(D{1});
Psi{2} = cpdgen(D{2});
if K == 3
    D{3}   = {sol.variables.D31, sol.variables.D32, sol.variables.D33};
    Psi{3} = cpdgen(D{3});
end

% compute reconstruction error --------------------------
costf = 0;
for k=1:K
    err = Y{k} - tmprod(Z, {P1{k}, P2{k}, P3{k}}, [1,2,3]) - Psi{k};
    costf = costf + norm(err(:))^2;
end

