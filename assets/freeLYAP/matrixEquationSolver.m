function X = matrixEquationSolver(A, B, F, X0, flag_use_bartelsStewart)
%MATRIXEQUATIONSOLVER 
%   X = MATRIXEQUATIONSOLVER(A, B, F) solves the generalized Sylvester equation:
%
%      sum_j A{j} X B{j}^T = F, 
%
%   where A and B are cell array of nxn and mxm matrices, respectively. 
% 
%   There is no known efficient and general algorithm when the number of terms
%   in the sum is more than 2. The resulting algorithm based on Kronecker
%   products costs O((m*n)^3) operations.

% Written by Alex Townsend, Nov 2014. (alex.townsend1987@gmail.com)

% Do not form matrices above this size: 
maxSize = 10000; 
tol = 1e-8; % numerical tolerance for some of the iterative algorithms

% Initialize starting vecor, if necessary.
if nargin < 4 || isempty(X0)
    X0 = zeros(size(A{1},2), size(B{1},2));
end
if nargin < 5
    flag_use_bartelsStewart = true;
end



if ( ~isa(A, 'cell') || ~isa(B, 'cell') ) 
    error('MATRIXEQUATIONSOLVER:INPUTS:CELL',...
                      'Defining matrices are not in a cell array.'); 
end

% Check that the sizes of the cell arrays are consistent: 
if ( size(A) ~= size(B) ) 
    error('MATRIXEQUATIONSOLVER:INPUTS:SQUARE', ...
                       'Ambigous number of terms in the matrix equation.'); 
end 

% All the matrices should be square and of the same size: 
[mA, nA] = size( A{ 1 } ); 
[mB, nB] = size( B{ 1 } ); 
if ( (mA ~= nA) || (mB ~= nB) )
    error('MATRIXEQUATIONSOLVER:INPUTS:SIZES', ...
                       'Rectangular matrices are not allowed.'); 
end

% flag to use an alternative solution method
flag_use_alternative_method = false;

if ( ( max(size(A)) == 1 ) )
    
	% Solving A X B^T = F
    Y = A{1} \ F;    % Y = XB^T
    X = Y / B{1}.';  % Solve for X
    
elseif ( ( max(size(A)) == 2 ) ) 
    % If we have two terms, then use bartelsStewart(): 
    if flag_use_bartelsStewart == true
        X = bartelsStewart(A{1}, B{1}, A{2}, B{2}, F ); 
    else
        flag_use_alternative_method = true;
    end
    % test if there is any problem with the solution
    if (flag_use_bartelsStewart == true) && ( any(isnan(X(:))) || any(isinf(X(:))) )
        warning('Bartels and Stewart method failed! Resortign to an alternative method...')
        flag_use_alternative_method = true;
    end
else
    flag_use_alternative_method = true;
end



% solve using the following method if BS faield or if there are more than 3
% terms in the generalized equation
if flag_use_alternative_method == true
    % if matrix is not too big, solve using Kronecker method (faster)
    if nA * nB <= maxSize 
        X = solve_by_kronecker_method(A, B, F, nA, nB, maxSize);

    else
        % if the matrix is too big, use a conjugate gradient or some iterative algorithm
        A_tens = cell2mat(permute(A,[1,3,2]));
        B_tens = cell2mat(permute(B,[1,3,2]));

        maxiter = 10000;
        X = solve_gen_cgkr(A_tens, permute(B_tens,[2,1,3]), F, tol, maxiter, X0);

        % X = solve_gen_gbia(A_tens, permute(B_tens,[2,1,3]), F, tol); %, maxiter, approximate) % slow
%         X = solve_gen_dfpm(A_tens, permute(B_tens,[2,1,3]), F, tol); %, maxiter, optpar, dt, eta);
    end  
end

end


function [X] = solve_by_kronecker_method(A, B, F, nA, nB, maxSize)
    % Form the large Kronecker matrix. 
    % Check that the matrices are not too large: 
    if ( nA * nB > maxSize ) 
        error('MATRIXEQUATIONSOLVER:MEMORY', ... 
                        'A very large matrix will be formed.' );
    end

    % Form the large nA x nB matrix: 
    C = zeros( nA*nB ); 
    for j = 1:max(size(A))
        C = C + kron( B{j}, A{j} ); 
    end
    % Solve the large linear system: 
    x = C \ F(:); 
    % Reshape solution vector to a matrix: 
    X = reshape( x, nA, nB ); 

end