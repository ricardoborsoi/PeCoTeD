function [X, iter] = solve_gen_gbia(A, B, C, tol, maxiter, approximate)
% SOLVE_GEN_GBIA    Solve a generalized Sylvester equation.
%    X = SOLVE_GEN_GBIA(A,B,C) uses the Gradient Based Iterative Algorithm of Ding
%    and Chen [1] to find the m-by-n matrix X such that
%
%       A(:,:,1)*X*B(:,:,1) + ... + A(:,:,p)*X*B(:,:,p) = C
%
%    where A are B are an m-by-m-by-p and n-by-n-p tensors, respectively.
%
%    [X,ITER] = SOLVE_GEN_GBIA(A,B,C) returns the number ITER of iterations
%    that the method used to produce the solution X.
%
%    X = SOLVE_GEN_GBIA(A,B,C,TOL) stops the iteration when the 1-norm of the
%    residual is below TOL. The default value for this parameter is 1e-15.
%
%    X = SOLVE_GEN_GBIA(A,B,C,TOL,MAXITER) stops after at most MAXITER
%    iterations. The default value for this parameter is 10000.
%
%    X = SOLVE_GEN_GBIA(A,B,C,TOL,MAXITER,APPROXIMATE) uses the optimal
%    parameter if APPROXIMATE is set to false, and an approximation to them if
%    APPROXIMATE is set to true, p > 1, and the equation is not a Sylvester
%    equation. This parmeter is set to true by default.
%
%    [1] Ding, F. and Chen, T., "Gradient Based Algorithms for Solving Class of
%        Matrix Equations." IEEE Trans. Autom. Control 50.8 (2005): 1216-1221.

  if nargin < 4
    tol = 1e-15;
  end
  if nargin < 5
    maxiter = 50000;
  end
  if nargin < 6
    approximate = true; % Use approximation of optimal parameters.
  end

  % Check that dimensions are consistent.
  [m, m1, ell] = size(A);
  [n, n1, ell1] = size(B);
  assert(m == m1);
  assert(n == n1);
  assert(ell == ell1);

  % Compute parameter mu in [Eq. (32), 1].
  mu = 0;
  opts.issym = true;
  opts.tol = 1e-5;
  for j = 1:ell
    fA = @(x)(A(:,:,j)*(A(:,:,j)'*x));
    fB = @(x)(B(:,:,j)'*(B(:,:,j)'*x));
    mu = mu + eigs(fA, m, 1, 'lm', opts) * eigs(fB, n, 1, 'lm', opts);
  end
  mu = 1 / mu;

  % Iterative method.
  X = 1e-6 * ones(m,n);
  normAB = 0;
  for i = 1:ell
    normAB = normAB + norm(A(:,:,i), 1)*norm(B(:,:,i),1);
  end
  normC = norm(C, 1);
  for iter = 1:maxiter
    F = residual(X, A, B, C);
    for i = 1:ell
      Y(:,:,i) = X + mu*A(:,:,i).' * F * B(:,:,i).';
    end
    X = sum(Y, 3) / ell;
    if (norm(F,1) / (normAB * norm(X, 1) + normC) < tol)
      break;
    end
  end

  function Y = residual(X, A, B, C)
    Y = C;
    for i = 1:ell
      Y = Y - A(:,:,i)*X*B(:,:,i);
    end
  end

end