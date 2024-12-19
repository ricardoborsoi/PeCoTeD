function [X,iter] = solve_gen_cgkr(A, B, C, tol, maxiter, X0)
% SOLVE_GEN_CGKR    Solve a generalized Sylvester equation.
%    X = SOLVE_GEN_CGKR(A,B,C) finds the m-by-n matrix X such that
%
%       A(:,:,1)*X*B(:,:,1) + ... + A(:,:,p)*X*B(:,:,p) = C
%
%    where A are B are an m-by-m-by-p and n-by-n-p tensors, respectively, using
%    the unpreconditioned Conjugate Gradient method [Alg. 6.17, 1] on the linear system
%
%       (kron(B(:,:,1), A(:,:,1) + ... + kron(B(:,:p), A(:,:,p)) = C(:).
%
%    [X,ITER] = SOLVE_GEN_CGKR(A,B,C) returns the number ITER of iterations
%    that the method used to produce the solution X.
%
%    X = SOLVE_GEN_CGKR(A,B,C,TOL) stops the iteration when the 1-norm of the
%    residual is below TOL. The default value for this parameter is 1e-15.
%
%    X = SOLVE_GEN_CGKR(A,B,C,TOL,MAXITER) stops after at most MAXITER
%    iterations. The default value for this parameter is 10000.
%
%    X = SOLVE_GEN_CGKR(A,B,C,TOL,MAXITER,X0)) uses the m-by-n matrix X0 as
%    initial guess for the conjugate gradient iteration. The default vlaue for
%    this parameter is a matrix of zeros.
%
%    [1] Saad, Y. "Iterative Method for Sparse Linear Systems", 2nd ed., SIAM, 2003.

  if nargin < 4
    tol = 1e-15;
  end
  if nargin < 5
    maxiter = 10000;
  end

  % Check that dimensions are consistent.
  [m, m1, ell] = size(A);
  [n, n1, ell1] = size(B);
  assert(m == m1);
  assert(n == n1);
  assert(ell == ell1);

  % Initialize starting vecor, if necessary.
  if nargin < 6
    X0 = zeros(m, n);
  else
    [m1, n1] = size(X0);
    assert(m == m1);
    assert(n == n1);
  end
  X = X0;

  % Conjugate Gradient iteration
  R = C - matvec(X0, A, B);
  res_init = norm(R,'fro');
  P = R;

  beta = 0;
  Snew = norm(R, 'fro')^2;

  for iter = 1:maxiter

    S = Snew;
    AP = matvec(P, A, B);
    alpha = S / (P(:)'*AP(:));
    X = X + alpha*P;
    R = R - alpha*AP;
    Snew = norm(R, 'fro')^2; % S=R.'*R;

    beta = Snew/S;
    P = R + beta*P;
    normA = norm(R(:)) / res_init;

    if normA < tol
      break
    end

  end

  function Y = matvec(X, A, B)
    Y=zeros(size(X));
    for i = 1:ell
      Y = Y + A(:,:,i)*X*B(:,:,i);
    end
  end

end
