function [X, iter] = solve_sylv_dfpm(A, B, C, tol, maxiter, optpar, dt, eta)
% SOLVE_SYLV_DFPM    Solve the Sylvester equation.
%    X = SOLVE_SYLV_DFPM(A,B,C) uses the discrete Dynamical Functional Particle
%    Method (DFPM) to solve the Sylvester equation A*X + X*B = C.
%
%    [X,ITER] = SOLVE_SYLV_DFPM(A,B,C) returns the number ITER of iterations
%    that the method used to produce the solution X.
%
%    [...] = SOLVE_SYLV_DFPM(A,B,C,TOL) stops the iteration when the 1-norm of
%    the relative residual is below TOL. The default value for this parameter is
%    1e-15.
%
%    [...] = SOLVE_SYLV_DFPM(A,B,C,TOL,MAXITER) stops after at most MAXITER
%    iterations. The default value for this parameter is 10000.
%
%    [...] = SOLVE_SYLV_DFPM(A,B,C,TOL,MAXITER,OPTPAR) chooses how the timestep
%    and damping that the algorithms uses are to be chosen. If OPTPAR=TRUE, the
%    function computes the optimal parameters using the eigenvalues of the
%    coefficient matrix of the linear system associated to the generalized
%    Sylvester matrix equation. This is the default option.
%
%    [...] = SOLVE_SYLV_DFPM(A,B,C,TOL,MAXITER,OPTPAR,DT,ETA) sets the timestep
%    to DT and the damping coefficient to ETA. The value of the seventh and
%    eight arguments is ignored if OPTPAR=TRUE. The default value of DT and
%    ETA are 0.1 and 1, respectively.

  if nargin < 4
    tol = 1e-15;
  end
  if nargin < 5
    maxiter = 50000;
  end
  if nargin < 6
    optpar = true; % Use optimal parameters.
  end
  if optpar < 0
    if nargin < 7
      dt = 0.1;
    end
    if nargin < 8
      eta = 1;
    end
  end

  % Check that dimensions are consistent.
  [m, m1, ell] = size(A);
  [n, n1, ell1] = size(B);
  assert(m == m1);
  assert(n == n1);
  assert(ell == ell1);
  assert(ell == 1);

  % Approximate extreme eigenvalues of the Kronecker matrix.
  if optpar
    eigtol = 1e-1;
    if m == 1 % Avoid bug in eigs.
      eigminA = A;
      eigmaxA = A;
    else
      eigminA = eigs(A, 1, 'sm', 'tolerance', eigtol);
      eigmaxA = eigs(A, 1, 'lm', 'tolerance', eigtol);
    end
    if isequal(A, B) ||... % AXA = C
          isequal(A, B.')  % AXA.' = C
      eigminB = eigminA;
      eigmaxB = eigmaxA;
    elseif isequal(A, B')  % AXA' = C with complex A
      eigminB = conj(eigminA);
      eigmaxB = conj(eigmaxA);
    else
      if n == 1 % Avoid bug in eigs.
        eigminB = B;
        eigmaxB = B;
      else
        eigminB = eigs(B, 1, 'sm', 'tolerance', eigtol);
        eigmaxB = eigs(B, 1, 'lm', 'tolerance', eigtol);
      end
    end
    eigmin = eigminA + eigminB;
    eigmax = eigmaxA + eigmaxB;

    lmin = sqrt(eigmin);
    lmax = sqrt(eigmax);
    dt = 2 / (lmin + lmax);
    eta = dt * lmin * lmax;
  end

  X = randn(m, n);
  Y = randn(m, n);
  F = residual(X, A, B, C);
  mynorm = 1;
  normAB = norm(A) + norm(B);
  normC = norm(C, mynorm);
  for iter = 1:maxiter
    Y = Y*(1-eta*dt) + dt*F;
    X = X + dt*Y;
    F = residual(X, A, B, C);
    if (norm(F, mynorm) / (normAB * norm(X, mynorm) + normC) < tol)
      break;
    end
  end

  function Y = residual(X, A, B, C)
    Y = C - A*X - X*B;
  end

end