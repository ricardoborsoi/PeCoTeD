function [X, iter] = solve_gen_dfpm(A, B, C, tol, maxiter, optpar, dt, eta)
% SOLVE_GEN_DFPM    Solve a generalized Sylvester equation.
%    X = SOLVE_GEN_DFPM(A,B,C) uses the discrete Dynamical Functional Particle
%    Method (DFPM) to find the m-by-n matrix X such that
%
%       A(:,:,1)*X*B(:,:,1) + ... + A(:,:,p)*X*B(:,:,p) = C
%
%    where A are B are an m-by-m-by-p and n-by-n-p tensors, respectively.
%
%    [X,ITER] = SOLVE_GEN_DFPM(A,B,C) returns the number ITER of iterations
%    that the method used to produce the solution X.
%
%    [...] = SOLVE_GEN_DFPM(A,B,C,TOL) stops the iteration when the 1-norm of
%    the residual is below TOL. The default value for this parameter is 1e-15.
%
%    [...] = SOLVE_GEN_DFPM(A,B,C,TOL,MAXITER) stops after at most MAXITER
%    iterations. The default value for this parameter is 10000.
%
%    [...] = SOLVE_GEN_DFPM(A,B,C,TOL,MAXITER,OPTPAR) chooses how the timestep
%    and damping that the algorithms uses are to be chosen. The value of this
%    parameter is interpreted as follows.
%
%      * OPTPAR = 0 - The algorithm computes the optimal parameters using the
%      eigenvalues of the coefficient matrix of the linear system associated to
%      the generalized Sylvester matrix equation. The matrix is not explicitly
%      constructed when not necessary, that is, for the matrix equation in the
%      form A*X*B = C, A*X+X*B = C, and A*X*B+D1*X*D2 = C, where D1 an D2 are
%      digaonal matrices.
%
%      * OPTPAR > 0 - The algorithm computes the optimal parameters using an
%      approximation of the eigenvalues of the coefficient matrix of the linear
%      system associated to the generealized Sylveter matrix equation. This
%      approximation only requires the extreme eigenvalues of the coefficient
%      matrices on the left-hand side, and in the general case is more efficient
%      than the algorithm for OPTPAR = 0. This is the default behavior if the
%      value of OPTPAR is not specified.
%
%      * OPTPAR < 0 - The algorithms use the timestep and damping coefficient
%      specified by the user (see below).
%
%    [...] = SOLVE_GEN_DFPM(A,B,C,TOL,MAXITER,OPTPAR,DT,ETA) sets the timestep
%    to DT and the damping coefficient to ETA. The value of the seventh and
%    eight arguments is ignored unless OPTPAR < 0. The default value of DT and
%    ETA are 0.1 and 1, respectively.

  if nargin < 4
    tol = 1e-15;
  end
  if nargin < 5
    maxiter = 50000;
  end
  if nargin < 6
    optpar = 1; % Use approximation of optimal parameters.
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

  % Approximate extreme eigenvalues of the Kronecker matrix.
  if optpar >= 0
    issylvester = false;
    isstein = false;
    eigtol = 1e-1;
    if ell == 1
      i1 = 1;
      i2 = 1;
    elseif ell == 2
      % Special cases: Sylvester and Stein equation.
      iseye = @(A)(isdiag(A) && all(diag(A) == 1));
      if iseye(A(:,:,1)) && iseye(B(:,:,2))
        i1 = 2;
        i2 = 1;
        issylvester = true;
      elseif iseye(A(:,:,1)) && iseye(B(:,:,2))
        i1 = 1;
        i2 = 2;
        issylvester = true;
      elseif iseye(A(:,:,1)) && iseye(B(:,:,1))
        i1 = 1;
        i2 = 1;
        isstein = true;
      elseif iseye(A(:,:,2)) && iseye(B(:,:,2))
        i1 = 2;
        i2 = 2;
        isstein = true;
      end
    end
    if ell == 1 || issylvester || isstein % Specil cases
      if m == 1
        eigminA = A(:,:,i1);
        eigmaxA = A(:,:,i1);
      else
        eigminA = eigs(A(:,:,i1), 1, 'sm', 'tolerance', eigtol);
        eigmaxA = eigs(A(:,:,i1), 1, 'lm', 'tolerance', eigtol);
      end
      if isequal(A(:,:,i1), B(:,:,i2)) ||... % AXA = C
            isequal(A(:,:,i1), B(:,:,i2).')  % AXA.' = C
        eigminB = eigminA;
        eigmaxB = eigmaxA;
      elseif isequal(A(:,:,i1), B(:,:,i2)')  % AXA' = C with complex A
        eigminB = conj(eigminA);
        eigmaxB = conj(eigmaxA);
      else
        if n == 1
          eigminB = B(:,:,i1);
          eigmaxB = B(:,:,i1);
        else
          eigminB = eigs(B(:,:,i1), 1, 'sm', 'tolerance', eigtol);
          eigmaxB = eigs(B(:,:,i1), 1, 'lm', 'tolerance', eigtol);
        end
      end
      if ell == 1
        eigmin = eigminA * eigminB;
        eigmax = eigmaxA * eigmaxB;
      elseif issylvseter
        eigmin = eigminA + eigminB;
        eigmax = eigmaxA + eigmaxB;
      elseif isstein
        eigmin = eigminA * eigminB + 1;
        eigmax = eigmaxA * eigmaxB + 1;
      end
    else
      if optpar == 0
        M = zeros(m*n, m*n);
        for i = 1:ell
          M = M + kron(B(:, :, i), A(:, :, i));
        end
        if m*n == 1
          eigmin = M;
          eigmax = M;
        else
          eigmin = eigs(M, 1, 'sm', 'tolerance', eigtol);
          eigmax = eigs(M, 1, 'lm', 'tolerance', eigtol);
        end
      else % optpar > 0
        eigmin = 0;
        eigmax = 0;
        for i = 1:ell
          if m == 1
            lminA = A(:,:,i);
            lmaxA = A(:,:,i);
          else
            lminA = eigs(A(:,:,i), 1, 'sm', 'tolerance', eigtol);
            lmaxA = eigs(A(:,:,i), 1, 'lm', 'tolerance', eigtol);
          end
          if n == 1
            lminB = B(:,:,i);
            lmaxB = B(:,:,i);
          else
            lminB = eigs(B(:,:,i), 1, 'sm', 'tolerance', eigtol);
            lmaxB = eigs(B(:,:,i), 1, 'lm', 'tolerance', eigtol);
          end
          lmin = lminA * lminB;
          lmax = lmaxA * lmaxB;
          eigmin = eigmin + lmin;
          eigmax = eigmax + lmax;
        end
      end
    end
    lmin = sqrt(eigmin);
    lmax = sqrt(eigmax);
    dt = 2 / (lmin + lmax);
    eta = dt * lmin * lmax;
  end

  X = randn(m, n);
  Y = randn(m, n);
  F = residual(X, A, B, C);
  mynorm = 1;
  normAB = 0;
  for i = 1:ell
    normAB = normAB + norm(A(:,:,i), mynorm)*norm(B(:,:,i), mynorm);
  end
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
    Y = C;
    for i = 1:ell
      Y = Y - A(:,:,i)*X*B(:,:,i);
    end
  end

end