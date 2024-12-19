function outstruct = test_gensylv_case_fun(m, n, p, tol, maxiter, etaA, etaB)
% TEST_GENSYLV_CASE_FUN    Generate generalized Sylvester equation and solve it.
%    TEST_GENSYLV_CASE_FUN(M,N,P) generates the right-hand side and the matrix
%    coefficients of a generalized Sylvester matrix equation
%
%        A(:,:,1)*X*B(:,:,1) + ... + A(:,:,P)*X*B(:,:,P) = C
%
%    where A are B are an M-by-M-by-P and N-by-N-by-P tensors, respectively.
%
%    The function solves the equation using the Discrete Dynamical Functional
%    Particle Method and the Gradient Based Iterative Algorithm, and reports
%    some information about the generated problem as well as the solution
%    process.
%
%    TEST_GENSYLV_CASE_FUN(M,N,P,TOL) runs the iteravive methods with the value
%    VAL used as tolerance threshold for the stopping criterion. The default
%    value of this parameter is 1e-15.
%
%    TEST_GENSYLV_CASE_FUN(M,N,P,TOL,MAXITER) specifies the maximum number of
%    iterations the iterative algorithms are allowed to perform. The default
%    value of this parameter is 10000.
%
%    TEST_GENSYLV_CASE_FUN(M,N,P,TOL,MAXITER,MAXEIGA,MAXEIGB) allows the user to
%    choose the distribution of the positive real eigenvalues of the coefficints
%    of the matrix equation. The slices of A and B will have eigenvalues between
%    1 and EIGRATA and between 1 and EIGRATB, respectively, and will be
%    simultaneously diagonalizable, which ensures that the matrix coefficeint of
%    the Kronecker linear system will have real positive eigenvalues as well.
%    Both EIGRATA and EIGRATB must be positive real numbers larger than 1. They
%    default to 10 if not specified.
%
%    See also SOLVE_GEN_DFPM, SOLVE_GEN_GBIA.

  if nargin < 3 || nargin > 7
    error('This function requires three to seven arguments.');
  end
  if nargin < 4
    tol = 1e-15;
  end
  if nargin < 5
    maxiter = 10000;
  end
  if nargin < 6
    etaA = 1e1;
  end
  if nargin < 7
    etaB = 1e1;
  end

  % Generate coefficients and reference solution.
  [Xtrue, A, B, C, M] = generate_test_case(m, n, p, etaA, etaB);

  tinit = tic;
  X_kron = solve_gen_kron(A, B, C);
  tkron = toc(tinit);
  if tkron < 1
    tkron = timeit(@()solve_gen_kron(A, B, C));
  end

  tinit = tic;
  [X_dfpm, iter_dfpm] = solve_gen_dfpm(A, B, C, tol, maxiter, 0);
  tdfpm = toc(tinit);
  if tdfpm < 1
    tdfpm = timeit(@()solve_gen_dfpm(A, B, C, tol, maxiter, 0));
  end

  dt = 0.1;
  eta = 3;
  tinit = tic;
  [X_dfpmw, iter_dfpmw] = solve_gen_dfpm(A, B, C, tol, maxiter, 1, dt, eta);
  tdfpmw = toc(tinit);
  if tdfpmw < 1
    tdfpmw = timeit(@()solve_gen_dfpm(A, B, C, tol, maxiter, 1, dt, eta));
  end

  tinit = tic;
  [X_gbia, iter_gbia] = solve_gen_gbia(A, B, C, tol, maxiter);
  tgbia = toc(tinit);
  if tgbia < 1
    tgbia = timeit(@()solve_gen_gbia(A, B, C, tol, maxiter));
  end

  normind = 1;
  outstruct.m = m;
  outstruct.n = n;

  outstruct.A = A;
  outstruct.B = B;

  outstruct.condM = cond(M);

  outstruct.res_kron = norm(residual(X_kron, A, B, C), normind);
  outstruct.res_dfpm = norm(residual(X_dfpm, A, B, C), normind);
  outstruct.res_dfpmw = norm(residual(X_dfpmw, A, B, C), normind);
  outstruct.res_gbia = norm(residual(X_gbia, A, B, C), normind);

  outstruct.err_kron = norm(Xtrue - X_kron, normind) / norm(Xtrue, normind);
  outstruct.err_dfpm = norm(Xtrue - X_dfpm, normind) / norm(Xtrue, normind);
  outstruct.err_dfpmw = norm(Xtrue - X_dfpmw, normind) / norm(Xtrue, normind);
  outstruct.err_gbia = norm(Xtrue - X_gbia, normind) / norm(Xtrue, normind);

  outstruct.iter_dfpm = iter_dfpm;
  outstruct.iter_dfpmw = iter_dfpmw;
  outstruct.iter_gbia = iter_gbia;

  outstruct.tkron = tkron;
  outstruct.tdfpm = tdfpm;
  outstruct.tdfpmw = tdfpmw;
  outstruct.tgbia = tgbia;

  fprintf(['\n----------------------------------\n',...
           '----------------------------------\n',...
           '----------------------------------\n',...
           'Problem size:        %5d x %5d\n',...
           'Condition number of M:    %.2e\n',...
           '----------------------------------\n',...
           'Relative residual:        %.2e   %.2e   %.2e   %.2e\n',...
           'Relative forward error:   %.2e   %.2e   %.2e   %.2e\n',...
           'Number of iterations:        ---     %8d   %8d   %8d\n',...
           'Time elapsed (seconds):   %.2e   %.2e   %.2e   %.2e\n'],...
          outstruct.m, outstruct.n, outstruct.condM,...
          outstruct.res_kron, outstruct.res_dfpm, outstruct.res_dfpmw, outstruct.res_gbia,...
          outstruct.err_kron, outstruct.err_dfpm, outstruct.err_dfpmw, outstruct.err_gbia,...
          outstruct.iter_dfpm, outstruct.iter_dfpmw, outstruct.iter_gbia,...
          outstruct.tkron, outstruct.tdfpm, outstruct.tdfpmw, outstruct.tgbia);

  function [Xtrue, A, B, C, M] = generate_test_case(m, n, p, cond1, cond2)
    A = zeros(m, m, p);
    B = zeros(n, n, p);
    Xtrue = randn(m, n);
    C = zeros(m, n);
    A = randeig(m, p, cond1);
    B = randeig(n, p, cond2);

    if nargout > 4
      M = zeros(m*n, m*n);
    end
    for i = 1:p
      if nargout > 4
        M = M + kron(B(:,:,i).', A(:,:,i));
      end
      C = C + A(:,:,i) * Xtrue * B(:,:,i);
    end
  end

  function Y = residual(X, A, B, C)
    [m,m1,ell] = size(A);
    [n,n1,ell1] = size(B);
    assert(m == m1);
    assert(n == n1);
    assert(ell == ell1);

    Y = C;
    for i = 1:ell
      Y = Y - A(:,:,i)*X*B(:,:,i);
    end
  end

end