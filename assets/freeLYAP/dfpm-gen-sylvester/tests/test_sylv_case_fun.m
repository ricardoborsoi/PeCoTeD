function outstruct = test_sylv_case_fun(m, n, condA, condB, tol, maxit)
% TEST_SYLV_CASE_FUN   Generate a Sylvester equation and solve it.
%    TEST_SYLV_CASE_FUN(M,N) generates the right-hand side and the matrix
%    coefficients of the matrix equation AX + XB = C, where the M-by-M matrix A
%    and the N-by-N matrix B are positive definite, and X and C are M-by-N
%    matrices. The function solves the Sylvester equation using the Discrete
%    Dynamical Functional Particle Method and reports some information about the
%    generated problem as well as the solution process.
%
%    TEST_SYLV_CASE_FUN(M,N,CONDA,CONDB) allows the user to choose the 2-norm
%    condition number of the matrices A and B using the real arguments CONDA and
%    CONDB, respectively. If the argument is positive, then it is the 2-norm
%    condition number of the corresponding matrix coefficient. Otherwise,
%    the corresponding matrix coefficient is a diagonally dominant matrix with
%    normally distributed non-diagonal entries.
%
%    See also RSVD_PSD, SOLVE_GEN_DFPM.

  if nargin < 2 || nargin > 6
    error('This function requires two to six arguments.');
  end
  if nargin < 3
    condA = 1e1; % -1;
  end
  if nargin < 4
    condB = 1e1; % -1;
  end
  if nargin < 5
    tol = 1e-15;
  end
  if nargin < 6
    maxit = 10000;
  end

  [Xtrue, A, B, C] = generate_sylv_test_case(m, n, condA, condB);

  tinit = tic;
  Xsylv = sylvester(A, B, C);
  tsylvester = toc(tinit);
  if tsylvester < 1
    tsylvester = timeit(@()sylvester(A, B, C));
  end

  tinit = tic;
  [Xdfpm, iter_dfpm] = solve_sylv_dfpm(A, B, C, tol, maxit);
  tdfpm = toc(tinit);
  if tdfpm < 1
    tdfpm = timeit(@()solve_sylv_dfpm(A, B, C, tol, maxit));
  end

  tinit = tic;
  [Xdfpmr, iter_dfpmr] = solve_sylv_dfpm_refinement(A, B, C, tol, maxit);
  tdfpmr = toc(tinit);
  if tdfpmr < 1
    tdfpmr = timeit(@()solve_sylv_dfpm_refinement(A, B, C, tol, maxit));
  end

  eigtol = 1e-5;
  if m == 1
    eigminA = A;
    eigmaxA = A;
  else
    eigminA = eigs(A, 1, 'sm', 'tolerance', eigtol);
    eigmaxA = eigs(A, 1, 'lm', 'tolerance', eigtol);
  end
  if n == 1
    eigminB = B;
    eigmaxB = B;
  else
    eigminB = eigs(B, 1, 'sm', 'tolerance', eigtol);
    eigmaxB = eigs(B, 1, 'lm', 'tolerance', eigtol);
  end

  eigmin =  eigminA + eigminB;
  eigmax = eigmaxA + eigmaxB;

  normind = 1;

  outstruct.m = m;
  outstruct.n = n;
  outstruct.condA = eigmaxA / eigminA;
  outstruct.condB = eigmaxB / eigminB;
  outstruct.condM = eigmax / eigmin;

  outstruct.err_sylv = norm(Xtrue - Xsylv, normind) / norm(Xtrue, normind);
  outstruct.err_dfpm = norm(Xtrue - Xdfpm, normind) / norm(Xtrue, normind);
  outstruct.err_dfpmr = norm(Xtrue - Xdfpmr, normind) / norm(Xtrue, normind);

  normsol = (norm(A, normind) + norm(B, normind)) * norm(Xsylv, normind) + norm(C, normind);
  outstruct.res_sylv = norm(A*Xsylv + Xsylv*B - C, normind) / normsol;
  outstruct.res_dfpm = norm(A*Xdfpm + Xdfpm*B - C, normind) / normsol;
  outstruct.res_dfpmr = norm(A*Xdfpmr + Xdfpmr*B - C, normind) / normsol;

  outstruct.iter_dfpm = iter_dfpm;
  outstruct.iter_dfpmr = iter_dfpmr;
  outstruct.time_sylv = tsylvester;
  outstruct.time_dfpm = tdfpm;
  outstruct.time_dfpmr = tdfpmr;

  fprintf(['\n----------------------------------\n',...
           '----------------------------------\n',...
           '----------------------------------\n',...
           'Problem size:        %5d x %5d\n',...
           'Condition number of A:    %.2e\n',...
           'Condition number of B:    %.2e\n',...
           'Condition number of M:    %.2e\n',...
           '----------------------------------\n',...
           'Forward error sylvester:  %.2e\n',...
           'Forward error DFPM:       %.2e\n',...
           'Forward error DFPM ref:   %.2e\n',...
           'Residual sylvester:       %.2e\n',...
           'Residual DFPM:            %.2e\n',...
           'Residual DFPM ref:        %.2e\n',...
           'Time sylvester (seconds): %.2e\n',...
           'Time DFPM (seconds):      %.2e\n',...
           '   # of iterations:       %8d\n',...
           'Time DFPM ref (seconds):  %.2e\n',...
           '   # of iterations:       %8d\n',...
           'time_dfpm/time_sylvester: %.2e\n'],...
          outstruct.m, outstruct.n,...
          outstruct.condA,...
          outstruct.condB,...
          outstruct.condM,...
          outstruct.err_sylv,...
          outstruct.err_dfpm,...
          outstruct.err_dfpmr,...
          outstruct.res_sylv,...
          outstruct.res_dfpm,...
          outstruct.res_dfpmr,...
          outstruct.time_sylv,...
          outstruct.time_dfpm,...
          outstruct.iter_dfpm,...
          outstruct.time_dfpmr,...
          outstruct.iter_dfpmr,...
          tdfpm/tsylvester);

  function [Xtrue, A, B, C] = generate_sylv_test_case(m, n, cond1, cond2)
    Xtrue = randn(m, n);
    A = randeig(m, 1, cond1);
    B = randeig(n, 1, cond2);
    C = A*Xtrue + Xtrue*B;
  end
end
