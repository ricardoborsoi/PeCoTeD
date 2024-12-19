%% Compare methods for the solution of generalized Sylvester equation.
% In thie experiment, the size of the matrices varies, while the condition
% number of the equation is fixed. The coefficients all have eigenvalues
% between 1 and maxeig, where maxeig is either 10, or 100,

% Set global input parameters.
maxeig = [1e1, 1e2];
p = 5;
tol = 8*eps();
maxit = 50000;

for k = 1:length(maxeig)

  rng(k);

  % Set local input parameters.
  condA = maxeig(k);
  condB = maxeig(k);
  sizes = [2, 5:5:100];
  nsizes = length(sizes);
  clear outstruct

  % Run experiments only if corresponding .mat file is not present.
  filename = sprintf('gensylv_%d_%d_%.2e', sizes(end), p, condA);
  matfile = sprintf('%s.mat', filename);
  if ~exist(matfile, 'file')
    for i = 1:nsizes
      m = sizes(i);
      n = sizes(i);
      titer = tic;
      outstruct(i) = test_gensylv_case_fun(m, n, p, tol, maxit, condA, condB);
      toc(titer)
    end
    save(matfile, 'outstruct');
  else
    load(matfile);
  end

  % Plot output.
  figure(1000+k)
  clf

  subplot(2,2,1) % Residual as condition number increases.
  semilogy(sizes, [outstruct.res_kron])
  hold on
  semilogy(sizes, [outstruct.res_dfpm])
  semilogy(sizes, [outstruct.res_dfpmw])
  semilogy(sizes, [outstruct.res_gbia])
  title('Relative residual')

  legend('kron', 'dfpm', 'dfpm-nonoptimal', 'gbia')

  subplot(2,2,2) % Residual as condition number increases.
  semilogy(sizes, [outstruct.err_kron])
  hold on
  semilogy(sizes, [outstruct.err_dfpm])
  semilogy(sizes, [outstruct.err_dfpmw])
  semilogy(sizes, [outstruct.err_gbia])
  title('Relative forward error')

  subplot(2,2,3) % Execution time as condition number increases.
  semilogy(sizes, [outstruct.tkron])
  hold on
  semilogy(sizes, [outstruct.tdfpm])
  semilogy(sizes, [outstruct.tdfpmw])
  semilogy(sizes, [outstruct.tgbia])
  title('Execution time')

  subplot(2,2,4) % Execution time as condition number increases.
  semilogy(sizes, [outstruct.iter_dfpm])
  hold on
  semilogy(sizes, [outstruct.iter_dfpmw])
  semilogy(sizes, [outstruct.iter_gbia])
  title('Iterations')

  % Save output to file.
  % Size(1) residuals(4) errors(4) condMu(1) timings(4) iterations(3)
  datfile = sprintf('%s.dat', filename);
  outfile = fopen(datfile, 'w');
  for i = 1:nsizes
    fprintf(outfile,...
            ['%6d ',...
             '%.2e %.2e %.2e %.2e ',...
             '%.2e %.2e %.2e %.2e ',...
             '%.2e ',...
             '%.2e %.2e %.2e %.2e ',...
             '%6d %6d %6d\n'],...
            outstruct(i).m,...
            outstruct(i).res_kron, outstruct(i).res_dfpm,...
            outstruct(i).res_dfpmw, outstruct(i).res_gbia,...
            outstruct(i).err_kron, outstruct(i).err_dfpm,...
            outstruct(i).err_dfpmw, outstruct(i).err_gbia,...
            outstruct(i).condM * eps(),...
            outstruct(i).tkron, outstruct(i).tdfpm,...
            outstruct(i).tdfpmw, outstruct(i).tgbia,...
            outstruct(i).iter_dfpm, outstruct(i).iter_dfpmw,...
            outstruct(i).iter_gbia);
  end
  fclose(outfile);
end
