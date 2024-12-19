%% Compare methods for the solution of generalized Sylvester equation.
% In thie experiment, the size of the matrices varies, while the condition
% number of the equation is fixed. The coefficients all have eigenvalues
% between 1 and maxeig, where maxeig is either 10, or 100,

% Set global input parameters.
maxeig = [1e1, 1e2];
p = 5;
tol = 8*eps();
maxit = 10000;

for k = 1:length(maxeig)

  rng(k);

  % Set local input parameters.
  condA = maxeig(k);
  condB = maxeig(k);
  sizes = [250, 500, 750, 1000, 1250, 1500, 1750, 2000];
  nsizes = length(sizes);
  clear outstruct

  % Run experiments only if corresponding .mat file is not present.
  filename = sprintf('large_gensylv_%d_%d_%.2e', sizes(end), p, condA);
  matfile = sprintf('%s.mat', filename);
  if ~exist(matfile, 'file')
    for i = 1:nsizes
      m = sizes(i);
      n = sizes(i);
      titer = tic;
      outstruct(i) = test_large_gensylv_case_fun(m, n, p, tol,...
                                                 maxit, condA, condB);
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
  semilogy(sizes, [outstruct.res_dfpmw])
  title('Relative residual')

  legend('dfpm-nonoptimal')

  subplot(2,2,2) % Residual as condition number increases.
  semilogy(sizes, [outstruct.err_dfpmw])
  title('Relative forward error')

  subplot(2,2,3) % Execution time as condition number increases.
  semilogy(sizes, [outstruct.tdfpmw])
  title('Execution time')

  subplot(2,2,4) % Execution time as condition number increases.
  semilogy(sizes, [outstruct.iter_dfpmw])
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
            outstruct(i).m,... % outstruct(i).res_dfpmw,...
            outstruct(i).err_dfpmw, outstruct(i).tdfpmw,...
            outstruct(i).iter_dfpmw);
  end
  fclose(outfile);
end