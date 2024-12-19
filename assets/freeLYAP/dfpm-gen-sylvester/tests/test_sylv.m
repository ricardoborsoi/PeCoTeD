% With easy matrices and n = m / 10.
% filename = 'sylvester_diag_dom.mat';
% if ~exist(filename, 'file')
%   counter = 0;
%   for i = [1,2,3,4,5,6,7,8,9,10]
%     counter = counter + 1;
%     outstruct_diag_dom(counter) = test_sylv_case_fun(i*100, i*10, -1, -1);
%   end
%   save(filename, 'outstruct_diag_dom')
% else
%   load(filename);
% end

%% With harder matrices and n = m / 10.
conds = [10, 100];

figcounter = 0;
for condnum = conds

  rng(condnum);
  figcounter = figcounter + 1;
  figure(figcounter)
  clf

  condA = condnum;
  condB = condnum;
  tol = 8 * eps();
  maxit = 10000;

  ms = [2, 10:10:500];
  n_ms = length(ms);
  ns = 500;
  n_ns = length(ns);

  for j = 1:n_ns
    n = ns(j);
    clear outstruct
    filename = sprintf('sylvester_%d_%d_%d',...
                       round(log10(condA)), round(log10(condB)), n);
    matfile = sprintf('%s.%s', filename, 'mat');
    if ~exist(matfile, 'file')
      for i = 1:n_ms
        m = ms(i);
        outstruct(i,j) = test_sylv_case_fun(m, n, condA, condB, tol, maxit);
      end
      save(matfile, 'outstruct');
    else
      load(matfile);
    end
  end

  m_vec = reshape([outstruct.m], [n_ms, n_ns]);
  time_dfpm = reshape([outstruct.time_dfpm], [n_ms, n_ns]);
  time_dfpmr = reshape([outstruct.time_dfpmr], [n_ms, n_ns]);
  time_sylv = reshape([outstruct.time_sylv], [n_ms, n_ns]);

  err_dfpm = reshape([outstruct.err_dfpm], [n_ms, n_ns]);
  err_dfpmr = reshape([outstruct.err_dfpmr], [n_ms, n_ns]);
  err_sylv = reshape([outstruct.err_sylv], [n_ms, n_ns]);
  condM = reshape([outstruct.condM], [n_ms, n_ns]);

  res_dfpm = reshape([outstruct.res_dfpm], [n_ms, n_ns]);
  res_dfpmr = reshape([outstruct.res_dfpmr], [n_ms, n_ns]);
  res_sylv = reshape([outstruct.res_sylv], [n_ms, n_ns]);

  subplot(1,3,1);
  plot(m_vec, time_dfpm)
  hold on
  plot(m_vec, time_dfpmr)
  plot(m_vec, time_sylv)
  legend('dfpm', 'sylvester')
  title('Timing')

  subplot(1,3,2);
  semilogy(m_vec, err_dfpm)
  hold on
  semilogy(m_vec, err_dfpmr)
  plot(m_vec, err_sylv)
  semilogy(m_vec, condM * eps())
  legend('dfpm', 'sylvester', '\kappa(M)u')
  title('Forward error')

  subplot(1,3,3);
  semilogy(m_vec, res_dfpm)
  hold on
  semilogy(m_vec, res_dfpmr)
  plot(m_vec, res_sylv)
  legend('dfpm', 'sylvester')
  title('Residual')

  datfile = sprintf('%s.%s', filename, 'dat');
  fid = fopen(datfile, 'w');
  for i = 1:length(time_dfpm)
    fprintf(fid, '%.4d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n',...
            m_vec(i),...
            res_dfpm(i), res_dfpmr(i), res_sylv(i),...
            err_dfpm(i), err_dfpmr(i),  err_sylv(i), condM(i) * eps(),...
            time_dfpm(i), time_dfpmr(i), time_sylv(i));
  end
  fclose(fid);
end
