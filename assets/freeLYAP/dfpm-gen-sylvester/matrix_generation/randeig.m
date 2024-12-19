function A = randeig(n, p, kappa)
% RANDEIG    Matrix with specified condition number and positive eigenvalues.
%    A = RANDEIG(N,P,KAPPA) generates a random N-by-N-P tensor A whose slices
%    have eigenvalues between 1 and KAPPA and are simultaneouly diagonalizable.

  assert(isreal(n) && isscalar(n) && round(n) == n);
  assert(isreal(p) && isscalar(p) && round(p) == p);
  assert(isreal(kappa) && isscalar(kappa) && kappa >= 1);

  A = zeros(n,n,p);
  R = randsvdfast(n, 2, 0, 3); % gallery('randsvd', n, 1);
  for i = 1:p
    D = rand(1, n) * (sqrt(kappa) - 1/sqrt(kappa)) + 1/sqrt(kappa); % eigenvalues
    [~, j] = max(D);
    D(j) = sqrt(kappa);
    [~, j] = min(D);
    D(j) = 1/sqrt(kappa);
    A(:,:,i) = (R .* D) / R;
  end
end
