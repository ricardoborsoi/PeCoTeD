function A = rsvd_psd(n, kappa)
% RSVD_PSD    Hermitian metrix with specified condition number.
%    A = RSVD_PSD(N,KAPPA) generates a random N-by-N positive definite matrix A.
%    If KAPPA is positive, then A is a randsvd matrix with 2-norm condition
%    number KAPPA. If KAPPA is negative, then A is a diagonally dominant matrix
%    with random normally distributed non-diagonal random entries.
  assert(isreal(n) && isscalar(n));
  if nargin < 2
    kappa = -1;
  else
    assert(isreal(kappa) && isscalar(kappa));
  end

  if kappa > 1
    A = gallery('randsvd', n, -kappa);
  else
    A = randn(n,n);
    A = A + A' + n*eye(n);
  end
end
