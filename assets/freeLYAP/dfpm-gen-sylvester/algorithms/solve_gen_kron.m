function X = solve_gen_kron(A, B, C)
% SOLVE_GEN_KRON    Solve a generalized Sylvester equation.
%    X = SOLVE_GEN_KRON(A,B,C) finds the m-by-n matrix X such that
%
%       A(:,:,1)*X*B(:,:,1) + ... + A(:,:,p)*X*B(:,:,p) = C
%
%    where A are B are an m-by-m-by-p and n-by-n-p tensors, respectively.
%    This is done by solving the Kronecker linear system associted with the
%    matrix equation using the MATLAB backslash operator.

  [m, m1, ell] = size(A);
  [n, n1, ell1] = size(B);

  % Check that dimensions are consistent.
  assert(m == m1);
  assert(n == n1);
  assert(ell == ell1);

  % Generate coefficient matrix and solve linear system.
  M = zeros(m*n, m*n);
  for i = 1:ell
    M = M + kron(B(:, :, i)', A(:, :, i));
  end
  X = reshape(M \ C(:), [m,n]);
end