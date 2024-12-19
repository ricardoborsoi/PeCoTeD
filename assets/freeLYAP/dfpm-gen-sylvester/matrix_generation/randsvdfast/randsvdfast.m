function A = randsvdfast(n, kappa, mode, method, matrix, classname, realout)
%RANDSVDFAST   Random matrix with pre-assigned singular values.
%   A = RANDSVDFAST(N,KAPPA,MODE,METHOD,MATRIX,CLASSNAME,REALOUT) generates a
%   matrix A of class CLASSNAME with condition number KAPPA and singular values
%   distributed according to MODE. The function generates a matrix of order N if
%   N is a positive integer, and of size N(1) by N(2) if N is a vector of length
%   2. By default, N and KAPPA are both set to 10.
%
%   The functions provides functionalities similar to those of the MATLAB
%   function GALLERY('RANDSVD', ...). The most notable difference is that this
%   routine allows the user to specify a custom distribution of the singular
%   values (see below), but does not implement the reduction to banded form.
%
%   The singular values can have one of the following distributions:
%      MODE = 0: one large singular value and one small singular value,
%      MODE = 1: one large singular value,
%      MODE = 2: one small singular value,
%      MODE = 3 (default): geometric distribution,
%      MODE = 4: arithmetic distribution,
%      MODE = 5: random singular values with uniformly distributed magnitude,
%      MODE = 6: the vector KAPPA contains the singular values.
%
%   The parameter METHOD selects the algorithm that will be used to generate the
%   test matrix. It can take any of the following values:
%      METHOD = 1 (default): [Alg. 3.1, 1],
%      METHOD = 2: [Alg. 3.2, 1],
%      METHOD = 3: [Alg. 4.1, 1] (only MODE = 0, 1, 2),
%      METHOD = 4: [Alg. 4.2, 1] (only MODE = 0, 1, 2).
%   This function is faster for METHOD = 3 or 4 than METHOD = 1 or 2.
%
%   The algorithm uses an orthogonal matrix Q that depends on the value of the
%   parameter MATRIX, which can take the following values:
%      MATRIX = 0 (default): Q is a Haar distributed random unitary
%               generated as the Q factor of the QR decomposition of the
%               matrix RANDN(N(1),N(2)).
%      MATRIX = an integer from 1 to 7: Q is the matrix
%               GALLERY('ORTHOG',N,MATRIX).
%      MATRIX is the function handle of a two-argument function that
%             generates an N(1)-by-N(2) matrix with orthonormal columns.
%
%   The output matrix will be of class CLASSNAME where CLASSNAME is either
%   'single' or 'double'. Constants are computed in double precision, whereas
%   operations at the scalar level are performed in precision CLASSNAME. The
%   entries of A will be real if REALOUT is TRUE, and complex otherwise. By
%   default the function generates a real matrix of doubles.
%
%   Reference:
%   [1] M. Fasi and N. J. Higham. Generating extreme-scale matrices with
%       specified singular values or condition numbers. SIAM J. Sci. Comput.,
%       43(1), 663–684, 2021.

% Copyright (c) 2020-2021 Massimiliano Fasi and Nicholas J. Higham
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%   * Redistributions of source code must retain the above copyright notice,
%     this list of conditions and the following disclaimer.
%
%   * Redistributions in binary form must reproduce the above copyright notice,
%     this list of conditions and the following disclaimer in the documentation
%     and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
% EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  %% Parse input and set default values.
  if nargin < 1
    n = 10;
  end
  if isempty(n) || length(n) > 2 || ~isnumeric(n) || ...
        any(~isfinite(n)) || any(round(n) ~= n) || any(n <= 0)
    error('randsvd:invalidN', ...
          'N must be a 1- or 2-dimensional array of positive integers.');
  end
  taketranspose = false;
  if length(n) == 2
    tmp = n;
    if n(1) >= n(2)
      m = tmp(1);
      n = tmp(2);
    else
      n = tmp(1);
      m = tmp(2);
      taketranspose = true;
    end
  else
    m = n;
  end


  if nargin < 2
    kappa = 10;
  end
  p = min(m, n);
  if isempty(kappa) || (length(kappa) ~= 1 && length(kappa) ~= p) || ...
        ~isnumeric(kappa) || any(~isfinite(kappa)) || any(kappa < 0) || ...
        (length(kappa) == 1 && kappa < 1)
    error('randsvd:invalidKAPPA', ...
          'KAPPA must be a valid condition number or set of singular values.');
  end

  if nargin < 3
    mode = 3;
  else
    if isempty(mode) || length(mode) ~= 1 || ~isnumeric(mode) || ...
          any(round(mode) ~= mode) || mode < 0 || mode > 6
      error('randsvd:invalidMODE', ...
            'MODE must be an integer between 0 and 6 inclusive.');
    end
  end

  if mode == 6 && length(kappa) ~= p
    error('randsvd:invalidKappaMode', ...
          'KAPPA must be a valid set of singular values when MODE is 6.');
  end

  if nargin < 4
    method = 1;
  else
    if isempty(method) || length(method) ~= 1 || ~isnumeric(method) || ...
          round(method) ~= method || method < 1 || method > 4
      error('randsvd:invalidMETHOD', ...
            'METHOD must be an integer between 1 and 4 inclusive.');
    end
  end

  if method == 3 || method == 4
    if length(kappa) ~= 1
      error('randsvd:invalidKappaMethod', ...
            'KAPPA must be a scalar at least 1 for this choice of METHOD.');
    end
    if mode > 2
      error('randsvd:invalidModeMethod', ...
            'MODE must be 0, 1, or 2 for this choice of METHOD.');
    end
    if m ~= n
      error('randsvd:invalidNMethod', ...
            'N(1) and N(2) must be equal for this choice of METHOD.');
    end
  end

  if nargin < 5
      matrix = 0;
    else
      if isempty(matrix) || ~isa(matrix, 'function_handle') && ...
            (length(matrix) ~= 1 || round(matrix) ~= matrix || ...
             matrix < 0 || matrix > 7)
        error('randsvd:invalidMATRIX', ...
              'MATRIX must be an integer between 0 and 7 inclusive.');
      end
  end

  if nargin < 6
    classname = 'double';
  else
    if ~ischar(classname) || (~strcmp(classname, 'single') && ...
                              ~strcmp(classname, 'double'))
      error('randsvd:invalidCLASSNAME', ...
            'CLASSNAME must be either ''single'' or ''double''.');
    end
  end

  if nargin < 7
    realout = true;
  else
    if ~islogical(realout) || length(realout) ~= 1
      error('randsvd:invalidREALOUT', ...
            'REALOUT must be a boolean value.');
    end
  end

  if isa(matrix, 'function_handle')
    orthog = matrix;
  else
    if matrix == 0
      orthog = @(m,n)orthoghaar(m,n,realout,classname);
    else
      orthog = @(m,n)orthogmatlab(m,n,classname);
    end
  end

  %% Generate and return test matrix.
  if method < 3
    % Generate singular values.
    switch mode
      case 0
        sigma = [1; ones(p-2, 1)./sqrt(kappa); 1./(kappa)];
      case 1
        sigma = [kappa; ones(p-1, 1)]./kappa;
      case 2
        sigma = [ones(p-1, 1); 1/kappa];
      case 3
        sigma = nthroot(kappa, 1-p).^(0:p-1)';
      case 4
        h = (kappa-1)/(p-1);
        sigma = fliplr(1 + h*(0:p-1))' / kappa;
      case 5
        sigma = exp(-rand(p, 1)*log(kappa));
      case 6
        sigma = reshape(kappa, [p, 1]);
    end
    % Call subfunction.
    if method == 1
      A = randsvd_fwd(m, n, orthog, sigma, realout, classname);
    else
      A = randsvd_bwd(m, n, orthog, sigma, realout, classname);
    end
    if taketranspose
      A = A.';
    end
  else
    % Call subfunction.
    if method == 3
      A = svdcond_fwd(n, kappa, orthog, mode, realout, classname);
    else
      A = svdcond_bwd(n, kappa, orthog, mode, realout, classname);
    end
  end

  %% Local functions.
  function Q = orthoghaar(m, n, realout, classname)
  % ORTHOGDEFAULT Haar distributed orthogonal or unitary matrix.
  %    A = ORTHOGDEFAULT(M,N,REALOUT,CLASSNAME) generates a random M-by-N matrix
  %    A distributed according to the Haar measure on the corresponding Stiefel
  %    manifold. The matrix has real entries if REALOUT=TRUE and complex entries
  %    if REALOUT=FALSE. The elements are of class CLASSNAME.
    A = randn(m, n, classname);
    if ~realout
      A = A + 1i*randn(m, n, classname);
    end
    [Q, R] = qr(A, 0);
    Q = Q .* sign(diag(R)');
  end

  function Q = orthogmatlab(m, n, classname)
  % ORTHOGMATLAB Interface for GALLERY('ORTHOG', ...).
    Q = gallery('orthog', m, matrix);
    Q = cast(Q(:, 1:n), classname);
  end

  function A = randsvd_fwd(m, n, orthog, sigma, realout, classname)
  % RANDSVD_FWD Random matrix with pre-assigned singular values.
    p = min(m, n);
    Q = orthog(m, p);
    Q = Q .* sigma';
    u = randn(p, 1, classname);
    v = randn(n-p, 1, classname);
    if realout
      alpha = -2;
    else
      theta = 2*(rand(1, 1, classname)-1/2) * pi;
      alpha = -(exp(1i * theta) + 1);
      u = u + 1i*randn(p, 1, classname);
      v = v + 1i*randn(n-p, 1, classname);
    end
    zeta = norm(u)^2 + norm(v)^2;
    alpha = alpha / zeta;
    y = Q * u;
    A = [Q + conj(alpha)*y*u' conj(alpha)*y*v'];
  end

  function A = randsvd_bwd(m, n, orthog, sigma, realout, classname)
  % RANDSVD_BWD Random matrix with pre-assigned singular values.
    p = min(m, n);
    Q = orthog(n, p) .* sigma;
    u = randn(p, 1, classname);
    v = randn(m-p, 1, classname);
    if realout
      alpha = -2;
    else
      theta = 2*(rand(1, 1, classname)-1/2) * pi;
      alpha = -(exp(1i * theta) + 1);
      u = u + 1i*randn(p, 1, classname);
      v = v + 1i*randn(m-p, 1, classname);
    end
    zeta = norm(u)^2 + norm(v)^2;
    alpha = alpha / zeta;
    y = Q * u;
    A = [Q' + alpha*u*y'; alpha*v*y'];
  end

  function A = svdcond_fwd(n, kappa, orthog, mode, realout, classname)
  % SVDCOND_FWD Random matrix with pre-assigned condition number.
    Q = orthog(n, n);
    if realout
      alpha = -2;
    else
      theta = 2*(rand(1, 1)-1/2) * pi;
      alpha = -(exp(1i * theta) + 1);
    end
    ell = floor((rand(1)*n))+1;
    qell = Q(ell, :);
    switch mode
      case 0
        s1 = sqrt(kappa);
        sn = 1/sqrt(kappa);
      case 1
        s1 = kappa;
        sn = 1;
      case 2
        s1 = 1;
        sn = 1/kappa;
    end
    if s1 ~= 1
      y = Q(:, 1)*Q(ell, 1)'*(s1-1);
    else
      y = zeros(n, 1, classname);
    end
    if sn ~= 1
      y = y + Q(:, n)*Q(ell, n)'*(sn-1);
    end
    y(ell) = y(ell) + 1;
    if s1 ~= 1
      Q(:, 1) = s1 * Q(:, 1);
    end
    if sn ~= 1
      Q(:, n) = sn * Q(:, n);
    end
    A = Q + alpha'*y*qell;
  end

  function A = svdcond_bwd(n, kappa, orthog, mode, realout, classname)
  % SVDCOND_BWD Random matrix with pre-assigned condition number.
    Q = orthog(n, n);
    if realout
      alpha = -2;
    else
      theta = 2*(rand(1, 1)-1/2) * pi;
      alpha = -(exp(1i * theta) + 1);
    end
    ell = floor((rand(1)*n))+1;
    qell = Q(:, ell);
    switch mode
      case 0
        s1 = sqrt(kappa);
        sn = 1/sqrt(kappa);
      case 1
        s1 = kappa;
        sn = 1;
      case 2
        s1 = 1;
        sn = 1/kappa;
    end
    if s1 ~= 1
      y = Q(1, :)*Q(1, ell)'*(s1-1);
    else
      y = zeros(1, n, classname);
    end
    if sn ~= 1
      y = y + Q(n, :)*Q(n, ell)'*(sn-1);
    end
    y(ell) = y(ell) + 1;
    if s1 ~= 1
      Q(1, :) = s1 * Q(1, :);
    end
    if sn ~= 1
      Q(n, :) = sn * Q(n, :);
    end
    A = Q + alpha*qell*y;
  end

end