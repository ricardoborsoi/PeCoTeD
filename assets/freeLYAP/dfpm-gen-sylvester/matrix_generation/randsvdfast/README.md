[![View randsvdfast on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://mathworks.com/matlabcentral/fileexchange/74485)

# The MATLAB `randsvdfast` test matrix

This repository contains a MATLAB function to generates a matrix with specified singular values or 2-norm condition number and the corresponding unit tests.

The function is named after the `randsvd` matrix in the MATLAB `gallery`, as it provides similar functionalities but uses a faster algorithm. The method was designed to generate test matrices for extreme-scale benchmarks such as the [High-performance Linpack Benchmark (HPL)](https://www.netlib.org/benchmark/hpl/) or the [HPL-AI Mixed-Precision Benchmark](hpl-ai.org).

## Usage

The command
```
A = randsvdfast(n, kappa, mode, method, matrix, classname, realout)
```
generates a matrix `A` of class `classname` with condition number `kappa` and singular values distributed according to `mode`. The function generates a matrix of order `n` if `n` is a positive integer, and of size `n(1)` by `n(2)` if `n` is a vector of length 2. By default, `n` and `kappa` are both set to 10.

The functions provides functionalities similar to those of the MATLAB function `galley('randsvd', ...)`. The most notable difference is that this routine allows the user to specify a custom distribution of the singular values (see below), but does not implement the reduction to banded form.

The singular values can have one of the following distributions:
   + `mode` = 0: one large singular value and one small singular value,
   + `mode` = 1: one large singular value,
   + `mode` = 2: one small singular value,
   + `mode` = 3 (default): geometric distribution,
   + `mode` = 4: arithmetic distribution,
   + `mode` = 5: random singular values with uniformly distributed magnitude,
   + `mode` = 6: the vector `kappa` contains the singular values.

The parameter `method` selects the algorithm that will be used to generate the test matrix. It can take any of the following values:
   + `method` = 1 (default): [Alg. 3.1, 1],
   + `method` = 2: [Alg. 3.2, 1],
   + `method` = 3: [Alg. 4.1, 1] (only `mode` = 0, 1, 2),
   + `method` = 4: [Alg. 4.2, 1] (only `mode` = 0, 1, 2).

This function is faster for  `method` = 3 or 4 than for `method` = 1 or 2.

The algorithm uses an orthogonal matrix *Q* that depends on the value of the parameter `matrix`, which can take the following values:
   + `matrix` = 0 (default): *Q* is a Haar distributed random unitary generated as the *Q* factor of the QR decomposition of the matrix `randn(n(1),n(2))`.
   + `matrix` = an integer from 1 to 7: *Q* is the matrix `gallery('orthog',n,matrix)`.
   + `matrix` is the function handle of a two-argument function that generates an `n(1)`-by-`n(2)` matrix with orthonormal columns.

The output matrix will be of class `classname` where `classname` is either
`'single'` or `'double'`. Constants are computed in double precision, whereas
operations at the scalar level are performed in precision `classname`. The
entries of `A` will be real if `realout` is `true`, and complex otherwise. By
default the function generates a real matrix of doubles.

## Tests

The class-based unit tests for the `randsvdfast` function can be ran with the command `test_run`.

## Anymatrix

The repository can be downloaded as a remote group into the extensible matrix collection [Anymatrix](https://github.com/mmikaitis/anymatrix) with
```
anymatrix('g', 'randsvdfast', 'mfasi/randsvdfast-matlab')
```
and the test matrix can be generated with
```
anymatrix('randsvdfast/randsvdfast', n, kappa, mode, method, matrix, classname, realout)
```

## Reference

[1] M. Fasi & N. J. Higham. [Generating extreme-scale matrices with specified singular values or condition numbers](http://dx.doi.org/10.1137/20M1327938). SIAM J. Sci. Comput., 43(1), 663–684, 2021.

## License

The code is distributed under the terms of the 2-Clause BSD License, see [license.txt](license.txt)
