# Dynamical Functional Particle Method for the Generalized Sylvester Equation



### About

The MATLAB code in this repository was used to produce the figures and tables in [sect. 4, 1]. The test suite was developed in MATLAB 2020b, and the experiments were run in MATLAB 2021b on a machine equipped with an AMD Ryzen 7 PRO 5850U and 32 GiB of RAM.



### Dependencies

The code in this repository requires the MATLAB `randsvdfast` function from the [randsvdfast-matlab](https://github.com/mfasi/randsvdfast-matlab) GitHub repository. This is included as a git submodule for convenience and can be cloned with
```
git submodule init
```
from the command line.



### Running the tests

The experiments can be reproduced by running the MATLAB script `run_tests` which executes the following scripts in the `tests` folder:
* `test_sylv`, which produces the data used in [Figure 1, 1];
* `test_gensylv`, which produces the data used in [Figure 2, 1]; and
* `test_large_gensylv`, which produces the data used in [Table 2, 1].
The actual tests are performed by the functions with the `_case_fun` suffix, also in the `tests` folder, which generate the matrix equations for a given choice of parameters, solve them, and store the output in an array of `struct`s.

The algorithms we consider are implemented by the MATLAB functions in the `algorithms` folder:
* `solve_cgkr`, which solves a generalized Sylvester equation with Hermitian positive definite coefficients using the unpreconditioned conjugate gradient method on the Kronecker form of the equation;
* `solve_gen_dfpm`, which solves a generalized Sylvester equation using the dynamical functional particle method;
* `solve_gen_gbia`, which solves a generalized Sylvester equation using the gradient based iterative algorithm of Ding and Chen [2];
* `solve_gen_kron`, which solves a generalized Sylvester equation using MATLAB backslash operator on the Kronecker form of the matrix equation;
* `solve_sylv_dfpm`, which solves a Sylvester equation using the dynamical functional particle method (this is a stripped-down variant of `solve_gen_dfpm`); and
* `solve_sylv_dfpm_refinement`, which solves a Sylvester equation by using the dynamical functional particle method with the solution of the Bartels–Stewart algorithm as starting value for the iteration (this method is not reported as it is not competitive on our test test).

The `matrix_generation` folder contains utility functions for generating test matrices with properties of interest.



### References

[1] A. Dmytryshyn, M. Fasi, and M. Guliksson. [*The dynamical functional particle method for the generalized Sylvester equation*](http://eprints.maths.manchester.ac.uk/2804/). MIMS EPrint 2021.4, Manchester Institute for Mathematical Sciences, The University of Manchester, UK, February 2021. Revised December 2021.

[2] F. Ding and T. Chen. [*Gradient based iterative algorithms for solving a class of matrix equations*](https://doi.org/10.1109/tac.2005.852558). IEEE Transactions on Automatic Control, 50(8):1216-1221, 2005.

If you use the code in this repository, please cite the preprint [1].


### License

This software is distributed under the terms of the Expat License (see [LICENSE.md](./LICENSE.md)).
