# htcomp
Digital filtering for evaluating Hankel transform (HT) integrals.

### Description

Just an exercise in refactoring legacy Fortran code (Fortran IV) from the paper:
[Anderson W L 1979 Numerical integration of related HTs of orders 0 and 1 adaptive digital filtering; Geophysics 44 1287–1305.](https://library.seg.org/doi/abs/10.1190/1.1441007)

### Compiling

A [Fortran Package Manager](https://github.com/fortran-lang/fpm) manifest file is included, so that a test (driver program) can be compiled with FPM. For example:

```
fpm test --profile release --V
```
