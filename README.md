# Rsampl #

[![Build Status](https://travis-ci.org/Qadra/Rsampl.svg?branch=master)](https://travis-ci.org/Qadra/Rsampl)
[![codecov](https://codecov.io/gh/Qadra/Rsampl/branch/master/graph/badge.svg)](https://codecov.io/gh/Qadra/Rsampl)

Implementation of Weighted Random Sampling algorithms for R.

The current implementation of Weighted Random Sampling in R is an *n²*
algorithm, this repository currently contains two algorithms for _weighted
random sampling_:

### Binary solution ###

Uses a simple binary tree with

* Sampling time: _O(k*log(n))_
* Preprocessing time: _O(n)_

Can be accessed using `method='binary'`.

### RSTree Solution ###

A more complicated tree with with the following

* Sampling time: _O(k * log(log(n)))_ expected
* Preprocessing time: _O(n)_

Can be accessed using `method='rstree'`.

# Installation #

When the repository has been cloned you can use `R CMD INSTALL .` in the cloned
directory to install the package in your user **R** installation.

# Usage #

To load the package execute `library(Rsampl)`. Usage examples can be found
using e.g. `?Rsample::wrs_sample`. The following are some basic usage examples.

### Using without preprocessing ###

Weighted Random Sampling without preprocessing:
```R
library(Rsampl)

# The number of indices we sample from
n <- 100
# The number of indices we sample on each sample. Precondition: k <= n
k <- 10

# Indexes we sample from
I <- 1:n
# The weights on which elements are selected
W <- runif(n)
W <- W/sum(W) # Normalize input!

# Sampling using a binary tree
samples <- wrs(I, W, k)
# is equal to
samples <- wrs(I, W, k=k, method="binary")

# The same for the RSTree
samples <- wrs(I, W, k=k, method="rstree")
```

### Usage with preprocessing ###

This library supports preparing the data structure as a separate step from the
actual sampling. If many samples need to be drawn from the same data
distribution, this leads to performance improvements.

```R
library(Rsampl)

# The number of indices we sample from
n <- 100
# The number of indices we sample on each sample. Precondition: k <= n
k <- 10

# Indexes we sample from
I <- 1:n
# The weights on which elements are selected
W <- runif(n)
W <- W/sum(W) # Normalize input!

# Generate the data structure to sample from:
binary_structure <- wrs_preprocessing(I, W, k, method='binary')

# Draw samples using the structure
samples <- wrs(k=k, struct=binary_structure)

```

# Development #

The following packages are required: `devtools`, `roxygen2`.

Start **R** in the base of the package directory load `devtools` and run unit
tests using `test()` (or `devtools::test()`).

After writing code you can refresh the package using `R CMD INSTALL .` followed
by `reload()` from the `devtools` package.


Before committing packages be sure to run `roxygenize()` followed by `check()`.
Warnings may be treated as build errors by Travis.
The `roxygenize()` command creates manuals and the `NAMESPACE` file from
comments in the included **R** code.

To view code coverage use the package `covr` and any of the following:
```R
# For just percentages
package_coverage()

# For browsable webinterface with code coverage
report(package_coverage())
```

# FAQ #
