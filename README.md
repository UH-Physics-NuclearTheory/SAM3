# Subensemble Acceptance Method 3.0

## Description

This repository contains material related to the **Subensemble Acceptance Method 3.0 (SAM-3.0)**.

SAM-3.0 is a procedure to evaluate the effect of exact global conservation laws on cumulants of observables measured in a subsystem or acceptance. The method expresses cumulants in the constrained ensemble in terms of joint unconstrained cumulants of the measured observables and the globally conserved charges.

Compared to earlier versions of SAM, SAM-3.0 works directly with the joint cumulant-generating function of the observables and the total conserved charges. This makes it applicable to arbitrary observables, including non-conserved quantities such as net protons, and to multiple simultaneously conserved charges, such as baryon number, electric charge, strangeness, or total energy.

## Reference

The primary reference for this repository, which documents the method, is

- R. Poberezhniuk, V. A. Kuznietsov, G. Pihan, V. Vovchenko,  
  *Subensemble Acceptance Method 3.0: General Corrections to Cumulants from Exact Conservation Constraints*,  
  [arXiv:XXXX.XXXXX [hep-ph]](https://arxiv.org/abs/XXXX.XXXXX)


## Method

Let

$$
\mathbf{X}=(X_1,\ldots,X_d)
$$

denote the observables of interest, and let

$$
\mathbf{B}=(B_1,\ldots,B_s)
$$

denote the globally conserved charges. The input to SAM-3.0 is the set of joint unconstrained, or “grand-canonical,” cumulants

$$
\kappa^{\rm gce}_{\mathbf{n},\mathbf{m}},
$$

where the multiindex $\mathbf{n}$ counts derivatives with respect to the observable sources and the multiindex $\mathbf{m}$ counts derivatives with respect to the conserved-charge sources.

The constrained cumulants

$$
\kappa^{\rm ce}_{\mathbf{n}}
$$

are obtained by solving the saddle-point equations for auxiliary charge-conjugate fields

$$
\boldsymbol{\lambda}(\mathbf{t})
$$

and substituting the solution into the constrained cumulant-generating function. The recursion is algebraic and is implemented using multivariate Faà di Bruno relations, equivalently colored set partitions or multivariate Bell polynomials.

The implementation assumes that the relevant unconstrained distributions are sufficiently sharply peaked for the leading saddle-point approximation to be valid. The input cumulants need not originate from an actual thermal grand-canonical ensemble; the “gce” label denotes the unconstrained reference distribution.

## Material

This repository contains two complementary implementations.

### Wolfram Mathematica notebook

* [**Mathematica/SAM-3.0.nb**](Mathematica/SAM-3.0.nb)

The notebook derives and evaluates SAM-3.0 cumulant formulas symbolically. It supports arbitrary numbers of observables and conserved charges, and computes cumulants recursively to user-specified order.

The notebook contains two independent leading-order implementations:

1. **Colored partitions**
   Uses the explicit colored set-partition form of the multivariate Faà di Bruno formula. This is the most transparent implementation and follows the analytic derivation closely.

2. **Coefficient matching**
   Constructs a truncated cumulant-generating function, substitutes a formal power series for the saddle fields (\lambda_a(\mathbf t)), and solves the saddle equations by matching powers of the observable sources. This provides an independent cross-check of the colored-partition recursion.

The colored-partition implementation is the main symbolic implementation.

### C++ implementation

* [**cpp/sam3.h**](cpp/sam3.h)
* [**cpp/test_sam3.cpp**](cpp/test_sam3.cpp)

The C++ code implements the leading-order SAM-3.0 recursion for numerical use. The input is a map of joint unconstrained cumulants,

```cpp
std::map<std::pair<std::vector<int>, std::vector<int>>, double>
```

where the first multiindex labels observable orders and the second multiindex labels conserved-charge orders.

The main interface is

```cpp
CEMap ComputeSAM3CanonicalCumulants(
    int observableDim,
    int chargeDim,
    int desiredOrder,
    const KappaMap& gceCumulants
);
```

The function returns a map

```cpp
std::map<std::vector<int>, double>
```

containing all constrained observable cumulants up to the requested total order.

The C++ implementation uses Eigen for inversion of the conserved-charge covariance matrix.

## Input convention

A GCE cumulant is stored as

```cpp
kappa[{n, m}] = value;
```

where

```cpp
n = {n1, ..., nd}
m = {m1, ..., ms}
```

correspond to

$$
\kappa^{\rm gce}_{(n_1,\ldots,n_d),(m_1,\ldots,m_s)}.
$$

For example, with two observables and three conserved charges,

```cpp
kappa[{{2,1}, {1,0,2}}] = value;
```

represents

$$
\kappa^{\rm gce}_{(2,1),(1,0,2)}.
$$

Missing cumulants are interpreted as zero.

## Example

```cpp
#include "sam3.h"

int main() {
    KappaMap kappa;

    // Example: one observable X and one conserved charge B.
    // A key is { observable orders, charge orders }, each a list of integers.
    kappa[{{2}, {0}}] = 10.0; // kappa_X2
    kappa[{{1}, {1}}] = 3.0;  // kappa_XB
    kappa[{{0}, {2}}] = 5.0;  // kappa_B2

    CEMap ce = ComputeSAM3CanonicalCumulants(
        1,      // observable dimension
        1,      // charge dimension
        2,      // maximum observable cumulant order
        kappa
    );

    double C2 = ce[{2}];

    return 0;
}
```

For Gaussian input, the second constrained cumulant reduces to the Schur-complement form

$$
\kappa^{\rm ce}_{2}=\kappa^{\rm gce}_{2,0}-\frac{\left(\kappa^{\rm gce}_{1,1}\right)^2}{\kappa^{\rm gce}_{0,2}} .
$$

A runnable version of this example is in [`cpp/example.cpp`](cpp/example.cpp), and
a worked example — a toy acceptance model with non-trivial cumulants at every
order — is in [`cpp/example_binomial.cpp`](cpp/example_binomial.cpp). See also
[`cpp/README.md`](cpp/README.md) for C++ build and usage details.

## Requirements


### Mathematica

The notebook requires Wolfram Mathematica. No external Mathematica packages are required unless stated inside the notebook.

### C++

The C++ implementation requires

* C++17-compatible compiler
* Eigen (the CMake build can fetch this automatically; see below)

It can be built either directly or with CMake.

**Direct compilation** (no build system needed; requires a local Eigen):

```bash
cd cpp
g++ -std=c++17 -O2 test_sam3.cpp -I. -I/path/to/eigen -o test_sam3
./test_sam3
```

**CMake** — resolves Eigen automatically: it uses an installed Eigen if one is
found (or an explicit `-DEIGEN3_INCLUDE_DIR=/path/to/eigen`), and otherwise
**downloads Eigen 3.4.0** at configure time, so no Eigen installation is needed:

```bash
cd cpp
cmake -B build            # uses a local Eigen if present, otherwise downloads it
cmake --build build
ctest --test-dir build    # runs the test suite (test_sam3)
```

Pass `-DSAM3_FETCH_EIGEN=OFF` to forbid the download and require a local Eigen.

## Repository structure

The repository layout is

```text
.
├── README.md
├── AUTHORS.md
├── LICENSE
├── Mathematica/
│   └── SAM-3.0.nb
└── cpp/
    ├── CMakeLists.txt
    ├── README.md
    ├── sam3.h
    ├── example.cpp
    ├── example_binomial.cpp
    └── test_sam3.cpp
```

The repository also includes optional cross-check material (the C++ and
Mathematica implementations are compared on identical random inputs):
`cpp/bridge_sam3.cpp`, `Mathematica/test_sam3.wl`, `Mathematica/crosscheck_cpp.wl`,
and the driver `scripts/crosscheck.sh`.

## Notes and limitations

The current C++ implementation is intended as a compact reference implementation of the leading-order SAM-3.0 recursion. It assumes that the expansion point is chosen such that the reference mean charges coincide with the fixed charges,

$$
\langle B_a\rangle_{\rm ref}=B_{0,a}.
$$

Equivalently, the saddle satisfies

$$
\boldsymbol\lambda(\mathbf 0)=0.
$$

If the fixed charges do not coincide with the reference means, the input cumulants should first be evaluated at the appropriate nonzero saddle point or exponential tilt.

The implementation currently focuses on the leading saddle-point contribution. Finite-size corrections, if included, should be documented separately.

## Background and previous versions

SAM-3.0 builds on earlier versions of the subensemble acceptance method:

- V. Vovchenko, O. Savchuk, R. Poberezhnyuk, M. I. Gorenstein, V. Koch,  
  *Connecting fluctuation measurements in heavy-ion collisions with the grand-canonical susceptibilities*,  
  [Phys. Lett. B 811, 135868 (2020)](https://doi.org/10.1016/j.physletb.2020.135868),  
  [[arXiv:2003.13905 [hep-ph]](https://arxiv.org/abs/2003.13905)]

- V. Vovchenko, R. Poberezhnyuk, V. Koch,  
  *Cumulants of multiple conserved charges and global conservation laws*,  
  [JHEP 10, 089 (2020)](https://doi.org/10.1007/JHEP10(2020)089),  
  [[arXiv:2007.03850 [hep-ph]](https://arxiv.org/abs/2007.03850)]

- V. Vovchenko,  
  *Correcting event-by-event fluctuations in heavy-ion collisions for exact global conservation laws with the generalized subensemble acceptance method*,  
  [Phys. Rev. C 105, 014903 (2022)](https://doi.org/10.1103/PhysRevC.105.014903),  
  [[arXiv:2106.13775 [hep-ph]](https://arxiv.org/abs/2106.13775)]

## Attribution

Publications using this repository should cite the SAM-3.0 paper.  
Publications using results or formulas from earlier versions of the method should also cite the corresponding SAM-1.0 or SAM-2.0 papers listed above.

## Authors

See [AUTHORS.md](AUTHORS.md).

## License

This repository is distributed under the MIT License. See [LICENSE](LICENSE) for details.

## Copyright

Copyright (C) 2026
Roman Poberezhniuk, Volodymyr A. Kuznietsov, Grégoire Pihan, Volodymyr Vovchenko
