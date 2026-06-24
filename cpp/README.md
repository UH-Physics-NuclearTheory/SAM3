# SAM-3.0 — C++ reference implementation

Header-only (`sam3.h`); the only dependency is a C++17 compiler and
[Eigen](https://eigen.tuxfamily.org).

For the method, the input convention, limitations, and references, see the
[top-level README](../README.md). This file just covers building and using the
C++ code.

## Files

| file | purpose |
|------|---------|
| `sam3.h`              | the implementation (single header) |
| `example.cpp`         | minimal quick-start (below): one and two conserved charges |
| `example_binomial.cpp`| worked example: a toy acceptance model, non-trivial at every order |
| `test_sam3.cpp`       | correctness tests against analytic ground truth |
| `bridge_sam3.cpp`     | emits data for the C++ ↔ Mathematica cross-check (`../scripts/crosscheck.sh`) |
| `CMakeLists.txt`      | optional CMake build |

## Quick start

```cpp
#include "sam3.h"

KappaMap kappa;                 // joint unconstrained cumulants; missing entries = 0
kappa[{{2}, {0}}] = 10.0;       // Var(X)
kappa[{{1}, {1}}] =  3.0;       // Cov(X, B)
kappa[{{0}, {2}}] =  5.0;       // Var(B)

// Constrained cumulants of X (2nd order), under exact global conservation of B:
CEMap ce = ComputeSAM3CanonicalCumulants(/*observables=*/1, /*charges=*/1,
                                         /*max order=*/2, kappa);
double kappaCE2 = ce[{2}];      // = 8.2  (Var(X) - Cov(X,B)^2 / Var(B))
```

A key is `{ observable-orders, charge-orders }`, each a `std::vector<int>`; e.g.
`kappa[{{2,1}, {1,0,2}}]` is the joint cumulant with two observables (orders 2, 1)
and three charges (orders 1, 0, 2). `sam3.h` also aliases `std::vector<int>` as
`MI`, so `{MI{2,1}, MI{1,0,2}}` is the same thing. See [`example.cpp`](example.cpp)
for a complete runnable program.

## Build and run

**Direct (g++)** — requires a local Eigen:

```bash
g++ -std=c++17 -O2 example.cpp -I. -I/path/to/eigen -o example
./example
```

**CMake** — finds a local Eigen, or downloads Eigen 3.4.0 automatically:

```bash
cmake -B build                 # or -DEIGEN3_INCLUDE_DIR=/path/to/eigen
cmake --build build
./build/example                # minimal quick-start
./build/example_binomial       # worked acceptance-model example
ctest --test-dir build         # the test suite (test_sam3)
```

Pass `-DSAM3_FETCH_EIGEN=OFF` to forbid the automatic download.
