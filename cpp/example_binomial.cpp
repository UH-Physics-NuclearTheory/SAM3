// Worked SAM-3.0 example: a toy acceptance model (cumulants non-trivial at every
// order). See example.cpp for the minimal quick-start.
//
//   g++ -std=c++17 -O2 example_binomial.cpp -I. -I/path/to/eigen -o example_binomial
//   ./example_binomial
//
// A number of particles is produced, N ~ Poisson(mu), and each is independently
// accepted with probability alpha. The total N is globally conserved; the
// observable X is the accepted count. SAM-3.0 returns the cumulants of X at fixed
// total N -- here exactly those of a Binomial(N, alpha).

#include "sam3.h"

#include <iostream>

int main() {
    const double mu    = 100.0; // mean number of produced particles
    const double alpha = 0.3;   // single-particle acceptance probability
    const int    Nmax  = 4;     // highest cumulant order to compute

    // Joint unconstrained ("grand-canonical") cumulants of (X, N). For two
    // independent Poisson sources (accepted + rejected particles) they are
    //   kappa^gce_{n,m} = alpha*mu   if n >= 1   (n derivatives w.r.t. X),
    //                   = mu         if n == 0, m >= 1.
    KappaMap kappa;
    for (int n = 0; n <= Nmax; ++n)
        for (int m = 0; m <= Nmax + 1; ++m)
            if (n + m > 0)
                kappa[{{n}, {m}}] = (n > 0) ? alpha * mu : mu;

    // Cumulants of the accepted count X under exact conservation of the total N.
    CEMap ce = ComputeSAM3CanonicalCumulants(1, 1, Nmax, kappa);

    std::cout << "Accepted-count cumulants under exact conservation of the total:\n";
    for (int n = 1; n <= Nmax; ++n)
        std::cout << "  kappa^ce_" << n << " = " << ce[{n}] << "\n";

    // These are the cumulants of Binomial(mu, alpha):
    //   C1 = mu*a, C2 = mu*a(1-a), C3 = mu*a(1-a)(1-2a), C4 = mu*a(1-a)(1-6a(1-a))
    std::cout << "\nexpected (Binomial cumulants): 30, 21, 8.4, -5.46\n";
    return 0;
}
