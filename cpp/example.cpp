// Minimal quick-start for the SAM-3.0 C++ reference implementation (sam3.h).
//
//   g++ -std=c++17 -O2 example.cpp -I. -I/path/to/eigen -o example   # or use CMake
//   ./example
//
// Set the grand-canonical (unconstrained) cumulants by hand, get the constrained
// cumulants out -- shown for one conserved charge (B), then two (B and Q).

#include "sam3.h"

#include <iostream>

int main() {
    // --- One observable X, one conserved charge B --------------------------
    // A key is { observable orders, charge orders }; entries left unset are zero.
    // (Both lists are std::vector<int>, aliased MI in sam3.h, so {{2}, {0}} is
    //  the same as {MI{2}, MI{0}}.)
    KappaMap kappa;
    kappa[{{2}, {0}}] = 10.0;  // Var(X)     = kappa^gce_{2,0}
    kappa[{{1}, {1}}] =  3.0;  // Cov(X, B)  = kappa^gce_{1,1}
    kappa[{{0}, {2}}] =  5.0;  // Var(B)     = kappa^gce_{0,2}

    CEMap ce = ComputeSAM3CanonicalCumulants(/*observables=*/1, /*charges=*/1,
                                             /*max order=*/2, kappa);
    std::cout << "One conserved charge (B):\n";
    std::cout << "  kappa^ce_1 = " << ce[{1}] << "\n";   // 0   (no mean specified)
    std::cout << "  kappa^ce_2 = " << ce[{2}] << "\n";   // 8.2 = 10 - 3*3/5  (Schur complement)

    // Higher orders work the same way -- raise the order argument. For this
    // Gaussian input (only 2nd-order cumulants set) the higher ones vanish:
    CEMap ce4 = ComputeSAM3CanonicalCumulants(1, 1, /*max order=*/4, kappa);
    std::cout << "  kappa^ce_3 = " << ce4[{3}] << ",  kappa^ce_4 = " << ce4[{4}]
              << "   (zero: Gaussian input)\n";

    // --- One observable X, two conserved charges B and Q -------------------
    // The charge index now has two entries: { B-order, Q-order }.
    KappaMap kappaBQ;
    kappaBQ[{{2}, {0, 0}}] = 10.0;  // Var(X)
    kappaBQ[{{1}, {1, 0}}] =  2.0;  // Cov(X, B)
    kappaBQ[{{1}, {0, 1}}] =  1.0;  // Cov(X, Q)
    kappaBQ[{{0}, {2, 0}}] =  2.0;  // Var(B)
    kappaBQ[{{0}, {0, 2}}] =  2.0;  // Var(Q)
    kappaBQ[{{0}, {1, 1}}] =  1.0;  // Cov(B, Q)

    CEMap ceBQ = ComputeSAM3CanonicalCumulants(/*observables=*/1, /*charges=*/2,
                                               /*max order=*/2, kappaBQ);
    std::cout << "\nTwo conserved charges (B, Q):\n";
    std::cout << "  kappa^ce_2 = " << ceBQ[{2}] << "\n";
    return 0;
}
