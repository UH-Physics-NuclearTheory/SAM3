// Self-contained tests for the SAM-3.0 C++ reference implementation (sam3.h).
//
//   g++ -std=c++17 -O2 test_sam3.cpp -I. -I/path/to/eigen -o test_sam3
//   ./test_sam3
//
// All checks compare against analytic ground truth:
//   - Gaussian input            -> Schur complement,
//   - Poisson split             -> exact binomial cumulants,
//   - two observables, one charge,
//   - multiple charges          -> matrix Schur complement,
//   - random input              -> the 2nd cumulant is always the Schur
//                                  complement (higher-order GCE cumulants do
//                                  not enter the conditioned variance).

#include "sam3.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

static int g_failures = 0;

static bool approx(double x, double y, double relTol = 1e-9, double absTol = 1e-9) {
    double diff = std::abs(x - y);
    double scale = std::max({1.0, std::abs(x), std::abs(y)});
    return diff <= absTol + relTol * scale;
}

static void requireApprox(const std::string& name, double got, double expected) {
    if (approx(got, expected)) {
        std::cout << "  PASS  " << name << " = " << std::setprecision(12) << got << "\n";
    } else {
        std::cout << "  FAIL  " << name << "  got=" << std::setprecision(17) << got
                  << "  expected=" << expected << "  diff=" << (got - expected) << "\n";
        ++g_failures;
    }
}

static double getCE(const CEMap& ce, const MI& n) {
    auto it = ce.find(n);
    if (it == ce.end()) { std::cout << "  FAIL  missing CE cumulant\n"; ++g_failures; return 0.0; }
    return it->second;
}

// --------------------------------------------------------------------------
// Test 1: one observable, one charge, Gaussian input. C2 = Var(X) - kXB^2/kB2.
// --------------------------------------------------------------------------
static void testGaussianOneCharge() {
    std::cout << "\n=== one observable, one charge, Gaussian ===\n";
    KappaMap kappa;
    kappa[{MI{2}, MI{0}}] = 10.0; // Var(X)
    kappa[{MI{1}, MI{1}}] = 3.0;  // Cov(X, B)
    kappa[{MI{0}, MI{2}}] = 5.0;  // Var(B)

    CEMap ce = ComputeSAM3CanonicalCumulants(1, 1, 4, kappa);
    requireApprox("Gaussian C1", getCE(ce, MI{1}), 0.0);
    requireApprox("Gaussian C2", getCE(ce, MI{2}), 8.2); // 10 - 9/5
    requireApprox("Gaussian C3", getCE(ce, MI{3}), 0.0);
    requireApprox("Gaussian C4", getCE(ce, MI{4}), 0.0);
}

// --------------------------------------------------------------------------
// Test 2: Poisson split X ~ Poi(a mu), Y ~ Poi((1-a) mu), B = X + Y, fixed
// B = <B>. The conditional X is Binomial(mu, a) -- exact, non-Gaussian.
// --------------------------------------------------------------------------
static void testPoissonBinomial() {
    std::cout << "\n=== Poisson split / exact binomial ===\n";
    const int Nmax = 4;
    const double mu = 100.0, alpha = 0.30;
    KappaMap kappa;
    for (int n = 0; n <= Nmax; ++n)
        for (int m = 0; m <= Nmax + 1; ++m) {
            if (n == 0 && m == 0) continue;
            kappa[{MI{n}, MI{m}}] = (n > 0) ? alpha * mu : mu;
        }
    CEMap ce = ComputeSAM3CanonicalCumulants(1, 1, Nmax, kappa);
    const double C1 = mu * alpha;
    const double C2 = mu * alpha * (1 - alpha);
    const double C3 = mu * alpha * (1 - alpha) * (1 - 2 * alpha);
    const double C4 = mu * alpha * (1 - alpha) * (1 - 6 * alpha * (1 - alpha));
    requireApprox("Binomial C1", getCE(ce, MI{1}), C1);
    requireApprox("Binomial C2", getCE(ce, MI{2}), C2);
    requireApprox("Binomial C3", getCE(ce, MI{3}), C3);
    requireApprox("Binomial C4", getCE(ce, MI{4}), C4);
}

// --------------------------------------------------------------------------
// Test 3: two observables, one charge, Gaussian covariance projection.
// --------------------------------------------------------------------------
static void testTwoObservables() {
    std::cout << "\n=== two observables, one charge, Gaussian ===\n";
    KappaMap kappa;
    kappa[{MI{2, 0}, MI{0}}] = 10.0;
    kappa[{MI{0, 2}, MI{0}}] = 20.0;
    kappa[{MI{1, 1}, MI{0}}] = 4.0;
    kappa[{MI{1, 0}, MI{1}}] = 3.0;
    kappa[{MI{0, 1}, MI{1}}] = -2.0;
    kappa[{MI{0, 0}, MI{2}}] = 5.0;

    CEMap ce = ComputeSAM3CanonicalCumulants(2, 1, 3, kappa);
    requireApprox("TwoObs C20", getCE(ce, MI{2, 0}), 10.0 - 3.0 * 3.0 / 5.0);
    requireApprox("TwoObs C02", getCE(ce, MI{0, 2}), 20.0 - 4.0 / 5.0);
    requireApprox("TwoObs C11", getCE(ce, MI{1, 1}), 4.0 + 6.0 / 5.0);
    requireApprox("TwoObs C30", getCE(ce, MI{3, 0}), 0.0);
    requireApprox("TwoObs C03", getCE(ce, MI{0, 3}), 0.0);
}

// --------------------------------------------------------------------------
// Helpers for random / multi-charge tests.
// --------------------------------------------------------------------------
static void genWeak(int dim, int total, MI& cur, std::vector<MI>& out) {
    if (dim == 1) { cur.push_back(total); out.push_back(cur); cur.pop_back(); return; }
    for (int k = 0; k <= total; ++k) { cur.push_back(k); genWeak(dim - 1, total - k, cur, out); cur.pop_back(); }
}
static std::vector<MI> weakComps(int dim, int total) {
    std::vector<MI> out; MI cur; genWeak(dim, total, cur, out); return out;
}
static Eigen::MatrixXd randomSPD(int s, std::mt19937& rng) {
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    Eigen::MatrixXd A(s, s);
    for (int i = 0; i < s; ++i) for (int j = 0; j < s; ++j) A(i, j) = dist(rng);
    Eigen::MatrixXd K = A * A.transpose();
    for (int i = 0; i < s; ++i) K(i, i) += 1.0 + s;
    return K;
}
static KappaMap makeRandomKappa(int d, int s, int obsMax, int chgMax, std::mt19937& rng) {
    std::uniform_real_distribution<double> dist(-2.0, 2.0);
    KappaMap kappa;
    for (int o = 0; o <= obsMax; ++o)
        for (const MI& obs : weakComps(d, o))
            for (int c = 0; c <= chgMax; ++c)
                for (const MI& chg : weakComps(s, c)) {
                    if (o == 0 && c == 0) continue;
                    kappa[{obs, chg}] = dist(rng);
                }
    Eigen::MatrixXd K = randomSPD(s, rng);
    MI zeroObs(d, 0);
    for (int a = 0; a < s; ++a) for (int b = 0; b < s; ++b) {
        MI chg(s, 0); chg[a] += 1; chg[b] += 1; kappa[{zeroObs, chg}] = K(a, b);
    }
    return kappa;
}
// Exact conditioned 2nd cumulant of a single observable: Var(X) - C^T K^{-1} C.
static double schurC2(const KappaMap& kappa, int s) {
    auto get = [&](const MI& n, const MI& m) {
        auto it = kappa.find({n, m}); return it == kappa.end() ? 0.0 : it->second; };
    double varX = get(MI{2}, MI(s, 0));
    Eigen::VectorXd C(s);
    Eigen::MatrixXd K(s, s);
    for (int a = 0; a < s; ++a) {
        MI ea(s, 0); ea[a] = 1;
        C[a] = get(MI{1}, ea);
        for (int b = 0; b < s; ++b) { MI cab(s, 0); cab[a] += 1; cab[b] += 1; K(a, b) = get(MI{0}, cab); }
    }
    return varX - C.dot(K.inverse() * C);
}

// --------------------------------------------------------------------------
// Test 4: multiple charges, Gaussian -> matrix Schur complement (order 2).
// --------------------------------------------------------------------------
static void testMultiChargeGaussian() {
    std::cout << "\n=== multiple charges, Gaussian (Schur complement) ===\n";
    for (int s = 1; s <= 3; ++s) {
        std::mt19937 rng(2024u + s);
        // Gaussian: only second-order cumulants populated.
        KappaMap kappa;
        std::uniform_real_distribution<double> dist(-2.0, 2.0);
        kappa[{MI{2}, MI(s, 0)}] = std::abs(dist(rng)) + 1.0;
        for (int a = 0; a < s; ++a) { MI ea(s, 0); ea[a] = 1; kappa[{MI{1}, ea}] = dist(rng); }
        Eigen::MatrixXd K = randomSPD(s, rng);
        MI zeroObs(1, 0);
        for (int a = 0; a < s; ++a) for (int b = 0; b < s; ++b) {
            MI chg(s, 0); chg[a] += 1; chg[b] += 1; kappa[{zeroObs, chg}] = K(a, b);
        }
        CEMap ce = ComputeSAM3CanonicalCumulants(1, s, 2, kappa);
        requireApprox("MultiCharge s=" + std::to_string(s) + " C2",
                      getCE(ce, MI{2}), schurC2(kappa, s));
    }
}

// --------------------------------------------------------------------------
// Test 5: random input -> the conditioned 2nd cumulant is the Schur complement
// for *any* Nmax (independent of higher-order GCE cumulants and charge means).
// --------------------------------------------------------------------------
static void testRandomSchur() {
    std::cout << "\n=== random input: C2 == Schur complement (any order) ===\n";
    for (int Nmax = 2; Nmax <= 6; ++Nmax)
        for (int s = 1; s <= 3; ++s) {
            std::mt19937 rng(700u + 10u * Nmax + s);
            KappaMap kappa = makeRandomKappa(1, s, Nmax, Nmax + 1, rng);
            CEMap ce = ComputeSAM3CanonicalCumulants(1, s, Nmax, kappa);
            requireApprox("Random s=" + std::to_string(s) + " Nmax=" + std::to_string(Nmax) + " C2",
                          getCE(ce, MI{2}), schurC2(kappa, s));
        }
}

int main() {
    try {
        testGaussianOneCharge();
        testPoissonBinomial();
        testTwoObservables();
        testMultiChargeGaussian();
        testRandomSchur();
    } catch (const std::exception& e) {
        std::cerr << "\nException: " << e.what() << "\n";
        return 1;
    }
    std::cout << "\n";
    if (g_failures == 0) { std::cout << "All SAM-3.0 tests passed.\n"; return 0; }
    std::cout << g_failures << " test(s) FAILED.\n";
    return 1;
}
