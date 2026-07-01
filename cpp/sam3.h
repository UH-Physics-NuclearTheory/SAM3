/*
 * sam3.h  --  Subensemble Acceptance Method 3.0 (SAM-3.0).
 *
 * Leading-order SAM-3.0 recursion: the canonical (globally charge-conserved)
 * cumulants of subsystem observables, reconstructed from joint unconstrained
 * ("grand-canonical") cumulants of the observables and the conserved charges.
 *
 * Self-contained: requires only a C++17 compiler and Eigen.
 *
 * --------------------------------------------------------------------------
 * Reference
 * --------------------------------------------------------------------------
 *   R. Poberezhniuk, V. A. Kuznietsov, G. Pihan, V. Vovchenko,
 *   "Subensemble Acceptance Method 3.0: General Corrections to Cumulants from
 *    Exact Conservation Constraints", arXiv:XXXX.XXXXX [hep-ph].
 *
 * Earlier versions of the subensemble acceptance method:
 *   - SAM, single conserved charge:      arXiv:2003.13905 (Phys. Lett. B 811, 135868)
 *   - SAM, multiple conserved charges:   arXiv:2007.03850 (JHEP 10 (2020) 089)
 *   - Generalized SAM:                   arXiv:2106.13775 (Phys. Rev. C 105, 014903)
 *
 * --------------------------------------------------------------------------
 * Method
 * --------------------------------------------------------------------------
 * For observables X = (X_1,...,X_d) correlated with globally conserved charges
 * B = (B_1,...,B_s), given the joint unconstrained cumulants kappa^gce_{n,m},
 * the constrained cumulants kappa^ce_n are obtained by:
 *   1. Solving order by order for the saddle-point coefficients lambda_{a;n}
 *      of the auxiliary charge-conjugate fields that enforce the conservation
 *      constraint, via inversion of the charge covariance matrix
 *      K_{ab} = kappa^gce_{0; e_a + e_b}.
 *   2. Substituting into the constrained cumulant-generating function, expressed
 *      as a sum over colored set partitions of the n observable derivatives
 *      (multivariate Faa di Bruno / multivariate Bell polynomials), with blocks
 *      weighted by kappa^gce and the lambda coefficients.
 *
 * This mirrors the "ColoredPartitions" method of the companion Mathematica
 * notebook (Mathematica/SAM-3.0.nb), whose independent "CoefficientMatching"
 * method provides a symbolic cross-check. Leading saddle-point order only; the
 * expansion point is assumed to satisfy lambda(0) = 0 (reference mean charges
 * equal the fixed charges). See README.md for details and limitations.
 *
 * --------------------------------------------------------------------------
 * Copyright (c) 2026
 *   Roman Poberezhniuk, Volodymyr A. Kuznietsov, Gregoire Pihan, Volodymyr Vovchenko
 * Distributed under the MIT License (see LICENSE).
 * --------------------------------------------------------------------------
 */

#ifndef SAM3_H
#define SAM3_H


#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <set>

// ---------------------------------------------------------------------------
// Input/output types.
//   MI       : a multi-index (vector of non-negative integer orders).
//   KappaMap : joint unconstrained ("GCE") cumulants kappa^gce_{n,m}, keyed by
//              (observable multi-index n of length d, charge multi-index m of
//              length s). Missing entries are interpreted as zero.
//   CEMap    : constrained ("CE") observable cumulants kappa^ce_n, keyed by n.
// ---------------------------------------------------------------------------
using MI = std::vector<int>;
using KappaKey = std::pair<MI, MI>;
using KappaMap = std::map<KappaKey, double>;
using CEMap    = std::map<MI, double>;

// ---------------------------------------------------------------------------
// Colored set partitions -- the multivariate Faa di Bruno building block (see
// the Appendix of arXiv:2106.13775). Bundled here so this header is
// dependency-free; same names/algorithm as the sample-moments library.
// ---------------------------------------------------------------------------
namespace sam3detail {

  using Block = std::vector<int>;
  using ColoredBlock = std::pair<int, Block>;   // (charge type 0..m-1, element indices)
  using ColoredPartition = std::vector<ColoredBlock>;

  /// All partitions of the set {0,...,n-1} into blocks, each block painted with
  /// one of m colors. Recursive; cached per (n, m) for the whole thread (the
  /// partitions never depend on the cumulant values, so the cache is shared
  /// across all calls -- the dominant saving in a parameter scan).
  inline const std::vector<ColoredPartition>&
  ColoredPartitionsOfSet(int n, int m) {
    static thread_local std::map<std::pair<int, int>,
                                 std::vector<ColoredPartition>> cache;
    auto key = std::make_pair(n, m);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    std::vector<ColoredPartition> ret;
    if (n == 0) return cache[key] = ret;       // empty set
    if (n == 1) {
      for (int j = 0; j < m; ++j)
        ret.push_back({ ColoredBlock(j, Block{0}) });
      return cache[key] = ret;
    }
    const auto& sub = ColoredPartitionsOfSet(n - 1, m);
    for (const auto& el : sub) {
      for (int j = 0; j < m; ++j) {            // element n-1 in a new block of type j
        auto t = el;
        t.push_back({ j, Block{n - 1} });
        std::sort(t.begin(), t.end());
        ret.push_back(t);
      }
      auto t = el;                             // element n-1 added to an existing block
      for (std::size_t i = 0; i < el.size(); ++i) {
        t[i].second.push_back(n - 1);
        ret.push_back(t);
        std::sort(t.begin(), t.end());
        t[i] = el[i];
      }
    }
    return cache[key] = ret;
  }

} // namespace sam3detail

/// SAM-3.0 canonical cumulant calculator (colored partition-of-set recursion).
///
/// Conventions (types defined at the top of this header):
///   - MI       : a multi-index, std::vector<int>.
///   - KappaMap : GCE input, keyed by (observable MI of length d, charge MI of
///                length s); value is the joint cumulant kappa^gce_{n,m}.
///                E.g. kappa[{ {1}, {1,0} }] = Cov(X, Q_1),
///                     kappa[{ {0}, {2,0,0} }] = Var(Q_1).
///   - CEMap    : output, keyed by an observable MI n; value is kappa^ce_n.
///
/// Use the free function ComputeSAM3CanonicalCumulants() rather than this class
/// directly. The charge covariance matrix K_{ab} = kappa^gce_{0; e_a+e_b} must
/// be non-singular (an exception is thrown otherwise).
class SAM3Calculator {
public:
    /// \param observableDim  number of subsystem observables d (>= 1)
    /// \param chargeDim      number of globally conserved charges s (>= 1)
    /// \param maxOrder       highest total cumulant order Nmax to compute
    /// \param gceCumulants   grand-canonical input cumulants kappa^gce (see above)
    SAM3Calculator(int observableDim,
                   int chargeDim,
                   int maxOrder,
                   const KappaMap& gceCumulants)
        : d_(observableDim),
          s_(chargeDim),
          Nmax_(maxOrder),
          kappa_(gceCumulants) {
        if (d_ <= 0 || s_ <= 0 || Nmax_ < 0) {
            throw std::invalid_argument("Invalid dimensions or maximum order.");
        }
        buildIndex();
    }

    /// Run the calculation: build K, solve the lambda coefficients, then
    /// evaluate kappa^ce_n for every observable multi-index n of total order
    /// 1 .. Nmax. Returns the CE cumulants keyed by n.
    CEMap compute() {
        buildKMatrix();             // K_{ab} = kappa^gce_{0; e_a+e_b}
        solveAllLambdaCoefficients(); // lambda_{a;n} up to order Nmax-1

        CEMap result;

        for(int nord = 1; nord <= Nmax_; ++nord) {
            const auto& ns = weakMI(d_, nord);
            for(auto& n : ns)
                result[n] = computeKappaCEn(n);
        }

        return result;
    }

    // Tracks missing cumulants during computation
    mutable std::set<KappaKey> missingKeys_;

    /// Prints all unique missing GCE cumulants that were treated as 0.
    /// Written to stderr so that results captured from stdout stay clean.
    void printMissingCumulants() const {
        if (missingKeys_.empty()) return;

        std::cerr << "\nMissing GCE cumulants (treated as 0):\n";
        for (const auto& key : missingKeys_) {
            std::cerr << "  kappa(" << MIToString(key.first) << "; "
                << MIToString(key.second) << ")\n";
        }
    }

private:
    int d_;
    int s_;
    int Nmax_;
    const KappaMap& kappa_;

    Eigen::MatrixXd K_;
    Eigen::FullPivLU<Eigen::MatrixXd> Klu_;
    Eigen::MatrixXd Kinv_;
    std::vector<double> B0_;

    // Integer encoding of multi-indices (mixed radix base B_): observable
    // indices -> [0, powObs_), charge indices -> [0, powChg_). This lets every
    // hot map be keyed by a long long instead of a std::vector, removing the
    // per-lookup key allocation and the lexicographic vector comparison.
    long long B_ = 0;
    long long powObs_ = 1;   // B_^d_
    long long powChg_ = 1;   // B_^s_
    std::vector<long long> powA_; // powA_[a] = B_^a: encoding offset of charge unit e_a

    // flattened GCE cumulants: key = encObs(obs)*powChg_ + encChg(chg)
    std::unordered_map<long long, double> kappaIdx_;

    // lambda_{a;n}:  key = a*powObs_ + encObs(n)
    std::unordered_map<long long, double> lambda_;
    // R_{a;n}:       key = a*powObs_ + encObs(n)
    std::unordered_map<long long, double> R_;

    // Brq cache:     key = encObs(r)*powChg_ + encChg(q)
    std::unordered_map<long long, double> Brq_;
    std::unordered_map<long long, char> Brqscomputed_; // key = encObs(r)

    // caches for multiindex generation
    std::map<std::pair<int, int>, std::vector<MI>> weakCache_;
    std::unordered_map<long long, std::vector<MI>> subCache_; // key = encObs(n)

    // partitions of a set: cached thread-locally inside getSetPoS (shared
    // across instances), since they depend only on (set size, #charge types).
    using PoSmulti = std::vector<sam3detail::ColoredPartition>;

    // reused scratch buffers so getqpi()/getkb()/buildElm() don't allocate a
    // vector per call (qpiScratch_ must stay distinct from kbScratch_: a getqpi
    // result can stay live across getkb calls).
    MI qpiScratch_;
    MI kbScratch_;
    std::vector<int> elmScratch_;

private:
    static MI zeroMI(int dim) {
        return MI(dim, 0);
    }

    static MI unitMI(int dim, int a) {
        MI e(dim, 0);
        e[a] = 1;
        return e;
    }

    static int total(const MI& n) {
        return std::accumulate(n.begin(), n.end(), 0);
    }

    static MI addMI(const MI& a, const MI& b) {
        MI c(a.size());
        for (std::size_t i = 0; i < a.size(); ++i) c[i] = a[i] + b[i];
        return c;
    }

    static MI subMIvec(const MI& a, const MI& b) {
        MI c(a.size());
        for (std::size_t i = 0; i < a.size(); ++i) c[i] = a[i] - b[i];
        return c;
    }

    static bool isZeroMI(const MI& n) {
        for (int x : n) {
            if (x != 0) return false;
        }
        return true;
    }

    static bool isUnitMI(const MI& n) {
        return total(n) == 1;
    }

    static double binomInt(int n, int k) {
        if (k < 0 || k > n) return 0.0;
        k = std::min(k, n - k);

        double res = 1.0;
        for (int i = 1; i <= k; ++i) {
            res *= static_cast<double>(n - k + i);
            res /= static_cast<double>(i);
        }
        return res;
    }

    static double factorialInt(int n) {
        double res = 1.0;
        for (int i = 2; i <= n; ++i) res *= static_cast<double>(i);
        return res;
    }

    static double miBinom(const MI& n, const MI& p) {
        double res = 1.0;
        for (std::size_t i = 0; i < n.size(); ++i) {
            res *= binomInt(n[i], p[i]);
        }
        return res;
    }

    static double miFact(const MI& n) {
        double res = 1.0;
        for (int x : n) res *= factorialInt(x);
        return res;
    }

    const std::vector<MI>& weakMI(int dim, int N) {
        auto key = std::make_pair(dim, N);
        auto it = weakCache_.find(key);
        if (it != weakCache_.end()) return it->second;

        std::vector<MI> out;

        if (dim == 1) {
            out.push_back(MI{N});
        } else {
            for (int k = 0; k <= N; ++k) {
                const auto& tails = weakMI(dim - 1, N - k);
                for (const MI& tail : tails) {
                    MI v;
                    v.reserve(dim);
                    v.push_back(k);
                    v.insert(v.end(), tail.begin(), tail.end());
                    out.push_back(v);
                }
            }
        }

        weakCache_[key] = out;
        return weakCache_[key];
    }

    const std::vector<MI>& subMI(const MI& n) {
        const long long key = encObs(n);
        auto it = subCache_.find(key);
        if (it != subCache_.end()) return it->second;

        std::vector<MI> out;
        MI cur(n.size(), 0);
        subMIRec(n, 0, cur, out);

        subCache_[key] = out;
        return subCache_[key];
    }

    static void subMIRec(const MI& n, int pos, MI& cur, std::vector<MI>& out) {
        if (pos == static_cast<int>(n.size())) {
            out.push_back(cur);
            return;
        }

        for (int k = 0; k <= n[pos]; ++k) {
            cur[pos] = k;
            subMIRec(n, pos + 1, cur, out);
        }
    }

    // Mixed-radix encodings (base B_) of observable / charge multi-indices.
    long long encObs(const MI& n) const {
        long long id = 0, p = 1;
        for (int i = 0; i < d_; ++i) { id += static_cast<long long>(n[i]) * p; p *= B_; }
        return id;
    }
    long long encChg(const MI& m) const {
        long long id = 0, p = 1;
        for (int i = 0; i < s_; ++i) { id += static_cast<long long>(m[i]) * p; p *= B_; }
        return id;
    }

    // Choose a radix larger than any component that can appear (input entries,
    // and queries which reach charge order Nmax_+1), then flatten kappa_ once.
    void buildIndex() {
        int maxComp = Nmax_ + 1;
        for (const auto& kv : kappa_) {
            for (int x : kv.first.first)  maxComp = std::max(maxComp, x);
            for (int x : kv.first.second) maxComp = std::max(maxComp, x);
        }
        B_ = maxComp + 1;
        powObs_ = 1; for (int i = 0; i < d_; ++i) powObs_ *= B_;
        powChg_ = 1; for (int i = 0; i < s_; ++i) powChg_ *= B_;
        powA_.resize(s_);
        for (int a = 0; a < s_; ++a) powA_[a] = (a == 0) ? 1 : powA_[a - 1] * B_;

        kappaIdx_.reserve(kappa_.size());
        for (const auto& kv : kappa_) {
            const MI& obs = kv.first.first;
            const MI& chg = kv.first.second;
            if (static_cast<int>(obs.size()) != d_ ||
                static_cast<int>(chg.size()) != s_)
                continue; // dimensions don't match this calculator; never queried
            kappaIdx_[encObs(obs) * powChg_ + encChg(chg)] = kv.second;
        }
    }

    MI decObs(long long code) const {
        MI n(d_, 0);
        for (int i = 0; i < d_; ++i) {
            n[i] = code % B_;
            code /= B_;
        }
        return n;
    }

    MI decChg(long long code) const {
        MI m(s_, 0);
        for (int i = 0; i < s_; ++i) {
            m[i] = code % B_;
            code /= B_;
        }
        return m;
    }

    std::string MIToString(const MI& x) const {
        std::ostringstream ss;
        ss << "{";
        for (size_t i = 0; i < x.size(); ++i) {
            if (i) ss << ",";
            ss << x[i];
        }
        ss << "}";
        return ss.str();
    }

    // Look up a flattened cumulant directly from pre-computed encodings, so the
    // hot callers never materialize unit / sum / zero charge multi-indices.
    double getKappaEnc(long long obsEnc, long long chgEnc) const
    {
        auto it = kappaIdx_.find(obsEnc * powChg_ + chgEnc);

        if (it == kappaIdx_.end()) {
            MI obs = decObs(obsEnc);
            MI chg = decChg(chgEnc);

            // Insert into the set to maintain a unique list
            missingKeys_.insert({ obs, chg });
            return 0.0;
        }

        return it->second;
    }

    double getKappa(const MI& n, const MI& m) const {
        return getKappaEnc(encObs(n), encChg(m));
    }

    double getLambda(int a, const MI& n) const {
        auto it = lambda_.find(static_cast<long long>(a) * powObs_ + encObs(n));
        if (it == lambda_.end()) return 0.0;
        return it->second;
    }

    void setLambda(int a, const MI& n, double value) {
        lambda_[static_cast<long long>(a) * powObs_ + encObs(n)] = value;
    }

    /// Colored set partitions of an n-element set into blocks painted with one
    /// of m charge types (Mathematica: ColoredPartitions[n, m]). Forwards to the
    /// thread-local cache in sam3detail (shared across calculator instances).
    const PoSmulti& getSetPoS(int n, int m) {
        return sam3detail::ColoredPartitionsOfSet(n, m);
    }

    /// Charge covariance matrix K_{ab} = kappa^gce_{0; e_a+e_b} (and its inverse),
    /// plus the charge means B0_[a] = kappa^gce_{0; e_a}. K must be non-singular.
    void buildKMatrix() {
        MI zeroD = zeroMI(d_);

        K_ = Eigen::MatrixXd::Zero(s_, s_);
        B0_.assign(s_, 0.0);

        for (int a = 0; a < s_; ++a) {
            B0_[a] = getKappa(zeroD, unitMI(s_, a));

            for (int b = 0; b < s_; ++b) {
                MI mab = addMI(unitMI(s_, a), unitMI(s_, b));
                K_(a, b) = getKappa(zeroD, mab);
            }
        }

        Klu_ = Eigen::FullPivLU<Eigen::MatrixXd>(K_);

        if (Klu_.rank() < s_) {
            throw std::runtime_error("Conserved-charge covariance matrix K is singular.");
        }

        Kinv_ = K_.inverse();
    }

    /// Charge-color counts q of a colored partition: q[j] = number of blocks
    /// painted with charge type j (Mathematica: the per-partition vector q).
    const MI& getqpi(const sam3detail::ColoredPartition& part) {
        qpiScratch_.assign(s_, 0);
        for(const auto& block : part)
            qpiScratch_[block.first]++;
        return qpiScratch_;
    }

    // Map each set element index (0 .. total(r)-1) to its observable dimension.
    // Depends only on r; built once per r into a reused scratch buffer.
    const std::vector<int>& buildElm(const MI& r) {
        elmScratch_.clear();
        for (std::size_t dim = 0; dim < r.size(); ++dim)
            for (int i = 0; i < r[dim]; ++i)
                elmScratch_.push_back(static_cast<int>(dim));
        return elmScratch_;
    }

    /// Induced observable multi-index of one block: counts how many of the
    /// block's set elements belong to each observable dimension, using the
    /// precomputed element->dimension map elm (Mathematica: InducedMI).
    const MI& getkb(const sam3detail::ColoredBlock& part,
                    const std::vector<int>& elm) {
        kbScratch_.assign(d_, 0);
        for (const auto& el : part.second)
            kbScratch_[elm[el]]++;
        return kbScratch_;
    }

    /// B_{r,q}: colored-partition building block -- sum, over partitions of the
    /// |r|-element set whose charge-color signature equals q, of the product of
    /// lambda coefficients over the (induced) blocks. Lazy, memoized in Brq_.
    double computeAddBrq(const MI& r, const MI& q) {
        const long long key = encObs(r) * powChg_ + encChg(q);
        auto it = Brq_.find(key);
        if (it != Brq_.end()) return it->second;

        const auto& elm = buildElm(r);
        const long long qEnc = encChg(q);
        double ret = 0.;
        const auto& PoS = getSetPoS(total(r), s_);
        for(const auto& part : PoS) {
            if (encChg(getqpi(part)) != qEnc)
                continue;
            double tret = 1.;
            for(const auto& block : part) {
                int j = block.first;
                tret *= getLambda(j, getkb(block, elm));
            }
            ret += tret;
        }

        return Brq_[key] = ret;
    }

    /// Compute B_{r,q} for *all* q in a single pass over the partitions of r
    /// (the per-q computeAddBrq() above would rescan them once per q). The
    /// single whole-set block (q of total order 1) is skipped here: that term
    /// needs lambda_{.,r} of the *same* order as r, which is not yet finalized
    /// when this runs during the lambda solve, so it is deferred to the lazy
    /// computeAddBrq() at its (later, safe) point of use. Bucketed into Brq_.
    void computeAddBrqAllqs(const MI& r) {
        const long long rEnc = encObs(r);
        if (Brqscomputed_[rEnc])
            return;
        const auto& elm = buildElm(r);
        std::unordered_map<long long, double> rets;
        const auto& PoS = getSetPoS(total(r), s_);
        for(const auto& part : PoS) {
            if (part.size() == 1)
                continue;
            const long long qEnc = encChg(getqpi(part));
            double tret = 1.;
            for(const auto& block : part) {
                int j = block.first;
                tret *= getLambda(j, getkb(block, elm));
            }
            rets[qEnc] += tret;
        }

        for(const auto &el : rets)
            Brq_[rEnc * powChg_ + el.first] = el.second;

        Brqscomputed_[rEnc] = 1;
    }

    /// Remainder R_{a;n}: the inhomogeneous (already-known lower-order) term of
    /// the lambda recursion at order n for charge a (Mathematica: RemainderR).
    /// Memoized in R_.
    double computeAddRan(int a, const MI& n) {
        const long long key = static_cast<long long>(a) * powObs_ + encObs(n);
        auto it = R_.find(key);
        if (it != R_.end()) return it->second;

        double ret = 0.;
        const auto& allps = subMI(n);
        for(const auto& p : allps) {
            const MI tr = subMIvec(n, p);
            int rTot = total(tr);

            computeAddBrqAllqs(tr);

            const long long pObs = encObs(p);
            const long long aOff = powA_[a];
            double tret2 = 0.;
            for (int Q = 1; Q <= rTot; ++Q) {
                for (const MI& q : weakMI(s_, Q)) {
                    if (isZeroMI(p) && total(q) == 1)
                        continue;
                    tret2 += getKappaEnc(pObs, encChg(q) + aOff) * computeAddBrq(tr, q);
                }
            }
            ret += miBinom(n, p) * tret2;
        }

        return R_[key] = ret;
    }

    /// Saddle-point coefficient lambda_{j;n} = -[ K^{-1} ( kappa^gce_{n; e_a}
    /// + R_{a;n} ) ]_j (Mathematica: KappaSaddlePoint). Memoized in lambda_.
    double computeAddLambdajn(int j, const MI& n) {
        const long long key = static_cast<long long>(j) * powObs_ + encObs(n);
        auto it = lambda_.find(key);
        if (it != lambda_.end()) return it->second;

        const long long nObs = encObs(n);
        double ret = 0.;
        for(int a = 0; a < s_; ++a) {
            ret += -Kinv_(j, a) * (getKappaEnc(nObs, powA_[a]) + computeAddRan(a, n));
        }

        return lambda_[key] = ret;
    }

    void computeAllLambdasOrderN(int nabs) {
        const auto& allns = weakMI(d_, nabs);
        for(int j = 0; j < s_; ++j) {
            for(const auto& n: allns) {
                computeAddLambdajn(j, n);
            }
        }
    }

    /// Solve all lambda_{a;n} up to total order nmax. Order Nmax-1 suffices:
    /// the order-N canonical cumulant only consumes lambdas up to order N-1.
    void computeAllLambdas(int nmax) {
        for(int nabs = 1; nabs <= nmax; ++nabs) {
            computeAllLambdasOrderN(nabs);
        }
    }

    /// Canonical cumulant kappa^ce_n: the colored multivariate Faa di Bruno sum
    /// over partitions of the n observable derivatives, each partition weighted
    /// by kappa^gce of the unpartitioned part and the lambda coefficients of the
    /// blocks (Mathematica: CanonicalFromSaddle / BellSumKernel with zero offset).
    double computeKappaCEn(const MI& n) {
        double ret = getKappaEnc(encObs(n), 0);
        const auto& allps = subMI(n);
        for(const auto& p : allps) {
            MI nmp = subMIvec(n, p);
            auto absnmp = total(nmp);
            if (absnmp == 0)
                continue;
            const long long pObs = encObs(p);
            const auto& elm = buildElm(nmp);
            const auto& PoS = getSetPoS(absnmp, s_);
            double tret = 0.;
            for(const auto& part : PoS) {
                if (part.size() > absnmp)
                    continue;
                if (isZeroMI(p) && part.size() == 1)
                    continue;
                const auto& qpi = getqpi(part);
                double tret2 = getKappaEnc(pObs, encChg(qpi));
                for(const auto& block : part) {
                    int j = block.first;
                    tret2 *= getLambda(j, getkb(block, elm));
                }
                tret += tret2;
            }
            ret += miBinom(n, p) * tret;
        }
        return ret;
    }

    void solveAllLambdaCoefficients() {
        computeAllLambdas(Nmax_ - 1);
    }
};

/// Compute the SAM-3.0 canonical cumulants of `observableDim` subsystem
/// observables under global conservation of `chargeDim` correlated charges,
/// from grand-canonical input (see the file header for the method and
/// references). This is the main public interface.
///
/// \param observableDim  number of observables d (>= 1)
/// \param chargeDim      number of conserved charges s (>= 1)
/// \param desiredOrder   highest total cumulant order Nmax
/// \param gceCumulants   grand-canonical cumulants kappa^gce, keyed by
///                       (observable MI of length d, charge MI of length s);
///                       missing entries are treated as zero (and reported
///                       on stderr after the computation)
/// \return  canonical cumulants kappa^ce_n keyed by the observable MI n, for
///          every n with total order 1 .. Nmax
inline CEMap ComputeSAM3CanonicalCumulants(int observableDim,
    int chargeDim,
    int desiredOrder,
    const KappaMap& gceCumulants) {
    SAM3Calculator calc(observableDim, chargeDim, desiredOrder, gceCumulants);
    CEMap result = calc.compute();

    // Print the unique missing cumulants after processing
    calc.printMissingCumulants();

    return result;
}

#endif // SAM3_H