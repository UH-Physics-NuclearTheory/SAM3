// Cross-check bridge: emit a deterministic random GCE cumulant set AND the C++
// CE cumulants it produces (via sam3.h), as a Mathematica data file, so the
// notebook's KappaCE can be evaluated on the *same* input and compared.
//
//   g++ -std=c++17 -O2 bridge_sam3.cpp -I. -I/path/to/eigen -o bridge_sam3
//   ./bridge_sam3 d s Nmax seed out.m
//
// The file defines: dim, chg, nmax, kGCE[...] downvalues, and
//   ceCpp = { {observableIndex, value}, ... }.

#include "sam3.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

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
static void writeMI(std::ostream& os, const MI& v) {
    os << "{";
    for (std::size_t i = 0; i < v.size(); ++i) { if (i) os << ", "; os << v[i]; }
    os << "}";
}

int main(int argc, char** argv) {
    int d    = argc > 1 ? std::atoi(argv[1]) : 1;
    int s    = argc > 2 ? std::atoi(argv[2]) : 2;
    int Nmax = argc > 3 ? std::atoi(argv[3]) : 3;
    unsigned seed = argc > 4 ? static_cast<unsigned>(std::atol(argv[4])) : 20260624u;
    std::string outfile = argc > 5 ? argv[5] : "bridge_sam3.m";

    std::mt19937 rng(seed);
    KappaMap kappa = makeRandomKappa(d, s, Nmax, Nmax + 1, rng);
    CEMap ce = ComputeSAM3CanonicalCumulants(d, s, Nmax, kappa);

    std::ofstream os(outfile);
    if (!os) { std::cerr << "cannot open " << outfile << "\n"; return 1; }
    os << std::setprecision(17);
    os << "(* SAM-3.0 bridge data (C++ sam3.h -> Mathematica). seed=" << seed << " *)\n";
    os << "dim = " << d << "; chg = " << s << "; nmax = " << Nmax << ";\n";
    os << "ClearAll[kGCE];\n";
    for (const auto& kv : kappa) {
        os << "kGCE[";
        writeMI(os, kv.first.first);
        os << ", ";
        writeMI(os, kv.first.second);
        os << "] = " << kv.second << ";\n";
    }
    os << "kGCE[_List, _List] := 0;\n";
    os << "ceCpp = {\n";
    bool first = true;
    for (int tot = 1; tot <= Nmax; ++tot)
        for (const MI& n : weakComps(d, tot)) {
            auto it = ce.find(n);
            if (it == ce.end()) continue;
            if (!first) os << ",\n";
            first = false;
            os << "  {"; writeMI(os, n); os << ", " << it->second << "}";
        }
    os << "\n};\n";
    os.close();

    std::cout << "wrote " << outfile << "  (d=" << d << " s=" << s
              << " Nmax=" << Nmax << ")\n";
    return 0;
}
