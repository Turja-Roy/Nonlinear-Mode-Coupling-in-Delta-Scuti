/* Integrate a mixed direct/parametric network to a bounded state.

   Direct and parametric coupling act on the same modes (MW25 sec 5). Which mode
   carries the parametric leg decides whether the network has a bounded state:

       daughter   a + b -> c,  c -> d + d
       parent     a + b -> c,  a -> d + d,  b -> d + d
       pair       a + b -> c,  a -> d1 + d2,  b -> d1 + d2

   Parents are self-excited and nothing else in `daughter` can absorb that flux,
   so its parents run away and are only integrable as a fixed-amplitude pump --
   the physical picture for modes whose saturation lies outside the network. In
   `parent` and `pair` the parametric leg drains the parents themselves, which
   is MW25's stabilisation mechanism; those run unpumped.

   Topology is data, not equations: each entry of TOPOLOGIES returns a mode list
   and a list of Triplets, and Network integrates any number of them.

   Selection is not on growth rate alone. Gamma/gamma_d ranks the strength of
   the instability, but the cost of the integration is set by Delta: the
   envelope oscillates at Delta while the amplitudes evolve on 1/gamma, so a
   system with Gamma << |Delta| needs ~|Delta|/Gamma oscillations per
   e-folding. Only near-resonant systems are integrable. */

#include "amplitude.hpp"
#include "csv.hpp"
#include "stability.hpp"

#include <CLI/CLI.hpp>

#include <cstdio>
#include <functional>
#include <numeric>

using amp::Network;
using amp::State;
using amp::Triplet;

namespace {

constexpr double YR = 3.15576e7;

double cell(const csv::Table& df, const std::string& c, size_t i) {
    return df.numbers(c)[i];
}

// Modes for the (l, n) columns named by `labels`, frequencies in c/d. omega
// is signed by this mode's own gamma (negative = self-excited), which the
// amplitude equations need and the CSV's plain-magnitude f_* columns don't
// carry. Right here because every leg is direct-sum or parametric, where the
// gamma-odd mode is the sum mode -- see amplitude.hpp.
std::vector<amp::Mode> modes_of(const csv::Table& df, size_t i, const std::string& labels) {
    std::vector<amp::Mode> out;
    for (char s : labels) {
        const std::string t(1, s);
        char buf[64];
        std::snprintf(buf, sizeof buf, "%c(%d,%+d)", s,
                      int(cell(df, "l_" + t, i)), int(cell(df, "n_" + t, i)));
        const double gamma = cell(df, "gamma_" + t, i);
        const double omega = (gamma < 0.0 ? -1.0 : 1.0) * std::abs(cell(df, "f_" + t, i) * CD);
        out.push_back(amp::Mode{buf, omega, gamma});
    }
    return out;
}

// a + b -> c. The triplet's 3 slots are symmetric now (no sum-slot role), so
// which one holds which mode no longer matters to the equations.
Triplet direct_leg(const csv::Table& df, size_t i) {
    return Triplet{2, 0, 1, cell(df, "kappa_direct", i)};
}

struct Built { std::vector<amp::Mode> modes; std::vector<Triplet> cps;
               std::vector<int> pumped; };

Built daughter(const csv::Table& df, size_t i) {
    return {modes_of(df, i, "abcd"),
            {direct_leg(df, i), Triplet{2, 3, 3, cell(df, "kappa_param", i)}},
            {0, 1}};
}
Built parent(const csv::Table& df, size_t i) {
    return {modes_of(df, i, "abcd"),
            {direct_leg(df, i), Triplet{0, 3, 3, cell(df, "kappa_param_a", i)},
                                Triplet{1, 3, 3, cell(df, "kappa_param_b", i)}},
            {}};
}
Built pair_(const csv::Table& df, size_t i) {
    return {modes_of(df, i, "abcde"),
            {direct_leg(df, i), Triplet{0, 3, 4, cell(df, "kappa_param_a", i)},
                                Triplet{1, 3, 4, cell(df, "kappa_param_b", i)}},
            {}};
}

// name -> (builder, modes frozen by default under --pump)
const std::map<std::string, std::function<Built(const csv::Table&, size_t)>> TOPOLOGIES = {
    {"daughter", daughter}, {"parent", parent}, {"pair", pair_}};

/* Parent amplitude to start from, in units where E = |q|^2 E_star.

   daughter: E_c scales as q^4, so the amplitude that puts the direct daughter
   at its parametric threshold follows the quarter power. Driving above it puts
   Gamma above |Delta_p|/2, which is what makes the integration affordable.

   parent/pair: the parent crosses its own threshold under kappa-mechanism
   growth, so start below it and let the linear driving carry it across. */
double drive_amplitude(const csv::Table& df, size_t i, const std::string& topology) {
    if (topology == "daughter")
        return 1e-6 * std::pow(cell(df, "E_c_over_E_th", i), -0.25);
    return std::sqrt(cell(df, "E_th_par", i));
}

/* Slaving is a sub-threshold statement: above E_th the parametric leg drains c
   and |q_c| sits below the slaved value by construction. Verify it where it
   should hold, at half the threshold amplitude, where d decays and c relaxes
   cleanly. The response is s_c mu q_a q_b with s_c = 2 for distinct parents --
   the bare mu misses the combinatorial factor of the amplitude equation, which
   is what this measurement caught. */
void check_slaving(const Network& net, double q_thr) {
    const Triplet& d = net.triplets()[0];
    const double d1 = net.detuning(d);
    const double gc = net.modes()[2].gamma;
    const double slaved = 2.0 * stab::mu(d.kappa, net.modes()[2].omega, d1, gc)
                        * (0.5 * q_thr) * (0.5 * q_thr);
    State q0(net.modes().size(), 1e-4 * q_thr);
    for (size_t i = 0; i < net.modes().size(); ++i)
        if (net.modes()[i].gamma < 0.0) q0[i] = 0.5 * q_thr;
    const double t_end = std::max(8.0 / std::hypot(d1, gc), 6.0 / gc);
    amp::Options opt;
    opt.n_out = 4000;
    opt.rtol = 1e-10;
    opt.frozen = {0, 1};
    const amp::Solution sol = net.integrate(q0, t_end, opt);
    const auto& A = sol.y;

    std::vector<double> err;
    for (size_t j = A.size() / 2; j < A.size(); ++j)
        err.push_back(std::abs(std::abs(A[j][2]) - slaved) / slaved);
    std::printf("\n  q_c vs 2 mu q_a q_b at 0.5x threshold (sub-threshold, slaved): "
                "median %.3e\n", csv::median(err));
    std::printf("  E_d there: start %.3e -> end %.3e (decays, as it must below "
                "threshold)\n", std::norm(A.front()[3]), std::norm(A.back()[3]));
}

}  // namespace

int main(int argc, char** argv) {
    CLI::App app{"Integrate a mixed direct/parametric network to a bounded state"};
    std::string src = "out/four_mode_dsct_M2.0.csv", out = "out/mixed_network.csv";
    std::string topology = "auto";
    double q_parent = 0.0, drive = 0.0, e_folds = 0.0;
    int pump = -1;                                   // -1 = topology default
    app.add_option("--csv", src);
    app.add_option("--out", out);
    app.add_option("--topology", topology)
        ->check(CLI::IsMember({"auto", "daughter", "parent", "pair"}));
    app.add_option("--q-parent", q_parent, "parent amplitude; default set by the topology");
    app.add_option("--drive", drive, "multiple of the default parent amplitude "
                   "(3 for daughter, 0.5 for the self-saturating topologies)");
    app.add_flag("--pump{1},!--no-pump{0}", pump, "freeze the parents at their initial "
                 "amplitude; on by default for daughter, whose parents cannot saturate");
    app.add_option("--e-folds", e_folds, "integration length in e-folds of the slowest rate");
    CLI11_PARSE(app, argc, argv);

    const csv::Table df = csv::read(src);
    const bool has_branch = df.has("branch");
    if (topology == "auto") topology = has_branch ? df.text("branch")[0] : "daughter";

    std::vector<size_t> keep;
    for (size_t i = 0; i < df.rows(); ++i)
        if (!has_branch || df.text("branch")[i] == topology) keep.push_back(i);
    if (keep.empty()) {
        std::printf("no %s candidates in %s\n", topology.c_str(), src.c_str());
        return 1;
    }

    // Rank on the detuned threshold, not on Gamma/gamma_d: the latter drops the
    // [1 + Delta^2/(2 gamma_d)^2] factor, which is ~89^2 at the median.
    const std::string rank = topology == "daughter" ? "E_c_over_E_th" : "E_par_over_E_th";
    std::stable_sort(keep.begin(), keep.end(),
                     [&](size_t a, size_t b) { return cell(df, rank, a) > cell(df, rank, b); });
    long above = 0;
    for (size_t i : keep) above += cell(df, rank, i) > 1.0;
    std::printf("%zu %s candidates, %ld above threshold at q = 1e-06\n",
                keep.size(), topology.c_str(), above);
    const size_t r = keep.front();

    Built b = TOPOLOGIES.at(topology)(df, r);
    const Network net(b.modes, b.cps);
    std::printf("\nsystem: ");
    for (const auto& m : b.modes) std::printf("%s  ", m.name.c_str());
    std::printf("\n");
    for (const auto& m : b.modes)
        std::printf("  %-14s f %8.4f c/d   gamma %+.3e s^-1\n",
                    m.name.c_str(), m.omega / CD, m.gamma);
    for (const Triplet& tr : b.cps) {
        std::printf("  leg [");
        for (int i : {tr.a, tr.b, tr.c}) std::printf("%s  ", b.modes[i].name.c_str());
        std::printf("]  kappa %+.3f  Delta %+.3e c/d\n", tr.kappa, net.detuning(tr) / CD);
    }
    std::printf("  %s %.3e at q = 1e-6\n", rank.c_str(), cell(df, rank, r));

    const bool do_pump = pump < 0 ? topology == "daughter" : pump != 0;
    const std::vector<int> frozen = do_pump ? b.pumped : std::vector<int>{};
    if (drive == 0.0) drive = topology == "daughter" ? 3.0 : 0.5;
    const double q0_par = q_parent != 0.0 ? q_parent : drive * drive_amplitude(df, r, topology);
    std::printf("\nparent amplitude %.3e (%gx the default for %s)\n",
                q0_par, drive, topology.c_str());

    /* Seeds must sit far enough above atol that a sub-threshold decay does not
       underflow them: a daughter that reaches zero cannot be re-excited when the
       parent later crosses, and the parents then overshoot by orders of
       magnitude before roundoff noise regrows it. */
    State q0(b.modes.size(), 1e-3 * q0_par);
    for (size_t i = 0; i < b.modes.size(); ++i)
        if (b.modes[i].gamma < 0.0) q0[i] = q0_par;

    std::vector<double> rates;
    if (topology == "daughter")
        rates = {std::hypot(net.detuning(b.cps[0]), b.modes[2].gamma), b.modes[3].gamma};
    else
        for (const auto& m : b.modes) rates.push_back(std::abs(m.gamma));
    const double n_folds = e_folds != 0.0 ? e_folds : (topology == "daughter" ? 8.0 : 40.0);
    const double t_end = n_folds / *std::min_element(rates.begin(), rates.end());
    std::printf("%ssystem, t_end %.2f yr", do_pump ? "pumped " : "", t_end / YR);
    if (!frozen.empty()) {
        std::printf(", frozen: ");
        for (int i : frozen) std::printf("%s ", b.modes[i].name.c_str());
    }
    std::printf("\n");

    // A network with no way to absorb the kappa-mechanism flux has no bounded
    // state; the cap turns that into a finished run instead of a step size
    // shrinking without limit.
    amp::Options opt;
    opt.n_out = 4000;
    opt.rtol = 1e-10;
    opt.e_max = (30.0 * q0_par) * (30.0 * q0_par);
    opt.frozen = frozen;
    const amp::Solution sol = net.integrate(q0, t_end, opt);
    const auto& t = sol.t;
    const auto& A = sol.y;
    if (t.back() < t_end)
        std::printf("\n  RUNAWAY: amplitudes passed 30x the drive at t = %.3f yr -- "
                    "integration stopped\n", t.back() / YR);

    std::printf("\nintegrated %.3e yr\n", t.back() / YR);
    std::vector<std::vector<double>> E(b.modes.size(), std::vector<double>(t.size()));
    for (size_t j = 0; j < t.size(); ++j)
        for (size_t i = 0; i < b.modes.size(); ++i) E[i][j] = std::norm(A[j][i]);
    for (size_t i = 0; i < b.modes.size(); ++i)
        std::printf("  E_%-14s start %.3e  end %.3e  max %.3e\n", b.modes[i].name.c_str(),
                    E[i].front(), E[i].back(),
                    *std::max_element(E[i].begin(), E[i].end()));

    const size_t half = t.size() / 2;
    for (size_t i = 0; i < b.modes.size(); ++i) {
        if (b.modes[i].gamma <= 0.0) continue;
        const auto lo = E[i].begin() + long(half);
        const double mean = std::accumulate(lo, E[i].end(), 0.0) / double(t.size() - half);
        if (mean <= 0.0) continue;
        const auto mm = std::minmax_element(lo, E[i].end());
        std::printf("  %s energy swing over the second half: %.3f\n",
                    b.modes[i].name.c_str(), (*mm.second - *mm.first) / mean);
    }

    if (topology == "daughter" && do_pump) check_slaving(net, q0_par / drive);

    std::vector<std::string> hdr = {"t_yr"};
    for (const auto& m : b.modes) hdr.push_back("E_" + m.name);
    csv::Writer w(out, hdr);
    for (size_t j = 0; j < t.size(); ++j) {
        std::string line = csv::fmt(t[j] / YR);
        for (size_t i = 0; i < b.modes.size(); ++i) line += "," + csv::fmt(E[i][j]);
        w.row(line);
    }
    std::printf("\n-> %s\n", out.c_str());
    return 0;
}
