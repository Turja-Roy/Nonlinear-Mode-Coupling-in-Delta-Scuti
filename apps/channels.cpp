/* Direct vs parametric census of the resonant triplets.

       ./channels --tag dsct_M2.0

   Reads observables_<tag>.csv, writes channels_<tag>.csv: one row per coupling
   channel with counts at both the radial-triplet and the (triplet, m) level,
   and the frequency band of the daughter each channel drives. */

#include "csv.hpp"
#include "stability.hpp"

#include <CLI/CLI.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>

int main(int argc, char** argv) {
    CLI::App app{"Direct vs parametric census of the resonant triplets"};
    std::string tag = "dsct_M2.0", out = "out";
    app.add_option("--tag", tag);
    app.add_option("--out", out);
    CLI11_PARSE(app, argc, argv);

    const std::filesystem::path src = std::filesystem::path(out) / ("observables_" + tag + ".csv");
    const csv::Table df = csv::read(src);
    const size_t n = df.rows();

    const auto &ga = df.numbers("gamma_a"), &gb = df.numbers("gamma_b"),
               &gc = df.numbers("gamma_c");
    const auto &wa = df.numbers("omega_a"), &wb = df.numbers("omega_b"),
               &wc = df.numbers("omega_c");
    const auto &mu_max = df.numbers("mu_max"), &delta = df.numbers("delta");

    /* Direct channels have a single daughter; the parametric channel feeds two,
       so quote their band instead. */
    std::vector<stab::Channel> cls(n);
    std::vector<double> f_lo(n), f_hi(n), frac(n);
    for (size_t i = 0; i < n; ++i) {
        cls[i] = stab::channel(ga[i], gb[i], gc[i]);
        const double w[3] = {std::abs(wa[i]) / CD, std::abs(wb[i]) / CD, std::abs(wc[i]) / CD};
        const int d = stab::daughter_slot(ga[i], gb[i], gc[i]);
        if (cls[i] == stab::Channel::parametric) {
            f_lo[i] = std::min(w[1], w[2]);
            f_hi[i] = std::max(w[1], w[2]);
        } else if (d >= 0) {
            f_lo[i] = f_hi[i] = w[d];
        } else {
            f_lo[i] = f_hi[i] = std::nan("");
        }
        frac[i] = std::abs(delta[i] / wa[i]);
    }

    // The radial triplet is the row's six (l, n) keys; gamma is m-degenerate,
    // so the channel is a property of that, not of the m combination.
    const char* KEYS[6] = {"l_a", "n_a", "l_b", "n_b", "l_c", "n_c"};
    std::set<std::array<int, 6>> seen;
    std::vector<bool> first(n, false);
    for (size_t i = 0; i < n; ++i) {
        std::array<int, 6> k{};
        for (int j = 0; j < 6; ++j) k[j] = int(df.numbers(KEYS[j])[i]);
        first[i] = seen.insert(k).second;
    }

    struct Out { const char* name; long n_radial, n_m_rows;
                 double mu_med, mu_max_, frac_med, q10, q50, q90; };
    std::vector<Out> table;
    for (stab::Channel ch : stab::CHANNELS) {
        std::vector<double> mu, fr, lo, hi, both;
        long n_rad = 0;
        for (size_t i = 0; i < n; ++i) {
            if (cls[i] != ch) continue;
            mu.push_back(mu_max[i]);
            lo.push_back(f_lo[i]);
            hi.push_back(f_hi[i]);
            both.push_back(f_lo[i]);
            both.push_back(f_hi[i]);
            if (first[i]) { ++n_rad; fr.push_back(frac[i]); }
        }
        if (mu.empty()) continue;
        const bool any_f = std::any_of(lo.begin(), lo.end(),
                                       [](double v) { return !std::isnan(v); });
        table.push_back(Out{stab::channel_name(ch), n_rad, long(mu.size()),
                            csv::median(mu), *std::max_element(mu.begin(), mu.end()),
                            csv::median(fr),
                            any_f ? csv::quantile(lo, 0.10) : std::nan(""),
                            any_f ? csv::quantile(both, 0.50) : std::nan(""),
                            any_f ? csv::quantile(hi, 0.90) : std::nan("")});
    }
    std::stable_sort(table.begin(), table.end(),
                     [](const Out& a, const Out& b) { return a.n_radial > b.n_radial; });

    std::printf("%s: %zu radial triplets, %zu with m\n\n", src.c_str(), seen.size(), n);
    std::printf("%-12s %8s %9s %12s %12s %14s %10s %10s %10s\n", "channel", "n_radial",
                "n_m_rows", "mu_max_med", "mu_max_max", "frac_det_med",
                "f_q10", "f_q50", "f_q90");
    long dir_r = 0, dir_m = 0, par_r = 0, par_m = 0;
    for (const Out& o : table) {
        std::printf("%-12s %8ld %9ld %12.4g %12.4g %14.4g %10.4g %10.4g %10.4g\n",
                    o.name, o.n_radial, o.n_m_rows, o.mu_med, o.mu_max_, o.frac_med,
                    o.q10, o.q50, o.q90);
        if (std::string(o.name).rfind("direct", 0) == 0) { dir_r += o.n_radial; dir_m += o.n_m_rows; }
        if (std::string(o.name) == "parametric")         { par_r += o.n_radial; par_m += o.n_m_rows; }
    }
    std::printf("\ndirect     %6ld radial %7ld with m\n", dir_r, dir_m);
    std::printf("parametric %6ld radial %7ld with m\n", par_r, par_m);

    const std::filesystem::path dst = std::filesystem::path(out) / ("channels_" + tag + ".csv");
    csv::Writer w(dst, {"channel", "n_radial", "n_m_rows", "mu_max_median", "mu_max_max",
                        "frac_detuning_median", "f_daughter_q10", "f_daughter_q50",
                        "f_daughter_q90"});
    for (const Out& o : table)
        w.row(o.name, o.n_radial, o.n_m_rows, o.mu_med, o.mu_max_, o.frac_med,
              o.q10, o.q50, o.q90);
    std::printf("-> %s\n", dst.c_str());
    return 0;
}
