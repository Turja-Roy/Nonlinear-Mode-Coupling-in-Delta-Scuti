/* Three-mode coupling tables for one model: triplets -> kappa -> observables.

       ./three_mode_search --model models/dsct_M2.0

   Writes three tables under --out, tagged with the model directory name:

       triplets_<tag>.csv     one row per radial triplet
       kappa_<tag>.csv        one row per (triplet, m), with kappa
       observables_<tag>.csv  the same rows plus mu and E_th, ranked by mu_max

   kappa_all_m dominates the runtime and observables() calls it already, so
   both m-level tables are built from one pass rather than two.

   For the wide net the defaults do not scale: point --detail-dir at
   detail_wide, restrict the sum mode with --parents-only, and rank on --m000. */

#include "csv.hpp"
#include "stability.hpp"

#include <CLI/CLI.hpp>

#include <chrono>
#include <cstdio>

int main(int argc, char** argv) {
    CLI::App app{"Three-mode coupling tables for one model"};
    std::string model = "models/dsct_M2.0", tag, out = "out";
    std::string detail_dir = "detail", inlist = "gyre_ad.in", gamma_mode = "rad";
    std::vector<std::string> nad = {"summary_nad.h5"};
    int l_max = 3, jobs = 1;
    double cut_dimless = DETUNING_CUT_DIMLESS;
    bool parents_only = false, m000 = false, force = false;

    app.add_option("--model", model);
    app.add_option("--tag", tag, "defaults to the model dir name");
    app.add_option("--out", out);
    app.add_option("--l-max", l_max);
    app.add_option("--cut", cut_dimless, "detuning cut in units of sqrt(GM/R^3)");
    app.add_option("--detail-dir", detail_dir,
                   "eigenfunction dumps under <model>/gyre; detail_wide for the l <= 25 net");
    app.add_option("--inlist", inlist);
    app.add_option("--nad", nad)->expected(-1);
    app.add_option("--gamma", gamma_mode, "tot adds the turbulent rate to the radiative one, "
                   "which also moves the parent/daughter split")
        ->check(CLI::IsMember({"rad", "tot"}));
    app.add_flag("--parents-only", parents_only,
                 "restrict the sum mode to the self-excited modes; the pair still "
                 "ranges over the whole net");
    app.add_flag("--m000", m000, "kappa at m = (0, 0, 0) only");
    app.add_option("-j,--jobs", jobs, "threads for the kappa quadrature; loading and "
                   "enumeration stay serial");
    app.add_flag("--force", force, "overwrite existing tables for this tag");
    CLI11_PARSE(app, argc, argv);

    const std::filesystem::path root(model), outdir(out);
    if (tag.empty()) tag = root.filename().string();
    const std::filesystem::path p_trip = outdir / ("triplets_" + tag + ".csv");
    const std::filesystem::path p_kap  = outdir / ("kappa_" + tag + ".csv");
    const std::filesystem::path p_obs  = outdir / ("observables_" + tag + ".csv");
    if (!force)
        for (const auto& p : {p_trip, p_kap, p_obs})
            if (std::filesystem::exists(p)) {
                std::fprintf(stderr, "refusing to overwrite %s  (pass --force)\n",
                             p.c_str());
                return 1;
            }

    const auto t0 = std::chrono::steady_clock::now();
    auto secs = [&] {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };

    std::printf("loading %s ...\n", model.c_str());
    Model m = load_model(root, detail_dir, inlist, nad);
    const Star& star = *m.star;
    const double cut = cut_dimless * star.omega_dyn();
    std::printf("  %zu modes, sqrt(GM/R^3) = %.6e rad/s, detuning cut = %.6e rad/s\n",
                m.eigfuncs.size(), star.omega_dyn(), cut);

    if (gamma_mode == "tot") {
        int before = 0, after = 0;
        for (auto& [k, ef] : m.eigfuncs) {
            before += ef.gamma < 0.0;
            ef.gamma += gamma_turb(ef);          // one double, not five arrays
            after += ef.gamma < 0.0;
        }
        std::printf("  gamma_rad + gamma_turb: driven modes %d -> %d\n", before, after);
    }

    std::set<Key> sum_keys;
    if (parents_only) {
        for (const auto& [k, ef] : m.eigfuncs) if (ef.gamma < 0.0) sum_keys.insert(k);
        std::printf("  sum mode restricted to %zu self-excited modes\n", sum_keys.size());
    }

    std::printf("enumerating triplets ...\n");
    auto trips = enumerate_triplets(m.eigfuncs, cut, l_max,
                                    parents_only ? &sum_keys : nullptr);
    std::printf("  %zu radial triplets  [%.1fs]\n", trips.size(), secs());
    if (!m000) std::printf("  %ld with m\n", count_with_m(trips));

    {
        std::vector<std::string> hdr = {"l_a", "n_a", "l_b", "n_b", "l_c", "n_c",
                                        "omega_a", "omega_b", "omega_c", "delta",
                                        "abs_delta_over_omega_a"};
        if (!m000) hdr.push_back("n_m");
        csv::Writer w(p_trip, hdr);
        for (const auto& r : triplet_rows(trips, !m000)) {
            if (m000) w.row(r.l_a, r.n_a, r.l_b, r.n_b, r.l_c, r.n_c,
                            r.omega_a, r.omega_b, r.omega_c, r.delta, r.abs_delta_over_omega_a);
            else      w.row(r.l_a, r.n_a, r.l_b, r.n_b, r.l_c, r.n_c,
                            r.omega_a, r.omega_b, r.omega_c, r.delta, r.abs_delta_over_omega_a,
                            r.n_m);
        }
    }
    std::printf("  -> %s\n", p_trip.c_str());

    std::printf("kappa, mu, thresholds ...\n");
    const auto rows = stab::build_rows(trips, m.eigfuncs, star, m000, jobs);
    std::printf("  %zu rows on %d thread(s)  [%.1fs]\n", rows.size(), jobs, secs());

    stab::write_rows(p_obs, rows);
    std::printf("  -> %s\n", p_obs.c_str());

    {
        csv::Writer w(p_kap, {"l_a", "n_a", "m_a", "l_b", "n_b", "m_b", "l_c", "n_c", "m_c",
                              "omega_a", "omega_b", "omega_c", "delta", "kappa", "refine",
                              "max_abs_n", "abs_kappa", "mu"});
        for (const stab::Row& r : rows)
            w.row(r.l_a, r.n_a, r.m_a, r.l_b, r.n_b, r.m_b, r.l_c, r.n_c, r.m_c,
                  r.omega_a, r.omega_b, r.omega_c, r.delta, r.kappa, r.refine,
                  std::max({std::abs(r.n_a), std::abs(r.n_b), std::abs(r.n_c)}),
                  std::abs(r.kappa), r.mu_c);
    }
    std::printf("  -> %s\ndone in %.1fs\n", p_kap.c_str(), secs());
    return 0;
}
