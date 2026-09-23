/* Generate a model's six GYRE inlists from its MESA .GYRE profile.

   Every window follows from the model:

     nu_dyn = sqrt(GM/R^3)                 -- sets the detuning cut
     Dnu    = 1 / (2 int dr/c_s)           -- p-mode large separation
     nu_ac  = c_s / (4 pi H_p) at the top  -- acoustic cutoff, caps the p scan
     I      = int (N/r) dr                 -- g-mode asymptotic integral
     K      = I * 86400 / (2 pi^2)         -- omega[c/d] = K sqrt(l(l+1)) / |n|

   K is the constant the wide net's per-l n_pg windows are built from.

       ./make_inlists models/dsct_M2.0_T7696 --profile profile255.data.GYRE

   The parametric-daughter band defaults to the dsct_M2.0 band scaled by
   nu_dyn. That is a scaling assumption, not a derivation: the true band is
   omega_c/2 over the direct daughter frequencies, known only after
   three_mode_search has enumerated triplets. Re-run with --daughter-band once
   they are, and diff the n_pg windows. */

#include "model.hpp"
#include "numeric.hpp"

#include <CLI/CLI.hpp>

#include <cstdarg>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>

namespace {

constexpr double DAY = 86400.0;
constexpr double MSUN = 1.988409870698051e33;
constexpr double RSUN = 6.957e10;

// dsct_M2.0 reference, for the scaled defaults
constexpr double REF_NU_DYN = 2.866;         // c/d
constexpr double REF_BAND_LO = 1.7, REF_BAND_HI = 5.2;   // c/d, parametric daughters
constexpr double REF_SCAN_LO = 0.8;          // c/d, bottom of the 'core' INVERSE scan
constexpr double REF_SPLIT = 25.0;           // c/d, INVERSE -> LINEAR handover
constexpr double REF_SCAN_HI = 95.0;         // c/d, top of the 'core' LINEAR scan
constexpr int CORE_N_PG = 20;                // |n_pg| for the parent / direct-daughter families

std::string sfmt(const char* f, ...) {
    va_list ap;
    va_start(ap, f);
    char buf[512];
    std::vsnprintf(buf, sizeof buf, f, ap);
    va_end(ap);
    return buf;
}

struct GyreModel { double M, R, nu_dyn, Dnu, nu_ac, I, K; int n; };

// M, R and the three integrals, straight off the GYRE input file.
GyreModel read_gyre_model(const std::filesystem::path& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("no such profile: " + path.string());
    std::string line;
    std::getline(f, line);
    GyreModel m{};
    { std::istringstream h(line); h >> m.n >> m.M >> m.R; }

    std::vector<double> r, P, rho, N2, G1;
    while (std::getline(f, line)) {
        std::istringstream s(line);
        std::vector<double> a;
        double v;
        while (s >> v) a.push_back(v);
        if (a.size() < 10 || !(a[1] > 0.0)) continue;    // drop r = 0, as the mask does
        r.push_back(a[1]); P.push_back(a[4]); rho.push_back(a[6]);
        N2.push_back(a[8]); G1.push_back(a[9]);
    }
    const eig::Index n = eig::Index(r.size());
    auto arr = [n](std::vector<double>& v) { return eig::Map<eig::ArrayXd>(v.data(), n); };
    const eig::ArrayXd R_ = arr(r), P_ = arr(P), rho_ = arr(rho), N2_ = arr(N2), G1_ = arr(G1);
    const eig::ArrayXd cs = (G1_ * P_ / rho_).sqrt();

    m.nu_dyn = std::sqrt(G * m.M / (m.R * m.R * m.R)) * DAY / (2.0 * M_PI);
    m.Dnu = 1.0 / (2.0 * num::trapz(1.0 / cs, R_)) * DAY;
    const double Hp = P_[n-1] / (rho_[n-1] * G * m.M / (m.R * m.R));
    m.nu_ac = cs[n-1] / (4.0 * M_PI * Hp) * DAY;
    m.I = num::trapz(N2_.max(0.0).sqrt() / R_, R_);
    m.K = m.I * DAY / (2.0 * M_PI * M_PI);
    return m;
}

// |n| in [K*Lam/hi, K*Lam/lo] per l, as (n_pg_min, n_pg_max), both < 0.
std::map<int, std::pair<int, int>> gnet_windows(double K, double lo, double hi,
                                                int l_lo, int l_hi, int floor_l, int floor_n) {
    std::map<int, std::pair<int, int>> out;
    for (int l = l_lo; l <= l_hi; ++l) {
        const double Lam = std::sqrt(double(l) * (l + 1));
        int n_min = int(std::ceil(K * Lam / hi));        // shallowest, highest frequency
        int n_max = int(std::floor(K * Lam / lo));       // deepest, lowest frequency
        if (l <= floor_l) n_min = std::max(n_min, floor_n);   // 'core' covers shallower
        if (n_max < n_min) continue;
        out[l] = {-n_max, -n_min};
    }
    return out;
}

std::string mode_block(int l, int n_lo, int n_hi, const std::string& tag) {
    return sfmt("\n&mode\n  l = %d\n  n_pg_min = %d\n  n_pg_max = %d\n  tag = '%s'\n/\n",
                l, n_lo, n_hi, tag.c_str());
}

void write(const std::filesystem::path& p, const std::string& text, bool force) {
    if (std::filesystem::exists(p) && !force)
        throw std::runtime_error("refusing to overwrite " + p.string() + " (pass --force)");
    std::ofstream(p) << text;
    std::printf("  wrote %s\n", p.c_str());
}

}  // namespace

int main(int argc, char** argv) {
    CLI::App app{"Generate a model's six GYRE inlists from its MESA .GYRE profile"};
    std::string model_dir, profile;
    std::vector<double> band_arg;
    int l_max_core = 6, l_max_pass1 = 15, l_max_pass2 = 25;
    bool force = false;
    app.add_option("model_dir", model_dir)->required();
    app.add_option("--profile", profile, "basename of the .GYRE file under mesa/LOGS/")
        ->required();
    app.add_option("--daughter-band", band_arg,
                   "parametric daughter band in c/d; default scales the dsct_M2.0 "
                   "band by nu_dyn")->expected(2);
    app.add_option("--l-max-core", l_max_core);
    app.add_option("--l-max-pass1", l_max_pass1);
    app.add_option("--l-max-pass2", l_max_pass2);
    app.add_flag("--force", force, "overwrite existing inlists");
    CLI11_PARSE(app, argc, argv);

    const std::filesystem::path root(model_dir), gyre = root / "gyre";
    std::filesystem::create_directories(gyre);
    // GYRE names these in detail_template but will not create them.
    std::filesystem::create_directories(gyre / "detail");
    std::filesystem::create_directories(gyre / "detail_wide");

    const GyreModel m = read_gyre_model(
        std::filesystem::absolute(root / "mesa" / "LOGS" / profile));
    const double scale = m.nu_dyn / REF_NU_DYN;
    const double band_lo = band_arg.empty() ? REF_BAND_LO * scale : band_arg[0];
    const double band_hi = band_arg.empty() ? REF_BAND_HI * scale : band_arg[1];

    /* The 'core' windows scale with nu_dyn off the hand-derived dsct_M2.0 ones,
       capped at the acoustic cutoff -- above it there is nothing to find. Dnu
       and nu_ac are reported in the header either way, so a bad window is
       visible without running GYRE. */
    const double core_lo = REF_SCAN_LO * scale;
    const double split = REF_SPLIT * scale;
    const double p_hi = std::min(REF_SCAN_HI * scale, 1.02 * m.nu_ac);
    if (p_hi <= split) {
        std::fprintf(stderr, "acoustic cutoff %.1f c/d is below the p-mode scan split "
                     "%.1f c/d -- check the model\n", m.nu_ac, split);
        return 1;
    }
    // 'gnet' scan is deliberately wider than the band the n_pg windows encode,
    // so the n_pg limits bind and an error in K cannot silently drop modes.
    const double g_lo = 0.88 * band_lo, g_hi = 1.16 * band_hi;

    const std::string hdr = sfmt(
        "! GENERATED by make_inlists -- edit that, not this.\n!\n"
        "! Model: %s   M = %.3f Msun   R = %.4f Rsun\n"
        "!   sqrt(GM/R^3) = %.4f c/d        detuning cut 0.15*nu_dyn = %.4f c/d\n"
        "!   large separation Dnu = %.3f c/d   acoustic cutoff = %.1f c/d\n"
        "!   int(N/r)dr = %.4e s^-1              K = %.3f c/d  (omega = K*Lam/|n|)\n",
        root.filename().c_str(), m.M / MSUN, m.R / RSUN, m.nu_dyn,
        0.15 * m.nu_dyn, m.Dnu, m.nu_ac, m.I, m.K);

    std::printf("\n%s:  nu_dyn=%.4f  Dnu=%.3f  nu_ac=%.1f  K=%.3f c/d\n",
                root.filename().c_str(), m.nu_dyn, m.Dnu, m.nu_ac, m.K);
    std::printf("  core scan %.2f-%.1f-%.1f c/d;  daughter band %.2f-%.2f c/d%s\n",
                core_lo, split, p_hi, band_lo, band_hi,
                band_arg.empty() ? "  (scaled)" : "  (given)");

    const std::string model_blk = sfmt(
        "\n&model\n  model_type = 'EVOL'\n  file = '../mesa/LOGS/%s'\n"
        "  file_format = 'MESA'\n/\n", profile.c_str());
    auto osc = [](const std::string& nad) {
        return "\n&osc\n  outer_bound = 'UNNO'\n  variables_set = 'DZIEM'\n" + nad + "/\n";
    };
    auto core_scans = [&](const std::string& tags) {
        return sfmt("\n&scan\n  grid_type = 'INVERSE'\n  freq_min = %.2f\n"
                    "  freq_max = %.1f\n  n_freq = 1000\n  freq_units = 'CYC_PER_DAY'\n%s/\n"
                    "\n&scan\n  grid_type = 'LINEAR'\n  freq_min = %.1f\n"
                    "  freq_max = %.1f\n  n_freq = 1000\n  freq_units = 'CYC_PER_DAY'\n%s/\n",
                    core_lo, split, tags.c_str(), split, p_hi, tags.c_str());
    };
    auto grid_core = [](const std::string& tags) {
        return "\n&grid\n  w_osc = 10\n  w_exp = 2\n  w_ctr = 10\n" + tags + "/\n";
    };

    const std::string det_items =
        "id,l,m,n_pg,n_p,n_g,omega,freq,x,xi_r,xi_h,eul_Phi,deul_Phi,lag_rho,rho,P,"
        "Gamma_1,As,c_1,V_2,M_r,M_star,R_star,L_star,E_norm";
    const std::string sum_items =
        "id,l,m,n_pg,n_p,n_g,omega,freq,E_norm,Delta_p,Delta_g,M_star,R_star,L_star";
    const std::string nad_items =
        "id,l,m,n_pg,n_p,n_g,omega,freq,E_norm,M_star,R_star,L_star";

    /* ------------------------------------------------------- narrow run */
    std::string narrow_modes;
    for (int l = 0; l < 4; ++l)
        narrow_modes += mode_block(l, -CORE_N_PG, CORE_N_PG, "l" + std::to_string(l));

    write(gyre / "gyre_ad.in",
          hdr + "! Narrow run: l <= 3, |n_pg| <= 20. Eigenfunctions for the coupling tables.\n"
          + "&constants\n/\n" + model_blk + narrow_modes + osc("") + "\n&rot\n/\n"
          + "\n&num\n  diff_scheme = 'COLLOC_GL4'\n/\n" + core_scans("") + grid_core("")
          + sfmt("\n&ad_output\n  summary_file = 'summary_ad.h5'\n"
                 "  detail_template = 'detail/detail.l%%l.n%%n.h5'\n"
                 "  summary_item_list = '%s'\n  detail_item_list = '%s'\n"
                 "  freq_units = 'CYC_PER_DAY'\n/\n\n&nad_output\n/\n",
                 sum_items.c_str(), det_items.c_str()), force);

    write(gyre / "gyre_nad.in",
          hdr + "! Nonadiabatic run for damping rates only: gamma = -Im(omega).\n"
                "! Eigenfunctions come from the adiabatic run, so no detail files here.\n"
          + "&constants\n/\n" + model_blk + narrow_modes
          + osc("  nonadiabatic = .TRUE.\n") + "\n&rot\n/\n"
          + "\n&num\n  diff_scheme = 'MAGNUS_GL2'\n  nad_search = 'AD'\n/\n"
          + core_scans("") + grid_core("") + "\n&ad_output\n/\n"
          + sfmt("\n&nad_output\n  summary_file = 'summary_nad.h5'\n"
                 "  summary_item_list = '%s'\n  freq_units = 'CYC_PER_DAY'\n/\n",
                 nad_items.c_str()), force);

    /* ---------------------------------------------------- wide net passes */
    const std::string gnet_scan = sfmt(
        "\n! Wider than the %.2f-%.2f c/d band the n_pg limits encode, so the\n"
        "! n_pg windows bind and a few percent of error in K cannot silently drop modes.\n"
        "&scan\n  grid_type = 'INVERSE'\n  freq_min = %.2f\n  freq_max = %.2f\n"
        "  n_freq = 3000\n  freq_units = 'CYC_PER_DAY'\n  tag_list = 'gnet'\n/\n",
        band_lo, band_hi, g_lo, g_hi);
    const std::string grid_gnet =
        "\n&grid\n  w_osc = 15\n  w_exp = 2\n  w_ctr = 10\n  tag_list = 'gnet'\n/\n";

    const auto w1 = gnet_windows(m.K, band_lo, band_hi, 2, l_max_pass1, l_max_core, 21);
    const auto w2 = gnet_windows(m.K, band_lo, band_hi, l_max_pass1 + 1, l_max_pass2,
                                 l_max_core, 21);
    auto count = [](const std::map<int, std::pair<int, int>>& w) {
        long s = 0;
        for (const auto& [l, r] : w) s += r.second - r.first + 1;
        return s;
    };
    std::printf("  gnet pass 1 (l 2-%d): %ld modes;  pass 2 (l %d-%d): %ld modes\n",
                l_max_pass1, count(w1), l_max_pass1 + 1, l_max_pass2, count(w2));

    struct Pass { const char* name; const std::map<int, std::pair<int, int>>* w;
                  bool core; const char* summ; std::string note; };
    const Pass passes[2] = {
        {"wide", &w1, true, "summary_ad_wide.h5",
         sfmt("! Wide daughter net, pass 1 -- l <= %d.\n", l_max_pass1)},
        {"wide_hi", &w2, false, "summary_ad_wide_hi.h5",
         sfmt("! Wide daughter net, pass 2 -- gnet only, l = %d-%d.\n",
              l_max_pass1 + 1, l_max_pass2)},
    };
    for (const Pass& p : passes) {
        std::string modes;
        if (p.core)
            for (int l = 0; l <= l_max_core; ++l)
                modes += mode_block(l, -CORE_N_PG, CORE_N_PG, "core");
        for (const auto& [l, r] : *p.w) modes += mode_block(l, r.first, r.second, "gnet");

        const std::string scans =
            (p.core ? core_scans("  tag_list = 'core'\n") : std::string()) + gnet_scan;
        const std::string grids =
            (p.core ? grid_core("  tag_list = 'core'\n") : std::string()) + grid_gnet;
        // eul_P/lag_P/eul_rho/U are already absent: ~30% smaller
        std::string nad_summ = p.summ;
        nad_summ.replace(nad_summ.find("_ad_"), 4, "_nad_");

        write(gyre / sfmt("gyre_ad_%s.in", p.name),
              hdr + p.note + "&constants\n/\n" + model_blk + modes + osc("") + "\n&rot\n/\n"
              + "\n&num\n  diff_scheme = 'COLLOC_GL4'\n/\n" + scans + grids
              + sfmt("\n&ad_output\n  summary_file = '%s'\n"
                     "  detail_template = 'detail_wide/detail.l%%l.n%%n.h5'\n"
                     "  summary_item_list = '%s'\n  detail_item_list = '%s'\n"
                     "  freq_units = 'CYC_PER_DAY'\n/\n\n&nad_output\n/\n",
                     p.summ, sum_items.c_str(), det_items.c_str()), force);

        std::string note2 = p.note;
        note2.replace(note2.find("--"), 2, "-- nonadiabatic,");
        write(gyre / sfmt("gyre_nad_%s.in", p.name),
              hdr + note2 + "&constants\n/\n" + model_blk + modes
              + osc("  nonadiabatic = .TRUE.\n") + "\n&rot\n/\n"
              + "\n&num\n  diff_scheme = 'MAGNUS_GL2'\n  nad_search = 'AD'\n/\n"
              + scans + grids + "\n&ad_output\n/\n"
              + sfmt("\n&nad_output\n  summary_file = '%s'\n  summary_item_list = '%s'\n"
                     "  freq_units = 'CYC_PER_DAY'\n/\n",
                     nad_summ.c_str(), nad_items.c_str()), force);
    }
    return 0;
}
