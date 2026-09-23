#include "model.hpp"
#include "numeric.hpp"

#include <highfive/H5File.hpp>

#include <fstream>
#include <numeric>
#include <regex>
#include <sstream>

namespace H5 = HighFive;

// GYRE writes every eigenfunction as an HDF5 compound {re, im}, adiabatic runs
// included. HighFive needs the layout spelled out once.
static H5::CompoundType create_ReIm() {
    return {{"re", H5::create_datatype<double>()},
            {"im", H5::create_datatype<double>()}};
}
HIGHFIVE_REGISTER_TYPE(ReIm, create_ReIm)


/* ---------------
    Raw files I/O
   ---------------
*/

static eig::ArrayXd read_real (const H5::File& f, const std::string& name) {
    auto v = f.getDataSet(name).read<std::vector<double>>();
    return eig::Map<eig::ArrayXd>(v.data(), v.size());
}

static eig::ArrayXcd read_complex (const H5::File& f, const std::string& name) {
    auto raw = f.getDataSet(name).read<std::vector<ReIm>>();
    eig::ArrayXcd out(raw.size());
    for (size_t i=0 ; i<raw.size() ; i++)
        out[i] = {raw[i].re, raw[i].im};
    return out;
}

DetailData read_detail (const std::filesystem::path& path) {
    H5::File f(path.string(), H5::File::ReadOnly);
    DetailData d;

    d.x   = read_real(f, "x");        d.rho = read_real(f, "rho");
    d.P   = read_real(f, "P");        d.Gamma_1 = read_real(f, "Gamma_1");
    d.c_1 = read_real(f, "c_1");      d.V_2 = read_real(f, "V_2");
    d.As  = read_real(f, "As");       d.M_r = read_real(f, "M_r");

    d.xi_r     = read_complex(f, "xi_r");
    d.xi_h     = read_complex(f, "xi_h");
    d.eul_Phi  = read_complex(f, "eul_Phi");
    d.deul_Phi = read_complex(f, "deul_Phi");
    d.has_lag_rho = f.exist("lag_rho");
    if (d.has_lag_rho) d.lag_rho = read_complex(f, "lag_rho");

    f.getAttribute("R_star").read(d.R_star);
    f.getAttribute("M_star").read(d.M_star);
    f.getAttribute("L_star").read(d.L_star);
    f.getAttribute("l").read(d.l);
    f.getAttribute("n_pg").read(d.n_pg);
    ReIm w; f.getAttribute("omega").read(w);  d.omega = {w.re, w.im};

    return d;
}

SummaryData read_summary (const std::filesystem::path& path) {
    H5::File f(path.string(), H5::File::ReadOnly);
    SummaryData s;

    s.l = f.getDataSet("l").read<std::vector<int>>();
    auto w = f.getDataSet("omega").read<std::vector<ReIm>>();
    s.omega.resize(w.size());
    for (size_t i = 0; i < w.size(); i++) s.omega[i] = {w[i].re, w[i].im};

    return s;
}


/* --------------
    MESA Profile
   --------------
*/

// 5 header lines, column names on line 6, whitespace-delimited rows after.
static std::pair<std::map<std::string,size_t>, std::vector<std::vector<double>>>
read_mesa_profile (const std::filesystem::path& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open MESA profile " + path.string());

    std::string line;
    for (int i=0 ; i<5 ; i++) std::getline(f, line);
    std::getline(f, line);

    std::map<std::string, size_t> names;
    {
        std::istringstream iss(line);
        std::string tok;
        size_t i = 0;
        while (iss >> tok) names[tok] = i++;
    }

    std::vector<std::vector<double>> rows;
    while (std::getline(f, line)) {
        std::istringstream iss(line);
        std::vector<double> row;
        double v;
        while (iss >> v) row.push_back(v);
        if (!row.empty()) rows.push_back(std::move(row));
    }
    return {names, rows};
}

// Matched in log P, which is monotonic in both grids and tracks ionisation.
std::pair<std::map<std::string, eig::ArrayXd>, int>
mesa_on_gyre_grid (const eig::ArrayXd& P, const std::filesystem::path& profile) {
    auto [names, rows] = read_mesa_profile(profile);
    size_t lp_col = names.at("logP");

    std::vector<size_t> order(rows.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(
        order.begin(), order.end(),
        [&](size_t a, size_t b) { return rows[a][lp_col] < rows[b][lp_col]; }
    );

    eig::ArrayXd lp_mesa(rows.size());
    for (size_t i=0 ; i<order.size() ; i++) lp_mesa[i] = rows[order[i]][lp_col];

    // GYRE's atmosphere points sit below MESA's surface pressure. num::interp
    // clamps them to the outermost MESA value rather than extrapolating into
    // the H ionisation zone, where dGamma1_dlnRho_s is largest.
    eig::ArrayXd lp = P.log10();

    static const std::array<const char*, 6> needed = {
        "dGamma1_dlnRho_s", "temperature", "cp", "luminosity", "conv_vel", "omega_conv"};

    std::map<std::string, eig::ArrayXd> out;
    for (auto c : needed) {
        size_t col = names.at(c);
        eig::ArrayXd y(rows.size());
        for (size_t i=0 ; i<order.size() ; i++) y[i] = rows[order[i]][col];
        out[c] = num::interp(lp, lp_mesa, y);
    }

    int n_clamped = int((lp < lp_mesa[0]).count());
    return {out, n_clamped};
}


/* -----------
    Structure
   -----------
*/

Star load_star (const DetailData& d,
                const std::optional<std::filesystem::path>& mesa_profile) {
    double R = d.R_star, M = d.M_star;
    eig::ArrayXd x = d.x, r = x * R;
    eig::ArrayXd rho = d.rho, G1 = d.Gamma_1;

    // c_1 = (r/R)^3 (M/M_r), so g = G M x / (c_1 R^2) stays finite at x = 0.
    eig::ArrayXd g = G * M * x / (d.c_1 * R * R);

    // r > 0: dg/dr = 4 pi G rho - 2 g / r.   r = 0: dg/dr = (4/3) pi G rho.
    eig::ArrayXd dg_dr = (r > 0).select(
        4.0 * M_PI * G * rho - 2.0 * g / r.max(1e-300),
        (4.0 / 3.0) * M_PI * G * rho);

    eig::ArrayXd dlnrho_dlnr = -d.V_2 * x.square() / G1 - d.As;

    eig::ArrayXd dGamma1_dlnrho_s, T, c_P, L_r, v_conv, omega_conv;
    int n_clamped = 0;

    if (mesa_profile) {
        auto [cols, nc] = mesa_on_gyre_grid(d.P, *mesa_profile);
        dGamma1_dlnrho_s = cols.at("dGamma1_dlnRho_s");
        T          = cols.at("temperature");
        c_P        = cols.at("cp");
        L_r        = cols.at("luminosity") * L_SUN;
        v_conv     = cols.at("conv_vel");
        omega_conv = cols.at("omega_conv");
        n_clamped  = nc;
    } else {
        // Polytrope path (testing): no MESA profile, so no thermal structure
        // dGamma1_dlnrho_s = 0 rather than NaN
        dGamma1_dlnrho_s = eig::ArrayXd::Zero(x.size());
        T = c_P = L_r = v_conv = omega_conv =
            eig::ArrayXd::Constant(x.size(), std::nan(""));
    }

    return Star{x, r, rho, d.P, G1, g, dg_dr, dlnrho_dlnr, dGamma1_dlnrho_s,
                d.M_r, T, c_P, L_r, v_conv, omega_conv,
                M, R, d.L_star, n_clamped};
}

bool Star::has_thermal_structure () const {
    return c_P.isFinite().all() && L > 0.0;
}

// Time to radiate the heat content above r. Sets where the adiabatic
// eigenfunctions stop being trustworthy: omega t_thermal < 2 pi.
eig::ArrayXd Star::t_thermal () const {
    if (!has_thermal_structure())
        throw std::runtime_error("no thermal structure in this model (polytrope?)");

    eig::ArrayXd u = c_P * T * rho * 4.0 * M_PI * r.square();
    eig::ArrayXd above(r.size());
    above[0] = 0.0;
    for (eig::Index i = 1; i < r.size(); i++)
        above[i] = above[i-1] + 0.5 * (u[i] + u[i-1]) * (r[i] - r[i-1]);

    return (above[above.size()-1] - above) / L;
}


/* -------
    Modes
   -------
*/

// int rho (xi_r^2 + Lambda^2 xi_h^2) r^2 dr, no 4 pi.
static double inertia (const Star& s, int l,
                       const eig::ArrayXd& xi_r, const eig::ArrayXd& xi_h) {
    return num::trapz(s.rho * (xi_r.square() + l*(l+1) * xi_h.square()) * s.r.square(), s.r);
}

double Eigenfunction::inertia_integral () const {
    return inertia(*starptr, l, xi_r, xi_h);
}

Eigenfunction load_one_ef (const std::filesystem::path& detail_file,
                        const std::optional<std::filesystem::path>& mesa_profile,
                        const DampingRates* gammas, GridCache& cache) {
    auto d = read_detail(detail_file);

    std::vector<double> key(d.x.data(), d.x.data() + d.x.size());
    auto it = cache.find(key);
    if (it == cache.end())
        it = cache.emplace(key, std::make_shared<const Star>(load_star(d, mesa_profile))).first;
    auto star = it->second;

    double omega_dimless = d.omega.real();
    double omega = omega_dimless * star->omega_dyn();
    double gm_r = G * star->M / star->R;

    eig::ArrayXd xi_r = d.xi_r.real() * star->R;
    eig::ArrayXd xi_h = d.xi_h.real() * star->R;
    eig::ArrayXd delta_Phi     = d.eul_Phi.real()  * gm_r;
    eig::ArrayXd ddelta_Phi_dr = d.deul_Phi.real() * gm_r / star->R;

    eig::ArrayXd div_xi;
    if (d.has_lag_rho) {
        // delta rho = -rho div.xi, GYRE reports lag_rho in units of rho.
        div_xi = -d.lag_rho.real();
    } else if (d.l > 0) {
        // Horizontal momentum equation. Undefined for l = 0, where xi_h = 0.
        div_xi = star->rho * (star->g * xi_r - omega*omega * star->r * xi_h + delta_Phi)
                 / (star->Gamma_1 * star->P);
    } else {
        throw std::runtime_error(detail_file.filename().string() +
            ": l = 0 needs lag_rho (xi_h is identically zero, so the horizontal "
            "Euler route does not exist). Re-run gyre_ad.in.");
    }

    // A = sqrt(E_star / (2 omega^2 I))
    double A = std::sqrt(star->E_star() / (2.0 * omega*omega * inertia(*star, d.l, xi_r, xi_h)));

    return Eigenfunction{
        d.l, d.n_pg, omega, omega_dimless,
        gammas ? (*gammas)(d.l, omega) : std::nan(""),
        xi_r * A, xi_h * A, div_xi * A, delta_Phi * A, ddelta_Phi_dr * A,
        star,
    };
}

std::vector<Mode> build_mode_list (const ModeMap& eigenfuncs,
                                   std::optional<int> l_max,
                                   std::optional<std::pair<int,int>> n_range,
                                   bool require_gamma) {
    std::vector<Mode> modes;
    // ModeMap already iterates sorted by (l, n_pg).
    for (auto& [key, ef] : eigenfuncs) {
        auto [l, n_pg] = key;
        if (l_max && l > *l_max) continue;
        if (n_range && !(n_range->first <= n_pg && n_pg <= n_range->second)) continue;
        // Drops modes the nonadiabatic run missed; a NaN gamma would otherwise
        // propagate silently into mu and E_th.
        if (require_gamma && !std::isfinite(ef.gamma)) continue;

        for (int m = -l; m <= l; m++)
            for (int s : {+1, -1})
                modes.push_back(Mode{n_pg, l, m, s, &ef});
    }
    return modes;
}


/* ---------------
    Damping Rates
   ---------------
*/

double DampingRates::operator()(int l, double q) const {
    auto it = omega.find(l);
    if (it == omega.end() || it->second.empty())
        return std::nan("");

    const std::vector<double>& w = it->second;

    // nearest index to q, snapping to whichever neighbor is closer
    size_t j = std::lower_bound(w.begin(), w.end(), q) - w.begin();
    if ( j>0 && (j==w.size() || q-w[j-1] < w[j]-q) ) j--;
    j = std::min(j, w.size()-1);

    std::vector<double> gaps;   // local mode spacing
    if (j>0) gaps.push_back(w[j] - w[j-1]);
    if (j+1 < w.size()) gaps.push_back(w[j+1] - w[j]);

    double spacing = gaps.empty() ?
        std::abs(q) : *std::min_element(gaps.begin(), gaps.end());

    return std::abs(w[j] - q) <= tol * spacing ? gamma.at(l)[j] : std::nan("");
}

DampingRates load_gammas (const std::vector<std::filesystem::path>& nad, double omega_dyn) {
    std::map<int, std::vector<std::pair<double,double>>> per_l;   // l -> (omega, gamma)

    for (auto& path : nad) {
        auto d = read_summary(path);
        for (size_t i = 0; i < d.l.size(); i++)
            per_l[d.l[i]].push_back({ d.omega[i].real() * omega_dyn,
                              -d.omega[i].imag() * omega_dyn });   // gamma = -Im(omega)
    }

    DampingRates out;
    for (auto& [l, pairs] : per_l) {
        // operator() binary-searches omega, so this sort is load-bearing.
        std::sort(pairs.begin(), pairs.end());
        for (auto& [w, g] : pairs) { out.omega[l].push_back(w); out.gamma[l].push_back(g); }
    }
    return out;
}


/* -------------------
    Turbulent Damping
   -------------------
*/

// Effective turbulent viscosity, cm^2/s, from Duguid et al. (2020) eq. 13.
// nu = u_mlt l_mlt f(w), w = |omega| / omega_c, omega_c = u_mlt / l_mlt, so
// l_mlt = v_conv / omega_conv and alpha_MLT never appears -- MESA's omega_conv
// already carries it. Zero outside the convection zones, NaN with no MESA profile.
eig::ArrayXd nu_turb (const Star& s, double omega) {
    auto cz = (s.v_conv > 0) && (s.omega_conv > 0);

    eig::ArrayXd denom = cz.select(s.omega_conv, 1.0);
    eig::ArrayXd w = cz.select(std::abs(omega) / denom, 1.0);

    // f is continuous at both breakpoints.
    eig::ArrayXd f = (w < 1e-2).select(
        5.0,
        (w <= 5.0).select(0.5 * w.pow(-0.5),
                          0.5 * std::pow(5.0, 1.5) * w.pow(-2.0)));

    eig::ArrayXd l_mlt = cz.select(s.v_conv / denom, 0.0);

    // else branch is 0 where v_conv is finite and NaN where it is not.
    return cz.select(s.v_conv * l_mlt * f, 0.0 * s.v_conv);
}

// Higgins & Kopal (1968) dissipation function in the form of Lai (1994) eq. 8.8,
// for xi = xi_r Y e_r + xi_h r grad Y with int |Y_lm|^2 dOmega = 1.
// F = 2 eps_ij eps_ij - (2/3)(div xi)^2 after the angular integral, so it is
// positive semi-definite. eps_rr = div xi - trace_h exactly, so the radial
// pieces cancel with no second estimate of dxi_r/dr to spoil it.
// No m dependence: Lai's (l+|m|)!/(l-|m|)! factor is an artefact of his
// Jackson-normalised Y_lm.
static eig::ArrayXd viscous_F (int l, const eig::ArrayXd& r, const eig::ArrayXd& xi_r,
                               const eig::ArrayXd& xi_h, const eig::ArrayXd& div_xi) {
    double L2 = double(l * (l + 1));
    eig::ArrayXd trace_h = 2.0 * xi_r / r - L2 * xi_h / r;
    eig::ArrayXd dxi_r = div_xi - trace_h;

    // GYRE gives dxi_r, but not dxi_h. Hence get derivative from cubic spline
    eig::ArrayXd dxi_h = num::CubicSpline(r, xi_h).deriv_at_knots();

    return 2.0 * dxi_r.square()
         + trace_h.square()
         + L2 * (dxi_h + (xi_r - xi_h) / r).square()
         + (l - 1) * L2 * (l + 2) * (xi_h / r).square()
         - (2.0 / 3.0) * div_xi.square();
}

// Turbulent damping rate, rad/s, from MW23 eq. 11:
//     gamma = (omega^2 / E) int dr rho r^2 nu_turb F(r),
// with E = E_star by the normalisation.
double gamma_turb (const Eigenfunction& ef) {
    const Star& s = *ef.starptr;

    // r > 0 && x <= 1
    eig::Index lo = 0;
    while (lo < s.r.size() && !(s.r[lo] > 0.0)) lo++;
    eig::Index hi = std::upper_bound(s.x.data(), s.x.data() + s.x.size(), 1.0) - s.x.data();
    eig::Index n = hi - lo;
    if (n <= 1) return 0.0;

    eig::ArrayXd r = s.r.segment(lo, n);
    eig::ArrayXd F = viscous_F(ef.l, r, ef.xi_r.segment(lo, n),
                               ef.xi_h.segment(lo, n), ef.div_xi.segment(lo, n));

    eig::ArrayXd integrand = s.rho.segment(lo, n)
                           * nu_turb(s, ef.omega).segment(lo, n)
                           * r.square() * F;

    return ef.omega * ef.omega * num::trapz(integrand, r) / s.E_star();
}


/* ---------------
    Model Loading
   ---------------
*/

static bool matches_detail_pattern (const std::filesystem::path& p) {
    static const std::regex re(R"(^detail\.l\d+\.n[+-]\d+\.h5$)");
    return std::regex_match(p.filename().string(), re);
}

// Get profileN.data if exist; nullopt if not (cases like test polytropes)
std::optional<std::filesystem::path> paired_mesa_profile (
        const std::filesystem::path& inlist, const std::filesystem::path& mesa_path) {
    std::ifstream f(inlist);
    if (!f) throw std::runtime_error("cannot open inlist " + inlist.string());

    std::stringstream ss; ss << f.rdbuf();
    std::string text = ss.str();

    std::smatch m;
    static const std::regex re(R"(^[ \t]*file[ \t]*=[ \t]*'([^']+)')",
                               std::regex::multiline);
    if (!std::regex_search(text, m, re))
        throw std::runtime_error("no model file in " + inlist.string());

    std::string name = std::filesystem::path(m[1].str()).filename().string();
    const std::string ext = ".GYRE";
    if (name.size() <= ext.size() || name.compare(name.size()-ext.size(), ext.size(), ext) != 0)
        return std::nullopt;

    name.erase(name.size() - ext.size());
    return mesa_path / "LOGS" / name;
}

Model load_model (const std::filesystem::path& model_path,
                  const std::string& detail_path, const std::string& inlist,
                  const std::vector<std::string>& nad) {
    auto gyre = model_path / "gyre", mesa = model_path / "mesa";
    auto detail = gyre / detail_path;

    std::vector<std::filesystem::path> files;
    for (auto& e : std::filesystem::directory_iterator(detail))
        if (matches_detail_pattern(e.path()))
            files.push_back(e.path());
    std::sort(files.begin(), files.end());
    if (files.empty()) throw std::runtime_error("no detail files in " + detail.string());

    auto profile = paired_mesa_profile(gyre / inlist, mesa);
    auto seed = read_detail(files.front());
    auto star = std::make_shared<const Star>(load_star(seed, profile));

    std::vector<std::filesystem::path> found;
    for (auto& n : nad)
        if (std::filesystem::exists(gyre / n))
            found.push_back(gyre / n);
    std::optional<DampingRates> gammas;
    if (!found.empty()) gammas = load_gammas(found, star->omega_dyn());

    ModeMap eigenfuncs;
    GridCache grid_cache;
    grid_cache[std::vector<double>(star->x.data(), star->x.data()+star->x.size())] = star;
    for (auto& f : files) {
        auto ef = load_one_ef(f, profile, gammas ? &*gammas : nullptr, grid_cache);
        eigenfuncs[{ef.l, ef.n_pg}] = std::move(ef);
    }

    return {star, std::move(eigenfuncs)};
}
