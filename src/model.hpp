#pragma once
#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace eig = Eigen;

constexpr double G = 6.67430e-8;
constexpr double L_SUN = 3.828e33;
constexpr double CD = 7.27220522e-5;        // rad/s per cycle/day


/* ---------------
    Raw files I/O
   ---------------
*/

struct ReIm { double re, im; };

struct DetailData {                                      
    eig::ArrayXd x, rho, P, Gamma_1, c_1, V_2, As, M_r;  
    eig::ArrayXcd xi_r, xi_h, eul_Phi, deul_Phi, lag_rho;
    double R_star, M_star, L_star;                       
    int l, n_pg;                                         
    std::complex<double> omega;                          
    bool has_lag_rho;                                    
};

struct SummaryData {
    std::vector<int> l;
    std::vector<std::complex<double>> omega;
};

DetailData read_detail (const std::filesystem::path& path);
SummaryData read_summary (const std::filesystem::path& path);


/* -----------
    Structure
   -----------
*/

struct Star {
    eig::ArrayXd x, r, rho, P, Gamma_1, g, dg_dr, dlnrho_dlnr,
                    dGamma1_dlnrho_s, M_r, T, c_P, L_r, v_conv, omega_conv;
    double M, R, L;
    int n_atm_clamped;  // number of points in the atmosphere that are clamped to the surface values. Just for checking the model, not used in calculations.

    double E_star () const { return G * M * M / R; }
    double omega_dyn () const {
        return std::sqrt( G*M / (R*R*R) );
    }
    eig::ArrayXd t_thermal () const;
    bool has_thermal_structure () const;
    eig::Index size () const { return x.size(); }
};

// dGamma1_dlnRho_s, temperature, cp, luminosity, conv_vel, omega_conv on the
// GYRE grid, matched in log P; second member is n_atm_clamped.
std::pair<std::map<std::string, eig::ArrayXd>, int>
mesa_on_gyre_grid (
    const eig::ArrayXd& P,
    const std::filesystem::path& profile
);

Star load_star (
    const DetailData& d,
    const std::optional<std::filesystem::path>& mesa_profile
);


/* -------
    Modes
   -------
*/

struct Eigenfunction {
    int l, n_pg;
    double omega, omega_dimless, gamma;     // rad/s
    eig::ArrayXd xi_r, xi_h, div_xi, delta_Phi, ddelta_Phi_dr;
    std::shared_ptr<const Star> starptr;

    double Lambda2 () const { return double(l * (l + 1)); }
    double inertia_integral () const;
};

using ModeMap = std::map<std::pair<int,int>, Eigenfunction>; // <(l, n_pg), Eigenfunction>

struct Mode {
    int n_pg, l, m, s;
    const Eigenfunction* ef;
};

struct DampingRates {
    std::map<int, std::vector<double>> omega, gamma;    // <l, []>, both sorted by omega
    double tol = 0.3;
    double operator()(int l, double q) const;           // return NaN if no match
};

struct Model {
    std::shared_ptr<const Star> star;
    ModeMap eigfuncs;
};

using GridCache = std::map<std::vector<double>, std::shared_ptr<const Star>>;


/* ---------
    Loaders
   ---------
*/

Model load_model (
    const std::filesystem::path& model_path,
    const std::string& detail_path = "detail",
    const std::string& inlist = "gyre_ad.in",
    const std::vector<std::string>& nad = {"summary_nad.h5"}
);

Eigenfunction load_one_ef (
    const std::filesystem::path& detail_file,
    const std::optional<std::filesystem::path>& mesa_profile,
    const DampingRates* gammas,     // nullptr if no nad found
    GridCache& cache
);

DampingRates load_gammas (
    const std::vector<std::filesystem::path>& nad,
    double omega_dyn
);

std::optional<std::filesystem::path> paired_mesa_profile (
    const std::filesystem::path& inlist,
    const std::filesystem::path& mesa_path
);

std::vector<Mode> build_mode_list (
    const ModeMap& eigenfuncs,
    std::optional<int> l_max = {},
    std::optional<std::pair<int,int>> n_range = {},
    bool require_gamma = false
);


/* -------------------
    Turbulent Damping
   -------------------
*/

eig::ArrayXd nu_turb (const Star& s, double omega);      // cm^2/s, Duguid+2020 eq. 13
double gamma_turb (const Eigenfunction& ef);             // rad/s, MW23 eq. 11, always > 0
