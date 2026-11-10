/* gamma_rad and gamma_turb (MW23 eq. 11) per mode, for MW23 Fig. 4.

       ./damping --model models/dsct_M2.0_T7696
       python3 scripts/damping_plot.py --data out/damping_dsct_M2.0_T7696.csv
*/

#include "csv.hpp"
#include "model.hpp"

#include <CLI/CLI.hpp>

#include <cstdio>

int main (int argc, char** argv) {
    CLI::App app{"gamma_rad and gamma_turb per mode"};
    std::string model = "models/dsct_M2.0_T7696", detail_dir = "detail_wide", inlist = "gyre_ad_wide.in";
    std::vector<std::string> nad = {"summary_nad_wide.h5", "summary_nad_wide_hi.h5"};
    std::string out;
    app.add_option("--model", model, "default models/dsct_M2.0_T7696");
    app.add_option("--detail-dir", detail_dir, "default detail_wide");
    app.add_option("--inlist", inlist, "default gyre_ad_wide.in");
    app.add_option("--nad", nad, "default summary_nad_wide.h5 summary_nad_wide_hi.h5")->expected(-1);
    app.add_option("--out", out, "default out/damping_<model tag>.csv");
    CLI11_PARSE(app, argc, argv);
    if (out.empty()) out = "out/damping_" + std::filesystem::path(model).filename().string() + ".csv";

    const Model m = load_model(model, detail_dir, inlist, nad);
    csv::Writer w(out, {"l", "n", "omega", "f", "gamma_rad", "gamma_turb"});
    for (const auto& [k, e] : m.eigfuncs) w.row(e.l, e.n_pg, e.omega, e.omega / CD, e.gamma, gamma_turb(e));
    std::printf("%zu modes -> %s\n", m.eigfuncs.size(), out.c_str());
    return 0;
}
