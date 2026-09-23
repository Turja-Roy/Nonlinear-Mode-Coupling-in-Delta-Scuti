#include "csv.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <stdexcept>

namespace csv {

std::string fmt (double v) {
    if (std::isnan(v)) return "";
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return buf;
}

std::string fmt (bool v) { return v ? "True" : "False"; }

Writer::Writer (const std::filesystem::path& p, const std::vector<std::string>& header,
                bool force) {
    if (!force && std::filesystem::exists(p))
        throw std::runtime_error("refusing to overwrite " + p.string() + " (pass --force)");
    std::filesystem::create_directories(p.parent_path());
    f_.open(p);
    if (!f_) throw std::runtime_error("cannot write " + p.string());
    for (size_t i=0 ; i<header.size() ; i++) f_ << (i ? "," : "") << header[i];
    f_ << '\n';
}

const std::vector<std::string>& Table::text (const std::string& name) const {
    auto it = text_.find(name);
    if (it == text_.end()) throw std::runtime_error("no column '" + name + "'");
    return it->second;
}

const std::vector<double>& Table::numbers (const std::string& name) const {
    auto cached = num_.find(name);
    if (cached != num_.end()) return cached->second;

    std::vector<double> v;
    v.reserve(rows());
    for (const std::string& s : text(name)) {
        if (s.empty()) { v.push_back(std::nan("")); continue; }   // pandas NaN
        size_t used = 0;
        double x = 0.0;
        try { x = std::stod(s, &used); } catch (...) { used = 0; }
        if (used != s.size())
            throw std::runtime_error("column '" + name + "' is not numeric: '" + s + "'");
        v.push_back(x);
    }
    return num_.emplace(name, std::move(v)).first->second;
}

Table read (const std::filesystem::path& p) {
    std::ifstream f(p);
    if (!f) throw std::runtime_error("cannot read " + p.string());

    Table out;
    std::string line, cell;
    if (!std::getline(f, line)) return out;
    std::istringstream head(line);
    while (std::getline(head, cell, ',')) out.header_.push_back(cell);

    std::vector<std::vector<std::string>> cols(out.header_.size());
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        for (size_t i=0 ; i<cols.size() ; i++) {
            if (!std::getline(ss, cell, ',')) cell.clear();
            cols[i].push_back(cell);
        }
    }
    for (size_t i=0 ; i<cols.size() ; i++) out.text_[out.header_[i]] = std::move(cols[i]);
    return out;
}

double quantile (std::vector<double> v, double q) {
    v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return std::isnan(x); }),
            v.end());
    if (v.empty()) return std::nan("");
    std::sort(v.begin(), v.end());

    const double pos = q * double(v.size() - 1);
    const size_t lo = size_t(std::floor(pos)), hi = size_t(std::ceil(pos));
    return v[lo] + (v[hi] - v[lo]) * (pos - double(lo));
}

double median (std::vector<double> v) { return quantile(std::move(v), 0.5); }

}  // namespace csv
