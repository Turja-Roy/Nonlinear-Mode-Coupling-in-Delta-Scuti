#pragma once
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <vector>

/* Minimal CSV I/O, pandas-compatible on the two things that bite: NaN is
   written as an empty field, and bool as True/False. */

namespace csv {

std::string fmt (double v);                      // NaN -> "", inf -> "inf"
std::string fmt (bool v);                        // "True" / "False"
inline std::string fmt (int v)  { return std::to_string(v); }
inline std::string fmt (long v) { return std::to_string(v); }
inline std::string fmt (const std::string& v) { return v; }
inline std::string fmt (const char* v) { return v; }

class Writer {
public:
    Writer (const std::filesystem::path& p, const std::vector<std::string>& header,
            bool force = true);

    template <class... T>
    void row (const T&... vals) {
        std::string line;
        int i = 0;
        ((line += (i++ ? "," : "") + fmt(vals)), ...);
        f_ << line << '\n';
    }

private:
    std::ofstream f_;
};

/* Every column is read as text. numbers() parses one column to double on
   first use and caches it, so a numeric column costs one parse and a string
   column is never parsed at all. */
class Table {
public:
    const std::vector<std::string>& header () const { return header_; }
    size_t rows () const { return header_.empty() ? 0 : text_.at(header_[0]).size(); }
    bool has (const std::string& name) const { return text_.count(name) != 0; }

    const std::vector<std::string>& text (const std::string& name) const;
    // Throws if any field fails to parse; an empty field is NaN, as in pandas.
    const std::vector<double>& numbers (const std::string& name) const;

    friend Table read (const std::filesystem::path& p);

private:
    std::vector<std::string> header_;
    std::map<std::string, std::vector<std::string>> text_;
    mutable std::map<std::string, std::vector<double>> num_;
};

Table read (const std::filesystem::path& p);

// Linear interpolation between order statistics, as numpy does it. NaNs dropped.
double quantile (std::vector<double> v, double q);
double median (std::vector<double> v);

}  // namespace csv
