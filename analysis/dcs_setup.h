#ifndef DCS_SETUP_H
#define DCS_SETUP_H
#include <vector>
#include <string>
#include <cmath>
#include <Rtypes.h>          // Color_t, kBlack, kRed, ...

namespace dcs {

// ---------------------------------------------------------------------------
//  Detector (global frame): faces the central-cluster target cask 0 along -x.
//  From generate_macros.sh: /dcs-monitor/det/setPosition 3025 1580 0 mm,
//                           /dcs-monitor/det/setRotation  0 0 -90.
// ---------------------------------------------------------------------------
constexpr double kDetX = 3025.0, kDetY = 1580.0, kDetZ = 0.0;

struct Point { double x, y; };

// ---------------------------------------------------------------------------
//  12 casks = three side-facing 2x2 CASTOR-440 clusters, symmetric about the
//  y-axis, inside the concrete storage hall. Add-order == caskNum:
//     central 0-3, +x 4-7, -x 8-11   (see generate_macros.sh).
//  Cask 0 (central, +x/+y) is the detector target -> SIGNAL.
// ---------------------------------------------------------------------------
constexpr int kNCask = 12;

inline const Point kCaskPos[kNCask] = {
    {  1580.0,  1580.0 },   // 0  central  (detector target = signal)
    {  1580.0, -1580.0 },   // 1  central
    { -1580.0,  1580.0 },   // 2  central
    { -1580.0, -1580.0 },   // 3  central
    {  6240.0,  1580.0 },   // 4  +x
    {  6240.0, -1580.0 },   // 5  +x
    {  9400.0,  1580.0 },   // 6  +x
    {  9400.0, -1580.0 },   // 7  +x
    { -6240.0,  1580.0 },   // 8  -x
    { -6240.0, -1580.0 },   // 9  -x
    { -9400.0,  1580.0 },   // 10 -x
    { -9400.0, -1580.0 },   // 11 -x
};

constexpr double kCaskR = 1270.0;   // emission-shell radius [mm]

// Spatially-faithful pad map for TCanvas::Divide(kPadCols,kPadRows):
//   column increases with x (left->right); top row = +y, bottom row = -y.
//   x-columns: -9400,-6240,-1580,1580,6240,9400  ->  cols 1..6
constexpr int kPadCols = 6, kPadRows = 2;
inline const int kCaskPad[kNCask] = { 4, 10, 3, 9, 5, 11, 6, 12, 2, 8, 1, 7 };

// ---- Signal / background policy (single switch point) ----------------------
//  Current: signal = cask 0, background = all 11 remaining casks.
inline bool IsSignalCask(int ic)     { return ic == 0; }
inline bool IsBackgroundCask(int ic) { return ic != 0; }

// ---- Per-cask colour: cask 0 stands out in BLACK; the three clusters are
//      grouped by hue family (central=reds, +x=blues, -x=greens) so the eye
//      can tell clusters apart at a glance on 12-way overlays. --------------
inline Color_t CaskColor(int ic) {
    switch (ic) {
        case 0:  return kBlack;        // signal target -- deliberately distinct
        // central cluster (1-3): reds/pinks
        case 1:  return kRed + 1;
        case 2:  return kRed - 4;
        case 3:  return kPink + 7;
        // +x cluster (4-7): blues/violets
        case 4:  return kAzure + 2;
        case 5:  return kAzure - 4;
        case 6:  return kBlue  - 7;
        case 7:  return kViolet + 6;
        // -x cluster (8-11): greens/teals
        case 8:  return kGreen + 2;
        case 9:  return kSpring - 6;
        case 10: return kTeal + 3;
        case 11: return kGreen - 9;
        default: return kGray + 1;
    }
}
inline Width_t CaskLineWidth(int ic) { return ic == 0 ? 2 : 1; }

// ===========================================================================
//  Config grammar / MakeTitle (retained for compatibility) — unchanged.
// ===========================================================================
namespace detail {
inline std::vector<std::string> Split(const std::string& s, char d) {
    std::vector<std::string> out; std::string cur;
    for (char c : s) { if (c == d) { out.push_back(cur); cur.clear(); } else cur += c; }
    out.push_back(cur);
    return out;
}
inline std::string PctFromPE(const std::string& tok) {
    std::string t = tok;
    if (t.rfind("PE", 0) == 0) t = t.substr(2);
    if (!t.empty() && t.back() == 'B') t.pop_back();
    for (auto& c : t) if (c == 'p') c = '.';
    double v = 0.0;
    try { v = t.empty() ? 0.0 : std::stod(t); } catch (...) { v = 0.0; }
    return std::to_string((int)std::lround(v * 100.0));
}
inline void MatBore(const std::string& tok, std::string& mat, std::string& bore) {
    auto dash = Split(tok, '-');
    std::string m = dash[0];
    size_t cm = m.find("cm");
    mat  = (cm == std::string::npos) ? m : m.substr(cm + 2);
    bore = "";
    if (dash.size() > 1)
        bore = (dash[1].rfind("0cm", 0) == 0) ? ", no bore" : "+bore";
}
} // namespace detail

inline std::string MakeTitle(const std::string& dir) {
    if (dir == "nominal") return "Nominal";
    auto tok = detail::Split(dir, '_');
    std::string geo    = (tok.size() > 1) ? tok[1] : "";
    std::string matTok = (tok.size() > 2) ? tok[2] : "";
    std::string peTok  = (tok.size() > 3) ? tok[3] : "";
    std::string geoS = (geo.rfind("hemiInflated", 0) == 0) ? "Infl"
                     : (geo == "hemi")                     ? "Hemi" : geo;
    std::string mat, bore; detail::MatBore(matTok, mat, bore);
    return geoS + " " + mat + bore + ", " + detail::PctFromPE(peTok) + "%B";
}

struct Config { std::string dir; std::string title; };

// Single simulated run: the 3-cluster hall with 2 cm Pb gamma layer.
// Add further 12-cask trees here to overlay/compare them.
inline const std::vector<Config>& Configs() {
    static const std::vector<Config> c = {
        { "roomReturn/hall_3cluster/nominal", "nominal" },
    };
    return c;
}

} // namespace dcs
#endif

