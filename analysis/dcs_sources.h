#ifndef DCS_SOURCES_H
#define DCS_SOURCES_H
#include <Rtypes.h>
#include <string>
#include <vector>
#include "style.h"

using namespace PetroffPalette;

namespace dcs {

struct SrcStyle {
    const char* tag;
    Color_t     color;
    Style_t     line;
    SrcStyle(const char* t, Color_t c, Style_t l) : tag(t), color(c), line(l) {}
};


// Single source-of-truth for the 5 sources + their plotting styles.
// inline const std::vector<SrcStyle>& Sources() {
//     static const std::vector<SrcStyle> s = {
//         { "gamma_662keV",  kRed     + 1, 1 },
//         { "gamma_1173keV", kGreen   + 2, 1 },
//         { "gamma_1332keV", kAzure   + 1, 1 },
//         { "isotope_Eu154", kMagenta + 1, 1 },
//         { "neutron_Watt",  kOrange  + 7, 1 },
//     };
//     return s;
// }
inline const std::vector<SrcStyle>& Sources() {
    static const std::vector<SrcStyle> s = {
        { "gamma_662keV",  getP6(0), 1 },
        { "gamma_1173keV", getP6(1), 1 },
        { "gamma_1332keV", getP6(2), 1 },
        { "isotope_Eu154", getP6(3), 1 },
        { "neutron_Watt",  getP6(4), 1 },
    };
    return s;
}

inline const char* PrettyLabel(const std::string& tag) {
    if (tag == "gamma_662keV")  return "#gamma 662 keV (^{137}Cs)";
    if (tag == "gamma_1173keV") return "#gamma 1173 keV (^{60}Co)";
    if (tag == "gamma_1332keV") return "#gamma 1332 keV (^{60}Co)";
    if (tag == "isotope_Eu154") return "^{154}Eu";
    if (tag == "neutron_Watt")  return "neutron (Watt)";
    return tag.c_str();
}

} // namespace dcs
#endif

