#ifndef DCS_MULTIDET_H
#define DCS_MULTIDET_H

#include <TH2F.h>

#include <cmath>
#include <string>

namespace dcs {

// -----------------------------------------------------------------------------
// Multidet geometry from generate_macros_multidet.sh:
//
//   for Z in -1050 -750 -450 -150 150 450 750 1050:
//     for I in 0..23:
//       PHI = 15*I
//       /dcs-monitor/det/clyc/add
//
// Therefore:
//   det = iz*24 + iphi
// -----------------------------------------------------------------------------
constexpr int    kNDetPhi       = 24;
constexpr int    kNDetZ         = 8;
constexpr int    kNMultidet     = kNDetPhi * kNDetZ;
constexpr double kDetPhiStepDeg = 15.0;

inline const double* MultidetZPlanes()
{
    static const double z[kNDetZ] = {
        -1050.0, -750.0, -450.0, -150.0,
          150.0,  450.0,  750.0, 1050.0
    };
    return z;
}

inline bool IsValidMultidetId(int det)
{
    return det >= 0 && det < kNMultidet;
}

inline int DetPhiIndex(int det)
{
    return IsValidMultidetId(det) ? (det % kNDetPhi) : 0;
}

inline int DetZIndex(int det)
{
    return IsValidMultidetId(det) ? (det / kNDetPhi) : 0;
}

inline double DetPhiDeg(int det)
{
    return kDetPhiStepDeg * DetPhiIndex(det);
}

inline double DetPhiRad(int det)
{
    return DetPhiDeg(det) * M_PI / 180.0;
}

inline double DetZMm(int det)
{
    return MultidetZPlanes()[DetZIndex(det)];
}

inline bool LooksLikeMultidetPath(const std::string& s)
{
    return s.find("nominal-all-positions") != std::string::npos ||
           s.find("all-positions")         != std::string::npos ||
           s.find("multidet")              != std::string::npos ||
           s.find("multi-det")             != std::string::npos;
}

// True detector-position heatmap:
//   x = detector phi [deg]
//   y = detector z [mm]
//
// Phi bin centers: 0, 15, ..., 345 deg.
// Z bin centers: -1050, -750, ..., 1050 mm.
inline TH2F* MakeDetectorPhiZHist(const char* name, const char* title)
{
    constexpr double phiLo = -0.5 * kDetPhiStepDeg;
    constexpr double phiHi = (kNDetPhi - 0.5) * kDetPhiStepDeg;

    constexpr double zLo = -1200.0;
    constexpr double zHi =  1200.0;

    auto* h = new TH2F(name, title,
                       kNDetPhi, phiLo, phiHi,
                       kNDetZ,   zLo,   zHi);

    h->GetXaxis()->SetTitle("detector #phi [deg]");
    h->GetYaxis()->SetTitle("detector z [mm]");
    h->GetZaxis()->SetTitle("ROI counts");
    h->SetStats(0);
    return h;
}

inline int DetectorPhiZBinX(const TH2F* h, int det)
{
    return h->GetXaxis()->FindBin(DetPhiDeg(det));
}

inline int DetectorPhiZBinY(const TH2F* h, int det)
{
    return h->GetYaxis()->FindBin(DetZMm(det));
}

inline void AddDetectorPhiZBin(TH2F* h, int det, double v)
{
    if (!h || !IsValidMultidetId(det)) return;

    const int bx = DetectorPhiZBinX(h, det);
    const int by = DetectorPhiZBinY(h, det);

    h->SetBinContent(bx, by, h->GetBinContent(bx, by) + v);
}

inline void SetDetectorPhiZBin(TH2F* h, int det, double v)
{
    if (!h || !IsValidMultidetId(det)) return;

    h->SetBinContent(DetectorPhiZBinX(h, det),
                     DetectorPhiZBinY(h, det),
                     v);
}

inline double GetDetectorPhiZBin(const TH2F* h, int det)
{
    if (!h || !IsValidMultidetId(det)) return 0.0;

    return h->GetBinContent(DetectorPhiZBinX(h, det),
                            DetectorPhiZBinY(h, det));
}

inline void DivideDetectorPhiZ(const TH2F* num,
                               const TH2F* den,
                               TH2F* out,
                               double zeroValue = 0.0)
{
    if (!num || !den || !out) return;

    for (int d = 0; d < kNMultidet; ++d) {
        const double n = GetDetectorPhiZBin(num, d);
        const double q = GetDetectorPhiZBin(den, d);
        SetDetectorPhiZBin(out, d, q > 0.0 ? n / q : zeroValue);
    }
}

} // namespace dcs

#endif

