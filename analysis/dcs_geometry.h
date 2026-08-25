#ifndef DCS_GEOMETRY_H
#define DCS_GEOMETRY_H
#include "geometry_constraints.h"
#include <TGraph.h>
#include <TH2Poly.h>
#include <TList.h>
#include <TVector2.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace dcs {
using namespace geo;

constexpr const char* kAnsiBold   = "\033[1m";
constexpr const char* kAnsiReset  = "\033[0m";
constexpr const char* kAnsiGreen  = "\033[32m";
constexpr const char* kAnsiYellow = "\033[33m";
constexpr const char* kAnsiRed    = "\033[31m";
constexpr const char* kAnsiCyan   = "\033[36m";

// Reproduces GeometryCASTOR440::GenerateFuelPositions().
inline std::vector<TVector2> BuildFuelPositions() {
    const double dx = kPitch, dy = kPitch * std::sqrt(3.0) / 2.0;
    const int R = kNRings;
    std::vector<TVector2> pts; pts.reserve(kNAssem);
    for (int q = -R; q <= R; ++q) {
        const int r1 = std::max(-R, -q - R), r2 = std::min(R, -q + R);
        for (int r = r1; r <= r2; ++r) {
            const double x = dx * (q + r / 2.0), y = dy * r;
            if (std::hypot(x, y) < 1e-6) continue;
            const double rmax = R * kPitch;
            if (std::abs(std::hypot(x, y) - rmax) < 1e-3) {
                double ang = std::atan2(y, x) * 180.0 / M_PI;
                if (ang < 0) ang += 360.0;
                bool corner = false;
                for (double a = 0; a < 360; a += 60)
                    if (std::abs(ang - a) < 1.0) { corner = true; break; }
                if (corner) continue;
            }
            pts.emplace_back(x, y);
        }
    }
    return pts;
}

inline TGraph* MakeHexPolygon(double cx, double cy, double scale = 1.0) {
    const double R = scale * kPitch / std::sqrt(3.0);
    auto* g = new TGraph(7);
    for (int i = 0; i < 6; ++i) {
        const double a = kPhiStart + i * 60.0 * M_PI / 180.;
        g->SetPoint(i, cx + R * std::cos(a), cy + R * std::sin(a));
    }
    g->SetPoint(6, g->GetX()[0], g->GetY()[0]);
    return g;
}

inline TH2Poly* MakeAssemblyHexMap(const char* name, const char* title,
                                   const std::vector<TVector2>& pts) {
    double xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9;
    for (const auto& p : pts) {
        xmin = std::min(xmin, p.X()); xmax = std::max(xmax, p.X());
        ymin = std::min(ymin, p.Y()); ymax = std::max(ymax, p.Y());
    }
    const double pad = 1.2 * kPitch;
    auto* h = new TH2Poly(name, title, xmin - pad, xmax + pad, ymin - pad, ymax + pad);
    for (const auto& p : pts) h->AddBin(MakeHexPolygon(p.X(), p.Y()));
    h->GetXaxis()->SetTitle("x [mm]");
    h->GetYaxis()->SetTitle("y [mm]");
    h->GetZaxis()->SetTitle("counts");
    return h;
}

inline void DrawMissingOverlay(const std::vector<TVector2>& pts,
                               const std::vector<bool>& missing) {
    for (size_t i = 0; i < pts.size(); ++i) {
        if (!missing[i]) continue;
        auto* g = MakeHexPolygon(pts[i].X(), pts[i].Y(), 0.97);
        g->SetFillColor(kWhite); g->SetFillStyle(1001);
        g->SetLineColor(kBlack); g->SetLineWidth(1);
        g->Draw("F SAME"); g->Draw("L SAME");
    }
}

// Rigid rotation of a whole TH2Poly about the cask centre (origin). Used to
// turn the diamond arrangement into a clean 2x2 grid. Bin contents preserved.
inline TH2Poly* RotateTH2Poly(TH2* src2, double deg, const char* newname) {
    auto* src = dynamic_cast<TH2Poly*>(src2);
    if (!src) return nullptr;
    const double th = deg * M_PI / 180., c = std::cos(th), s = std::sin(th);

    double xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30;
    TIter it0(src->GetBins()); TObject* o;
    while ((o = it0())) {
        auto* g = dynamic_cast<TGraph*>(((TH2PolyBin*)o)->GetPolygon());
        for (int i = 0; i < g->GetN(); ++i) {
            double X = c * g->GetX()[i] - s * g->GetY()[i];
            double Y = s * g->GetX()[i] + c * g->GetY()[i];
            xmin = std::min(xmin, X); xmax = std::max(xmax, X);
            ymin = std::min(ymin, Y); ymax = std::max(ymax, Y);
        }
    }
    auto* dst = new TH2Poly(newname, src->GetTitle(), xmin, xmax, ymin, ymax);
    dst->SetDirectory(nullptr);
    dst->GetXaxis()->SetTitle(src->GetXaxis()->GetTitle());
    dst->GetYaxis()->SetTitle(src->GetYaxis()->GetTitle());
    dst->GetZaxis()->SetTitle(src->GetZaxis()->GetTitle());

    TIter it(src->GetBins());
    while ((o = it())) {
        auto* b = (TH2PolyBin*)o;
        auto* g = dynamic_cast<TGraph*>(b->GetPolygon());
        auto* g2 = new TGraph(g->GetN());
        for (int i = 0; i < g->GetN(); ++i)
            g2->SetPoint(i, c * g->GetX()[i] - s * g->GetY()[i],
                            s * g->GetX()[i] + c * g->GetY()[i]);
        int nb = dst->AddBin(g2);
        dst->SetBinContent(nb, src->GetBinContent(b->GetBinNumber()));
    }
    dst->SetMinimum(src->GetMinimum());
    dst->SetMaximum(src->GetMaximum());
    return dst;
}

} // namespace dcs
#endif

