// optimize_gate.C
//   root -l 'optimize_gate.C("gamma","nominal")'
//   root -l 'optimize_gate.C("gamma","nominal",kCLYC,"SB",0.05,3,0.02,10, 1,1)'
//   root -l 'optimize_gate.C("gamma","nominal",kCLYC,"Q",0.05,3,0.02,10, 1,1, 50,200,200,0.35)'
//
// TWO-EDGE energy-gate optimizer. Scans EVERY window [lo,hi] and maximizes
//   FOM(lo,hi) = base(S,B) * (refPhi/RMSphi)^wPhi * (refZ/RMSz)^wZ
//
// base(S,B): "SB" = S/B (default, INTENSIVE), "Q" = 1 (quality-only: minimize
// the RMS product), or the extensive "SsqrtB"/"SsqrtSB" (kept for reference).
//   - S/B and RMS are INTENSIVE: unchanged if you measure 4x longer. With
//     effectively unlimited live-time the meaningful objective is data QUALITY
//     (contrast + localization), NOT raw counts, so extensive metrics like
//     S/sqrt(B) are the wrong objective and are not the default.
//   - MC-count floors (minS,minB) reject windows whose S/B and RMS are just
//     Poisson noise (e.g. the 11-signal/4-bkg razor bin). minSBR sets a hard
//     contrast floor so "Q" mode cannot pick a pure-but-empty sliver.
//
// RMSphi/RMSz are the FULL-RANGE std devs of cask-0 WALL hits with edep in
// [lo,hi], read from the energy-resolved maps h2_wall_{phi,z}_vs_edep. Built to
// equal ProjectionY(binLo,binHi)->GetRMS() exactly. These are NOT the FitFWHM
// core-fit widths in design_metrics.C (that fits only the peak core); RMS is
// tail-sensitive on purpose, which is what a gate optimizer needs.
//
// PHYSICS (verified from the data): the (phi,z) spot is NARROWEST in the Cs-137
// photopeak band (~0.5-0.74 MeV). A 662 keV photon cannot deposit >0.662 MeV,
// so everything above the Cs edge is penetrating Co-60/Eu-154 that reaches the
// detector from a BROAD range of wall positions -> wide. The <0.2 MeV tail is
// multiply-scattered/albedo -> also wide. Hence the localization lever is the
// HIGH edge (cut the Co/Eu tail), with a weak secondary gain from a small low
// edge. BOTH marginal sweeps are shown (canvases C and D).
#ifndef OPTIMIZE_GATE
#define OPTIMIZE_GATE 1

#include "dcs_setup.h"
#include "dcs_cuts.h"
#include "geometry_constraints.h"
#include "style.h"

#include <TAxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TStyle.h>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace PetroffPalette;
using namespace dcs;

namespace {

TString GateCacheName(const std::string& cfg, const std::string& part, int ic) {
    return Form("root_output/analysis_%s_%s_cask%d.root",
                cfg.c_str(), part.c_str(), ic);
}

enum FomKind { kSoverB, kSoverSqrtB, kSoverSqrtSpB, kQuality };
FomKind ParseFom(const std::string& s) {
    if (s == "SB"      || s == "S/B")          return kSoverB;
    if (s == "SsqrtSB" || s == "S/sqrt(S+B)")  return kSoverSqrtSpB;
    if (s == "Q"       || s == "quality")      return kQuality;
    return kSoverSqrtB;
}
const char* FomTitle(FomKind k) {
    switch (k) {
        case kSoverB:       return "S/B";
        case kSoverSqrtSpB: return "S/#sqrt{S+B}";
        case kQuality:      return "quality (narrowness only)";
        default:            return "S/#sqrt{B}";
    }
}
double BaseFom(FomKind k, double S, double B) {
    switch (k) {
        case kSoverB:       return (B > 0)     ? S / B              : 0.0;
        case kSoverSqrtSpB: return (S + B > 0) ? S / std::sqrt(S+B) : 0.0;
        case kQuality:      return 1.0;   // FOM = narrowness terms only
        default:            return (B > 0)     ? S / std::sqrt(B)   : 0.0;
    }
}

// Prefix sums of the 0th/1st/2nd y-moments per edep(x) bin of a TH2 (x=edep,
// y=spatial), summed over the FULL y-range so rms(lo,hi) == the value of
// ProjectionY(binLo,binHi)->GetRMS(). O(1) per window.
struct MomentScan {
    const TAxis* ax = nullptr;
    std::vector<double> cN, cS, cS2;
    bool ok = false;
    void build(TH2* h) {                         // full y-range, no band
        if (!h) { ok = false; return; }
        ax = h->GetXaxis();
        const int nx = ax->GetNbins(), ny = h->GetYaxis()->GetNbins();
        cN.assign(nx + 1, 0.0); cS.assign(nx + 1, 0.0); cS2.assign(nx + 1, 0.0);
        for (int ix = 1; ix <= nx; ++ix) {
            double n = 0, s = 0, s2 = 0;
            for (int iy = 1; iy <= ny; ++iy) {   // ALL y bins -> matches GetRMS
                const double c = h->GetBinContent(ix, iy);
                if (c == 0) continue;
                const double y = h->GetYaxis()->GetBinCenter(iy);
                n += c; s += c * y; s2 += c * y * y;
            }
            cN[ix]  = cN[ix - 1]  + n;
            cS[ix]  = cS[ix - 1]  + s;
            cS2[ix] = cS2[ix - 1] + s2;
        }
        ok = true;
    }
    double rms(double lo, double hi, double& nOut) const {
        nOut = 0;
        if (!ok) return 0;
        const int i = std::max(1, ax->FindBin(lo + 1e-9));
        const int j = std::min(ax->GetNbins(), ax->FindBin(hi - 1e-9));
        if (j < i) return 0;
        const double N  = cN[j]  - cN[i - 1];
        const double S1 = cS[j]  - cS[i - 1];
        const double S2 = cS2[j] - cS2[i - 1];
        nOut = N;
        if (N <= 0) return 0;
        const double var = S2 / N - (S1 / N) * (S1 / N);
        return var > 0 ? std::sqrt(var) : 0.0;
    }
};

} // anonymous

void optimize_gate(std::string  part    = "gamma",
                   std::string  cfgDir   = "nominal",
                   DetectorType det      = kCLYC,
                   std::string  fomName  = "SB",      // "SB"|"Q"|"SsqrtB"|"SsqrtSB"
                   double       eLo      = 0.0,       // scan low edge  [MeV]
                   double       eHi      = 2.0,       // scan high edge [MeV]
                   double       minWidth = 0.02,      // min gate width [MeV]
                   int          rebin    = 10,        // spectrum bin grouping
                   double       wPhi     = 1.0,       // phi-narrowness weight
                   double       wZ       = 1.0,       // z-narrowness weight
                   double       minSpatN = 50.0,      // min wall cts to trust RMS
                   double       minS     = 200.0,     // MC-count floor on signal
                   double       minB     = 200.0,     // MC-count floor on background
                   double       minSBR   = 0.0,       // reject windows below this S/B
                   double       tgtPhi   = 0.30,      // reference line [rad]
                   double       tgtZ     = 400.0)     // reference line [mm]
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(255);
    const FomKind fom   = ParseFom(fomName);
    const int     ncask = kNCask;

    // ---- build S (cask 0) & B (casks 1+2+3) + load cask-0 wall maps --------
    TH1D *S = nullptr, *B = nullptr;
    TH2  *phiMap = nullptr, *zMap = nullptr;
    for (int ic = 0; ic < ncask; ++ic) {
        TString fn = GateCacheName(cfgDir, part, ic);
        auto* f = TFile::Open(fn);
        if (!f || f->IsZombie()) {
            std::cerr << "[err] missing " << fn
                      << "\n       run plot_design_comp(\"" << part
                      << "\", ...) first to (re)generate caches.\n";
            if (f) f->Close();
            return;
        }
        auto* h = (TH1*)f->Get("h1_edep_sum_all_assemblies");
        if (!h) { std::cerr << "[err] no spectrum in " << fn << "\n"; f->Close(); return; }
        if      (ic == 0) { S = (TH1D*)h->Clone("S_spec"); S->SetDirectory(nullptr); }
        else if (!B)      { B = (TH1D*)h->Clone("B_spec"); B->SetDirectory(nullptr); }
        else              { B->Add(h); }

        if (ic == 0) {   // spatial maps only for the signal cask
            if (auto* p = (TH2*)f->Get("h2_wall_phi_vs_edep")) {
                phiMap = (TH2*)p->Clone("phiMap"); phiMap->SetDirectory(nullptr);
            }
            if (auto* z = (TH2*)f->Get("h2_wall_z_vs_edep")) {
                zMap = (TH2*)z->Clone("zMap"); zMap->SetDirectory(nullptr);
            }
        }
        f->Close();
    }
    if (!S || !B) { std::cerr << "[err] failed to build S/B spectra\n"; return; }

    const bool spatialOn = (phiMap && zMap && (wPhi > 0 || wZ > 0 || fom == kQuality));
    if (!spatialOn && (wPhi > 0 || wZ > 0 || fom == kQuality)) {
        std::cerr << "[warn] h2_wall_{phi,z}_vs_edep not found -- regenerate "
                     "caches. Falling back to pure-energy optimization.\n";
        if (fom == kQuality) {
            std::cerr << "[err] quality mode requires the wall maps. Aborting.\n";
            return;
        }
    }

    MomentScan Pm, Zm;
    if (spatialOn) { Pm.build(phiMap); Zm.build(zMap); }

    if (rebin > 1) { S->Rebin(rebin); B->Rebin(rebin); }
    const double wkeV = 1000.0 * S->GetBinWidth(1);

    // ---- current (hard-coded) gate + its reference widths ------------------
    CutConfig ref{ det, ParseParticle(part), kPidEnergy };
    double curLo, curHi; EnergyWindow(ref, curLo, curHi);
    auto integ = [](TH1D* h, double a, double b) {
        return h->Integral(h->GetXaxis()->FindBin(a + 1e-9),
                           h->GetXaxis()->FindBin(b - 1e-9));
    };
    const double curS = integ(S, curLo, curHi), curB = integ(B, curLo, curHi);
    double nTmp = 0;
    double refPhi = spatialOn ? Pm.rms(curLo, curHi, nTmp) : 1.0;
    double refZ   = spatialOn ? Zm.rms(curLo, curHi, nTmp) : 1.0;
    if (refPhi <= 0) refPhi = 1.0;
    if (refZ   <= 0) refZ   = 1.0;

    // ---- prefix sums for S/B (O(1) window integrals) -----------------------
    const int nb = S->GetNbinsX();
    std::vector<double> csS(nb + 1, 0.0), csB(nb + 1, 0.0);
    for (int b = 1; b <= nb; ++b) {
        csS[b] = csS[b - 1] + S->GetBinContent(b);
        csB[b] = csB[b - 1] + B->GetBinContent(b);
    }
    auto lowE = [&](int b) { return S->GetXaxis()->GetBinLowEdge(b); };
    auto upE  = [&](int b) { return S->GetXaxis()->GetBinUpEdge(b);  };
    const int bStart = std::max(1,  S->GetXaxis()->FindBin(eLo + 1e-9));
    const int bStop  = std::min(nb, S->GetXaxis()->FindBin(eHi - 1e-9));

    // combined FOM + physical outputs for one window. Returns 0 (rejected) for
    // windows below the MC-count floors or the contrast floor -> those cannot
    // win, killing both the "widest window" and "emptiest noise bin" failures.
    auto evalWindow = [&](double lo, double hi, double& S_, double& B_,
                          double& sPhi_, double& sZ_) {
        const int i = S->GetXaxis()->FindBin(lo + 1e-9);
        const int j = S->GetXaxis()->FindBin(hi - 1e-9);
        S_ = csS[j] - csS[i - 1];
        B_ = csB[j] - csB[i - 1];

        // widths always computed (used for the current-gate row + diagnostics)
        sPhi_ = sZ_ = 0;
        double Lphi = 1.0, Lz = 1.0;
        if (spatialOn) {
            double nP, nZ;
            sPhi_ = Pm.rms(lo, hi, nP);
            sZ_   = Zm.rms(lo, hi, nZ);
            if (sPhi_ > 0 && nP >= minSpatN) Lphi = std::pow(refPhi / sPhi_, wPhi);
            if (sZ_   > 0 && nZ >= minSpatN) Lz   = std::pow(refZ   / sZ_,   wZ);
        }

        // reliability + contrast vetoes
        if (S_ < minS || B_ < minB)          return 0.0;
        if (B_ > 0 && (S_ / B_) < minSBR)    return 0.0;

        return BaseFom(fom, S_, B_) * Lphi * Lz;
    };

    // ---- FULL TWO-EDGE brute-force scan over (lo,hi) + FOM map -------------
    const double axLo = lowE(bStart), axHi = upE(bStop);
    const int    nax  = bStop - bStart + 1;
    auto* hMap = new TH2D("h_fom_map",
        Form("gate FOM  (%s, w_{#phi}=%.2g, w_{z}=%.2g);"
             "gate LOW edge [MeV];gate HIGH edge [MeV]",
             FomTitle(fom), wPhi, wZ), nax, axLo, axHi, nax, axLo, axHi);

    double bestFom = 0.0, bestLo = 0, bestHi = 0, bS = 0, bB = 0, bP = 0, bZ = 0;
    for (int i = bStart; i <= bStop; ++i) {
        const double lo = lowE(i);
        for (int j = i; j <= bStop; ++j) {
            const double hi = upE(j);
            if (hi - lo < minWidth) continue;
            double s_, b_, sp_, sz_;
            const double fv = evalWindow(lo, hi, s_, b_, sp_, sz_);
            hMap->SetBinContent(i - bStart + 1, j - bStart + 1, fv);
            if (fv > bestFom) { bestFom = fv; bestLo = lo; bestHi = hi;
                                bS = s_; bB = b_; bP = sp_; bZ = sz_; }
        }
    }
    if (bestFom <= 0) {
        std::cerr << "[err] scan found no valid window -- loosen minS/minB/minSBR "
                     "(current: minS=" << minS << " minB=" << minB
                  << " minSBR=" << minSBR << ")\n";
        return;
    }

    double curSp = 0, curZ = 0, dS, dB;
    evalWindow(curLo, curHi, dS, dB, curSp, curZ);

    // ---- report ------------------------------------------------------------
    printf("\n=== TWO-EDGE gate optimisation (part=%s, cfg=%s, FOM=%s) ===\n",
           part.c_str(), cfgDir.c_str(), FomTitle(fom));
    printf("  scan [%.3f,%.3f] MeV | minWidth=%.3f | %.1f keV/bin | spatial=%s "
           "(w_phi=%.2g, w_z=%.2g)\n",
           eLo, eHi, minWidth, wkeV, spatialOn ? "ON" : "OFF", wPhi, wZ);
    printf("  floors: minS=%.0f minB=%.0f minSBR=%.3g minSpatN=%.0f\n",
           minS, minB, minSBR, minSpatN);
    printf("  %-8s %7s %7s %11s %11s %10s %11s %10s\n",
           "gate", "lo", "hi", "S(c0)", "B(1-11)", "S/B", "RMSphi[rad]", "RMSz[mm]");
    printf("  %-8s %7.3f %7.3f %11.4g %11.4g %10.3g %11.4g %10.4g\n",
           "current", curLo, curHi, curS, curB, curB>0?curS/curB:0, curSp, curZ);
    printf("  %-8s %7.3f %7.3f %11.4g %11.4g %10.3g %11.4g %10.4g\n",
           "OPTIMAL", bestLo, bestHi, bS, bB, bB>0?bS/bB:0, bP, bZ);
    printf("  base(%s): current=%.4g  optimal=%.4g\n",
           FomTitle(fom), BaseFom(fom, curS, curB), BaseFom(fom, bS, bB));
    if (spatialOn)
        printf("  RMS vs current:  RMSphi x%.3f   RMSz x%.3f\n",
               curSp>0?bP/curSp:0, curZ>0?bZ/curZ:0);
    printf("\n  paste into dcs_cuts.h EnergyWindow():\n");
    printf("     if (c.part == kGamma) { lo = %.3f; hi = %.3f; }\n\n",
           bestLo, bestHi);

    // ---- (A) S vs B spectra with the optimal band shaded -------------------
    auto* c1 = new TCanvas("c_gate_spectra",
        "signal vs background with optimal gate", 1100, 700);
    c1->SetLogy(); c1->SetGrid();
    const double ymax = 3.0 * std::max(1.0, std::max(S->GetMaximum(), B->GetMaximum()));
    S->SetTitle(Form("%s gate optimisation (cfg=%s);edep [MeV];counts",
                     part.c_str(), cfgDir.c_str()));
    S->SetLineColor(getP6(0)); S->SetLineWidth(2);
    B->SetLineColor(getP6(1)); B->SetLineWidth(2);
    S->GetXaxis()->SetRangeUser(eLo, eHi);
    S->SetMinimum(0.5); S->SetMaximum(ymax);
    S->Draw("HIST"); B->Draw("HIST SAME");
    auto* box = new TBox(bestLo, 0.5, bestHi, ymax);
    box->SetFillColorAlpha(kGreen + 2, 0.15);
    box->SetLineColor(kGreen + 2); box->SetLineStyle(2); box->Draw("l");
    S->Draw("HIST SAME"); B->Draw("HIST SAME");
    for (double x : { curLo, curHi }) {
        auto* l = new TLine(x, 0.5, x, ymax);
        l->SetLineColor(kGray + 2); l->SetLineStyle(3); l->Draw();
    }
    auto* leg = new TLegend(0.55, 0.68, 0.96, 0.90); leg->SetTextSize(0.026);
    leg->AddEntry(S, "signal (cask 0)", "l");
    leg->AddEntry(B, "background (casks 1-11)", "l");
    leg->AddEntry(box, Form("optimal [%.3f, %.3f] MeV", bestLo, bestHi), "f");
    leg->Draw();
    auto* tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.026);
    tx->DrawLatex(0.13, 0.16, Form("current [%.3f, %.3f] (grey dashed)",
                                   curLo, curHi));

    // ---- (B) THE two-way view: FOM map over (lo,hi) ------------------------
    auto* c2 = new TCanvas("c_gate_fom_map", "gate FOM map (two-edge)", 780, 700);
    c2->SetRightMargin(0.15);
    hMap->Draw("COLZ");
    auto* diag = new TLine(axLo, axLo, axHi, axHi);          // hi=lo boundary
    diag->SetLineStyle(2); diag->SetLineColor(kGray + 2); diag->Draw();
    auto* mk = new TMarker(bestLo, bestHi, 29);
    mk->SetMarkerColor(kRed + 1); mk->SetMarkerSize(2.4); mk->Draw();
    auto* mc = new TMarker(curLo, curHi, 24);
    mc->SetMarkerColor(kGray + 3); mc->SetMarkerSize(1.9); mc->Draw();
    auto* legM = new TLegend(0.16, 0.80, 0.55, 0.90); legM->SetTextSize(0.026);
    legM->AddEntry(mk, Form("optimum [%.3f,%.3f]", bestLo, bestHi), "p");
    legM->AddEntry(mc, Form("current [%.3f,%.3f]", curLo, curHi), "p");
    legM->Draw();

    // ---- marginal sweeps: hold one edge at the optimum, sweep the other ----
    //  mode=0 -> fix hi=bestHi, sweep lo ;  mode=1 -> fix lo=bestLo, sweep hi
    //  4 pads, each a REAL physical axis: S[counts], RMSphi[rad], RMSz[mm],
    //  FOM[arb.]. The high-edge sweep (D) is where localization moves most.
    if (spatialOn) {
        auto sweep = [&](int mode) {
            const bool low = (mode == 0);
            std::vector<double> vx, vS, vP, vZ, vF;
            if (low) {
                const int jHi = S->GetXaxis()->FindBin(bestHi - 1e-9);
                for (int i = bStart; i <= jHi; ++i) {
                    const double lo = lowE(i);
                    if (bestHi - lo < minWidth) break;
                    double s_, b_, sp_, sz_;
                    const double fv = evalWindow(lo, bestHi, s_, b_, sp_, sz_);
                    vx.push_back(lo); vS.push_back(s_);
                    vP.push_back(sp_); vZ.push_back(sz_); vF.push_back(fv);
                }
            } else {
                const int iLo = S->GetXaxis()->FindBin(bestLo + 1e-9);
                for (int j = iLo; j <= bStop; ++j) {
                    const double hi = upE(j);
                    if (hi - bestLo < minWidth) continue;
                    double s_, b_, sp_, sz_;
                    const double fv = evalWindow(bestLo, hi, s_, b_, sp_, sz_);
                    vx.push_back(hi); vS.push_back(s_);
                    vP.push_back(sp_); vZ.push_back(sz_); vF.push_back(fv);
                }
            }
            if (vx.size() < 2) return;

            const double fixedVal = low ? bestHi : bestLo;
            const char*  xlab     = low ? "gate LOW edge [MeV]"
                                        : "gate HIGH edge [MeV]";
            const char*  fixLbl   = low ? "hi" : "lo";
            const double markX    = low ? bestLo : bestHi;

            auto mkG = [&](std::vector<double>& y, Color_t col) {
                auto* g = new TGraph((int)vx.size());
                for (size_t k = 0; k < vx.size(); ++k) g->SetPoint((int)k, vx[k], y[k]);
                g->SetLineColor(col); g->SetLineWidth(3); g->SetMarkerColor(col);
                return g;
            };
            auto axMax = [&](std::vector<double>& v) {
                double m = 0; for (double x : v) m = std::max(m, x);
                return 1.15 * std::max(m, 1e-9);
            };
            auto vline = [&](double x, double y1) {
                auto* l = new TLine(x, 0, x, y1);
                l->SetLineColor(kRed + 1); l->SetLineStyle(2); l->SetLineWidth(2);
                l->Draw();
            };
            auto hline = [&](double y, const char* lbl) {
                auto* l = new TLine(vx.front(), y, vx.back(), y);
                l->SetLineColor(kGray + 2); l->SetLineStyle(3); l->SetLineWidth(2);
                l->Draw();
                auto* t = new TLatex(vx.front(), y, Form(" %s", lbl));
                t->SetTextSize(0.035); t->SetTextColor(kGray + 2); t->Draw();
            };

            auto* c = new TCanvas(
                Form("c_gate_sweep_%s", low ? "low" : "high"),
                low ? "sweep LOW edge (hi fixed)" : "sweep HIGH edge (lo fixed)",
                1200, 820);
            c->Divide(2, 2, 0.003, 0.003);

            // [1] signal S
            c->cd(1); gPad->SetGrid();
            auto* gS = mkG(vS, getP6(0));
            gS->SetTitle(Form("signal (%s = %.3f MeV);%s;S (cask 0) [counts]",
                              fixLbl, fixedVal, xlab));
            gS->GetYaxis()->SetRangeUser(0, axMax(vS));
            gS->Draw("AL"); vline(markX, axMax(vS));

            // [2] RMS_phi
            c->cd(2); gPad->SetGrid();
            auto* gP = mkG(vP, getP6(2));
            const double pMx = std::max(axMax(vP), 1.10 * tgtPhi);
            gP->SetTitle(Form("spot width #phi;%s;RMS_{#phi} [rad]", xlab));
            gP->GetYaxis()->SetRangeUser(0, pMx);
            gP->Draw("AL"); hline(tgtPhi, Form("ref %.2f rad", tgtPhi));
            vline(markX, pMx);

            // [3] RMS_z
            c->cd(3); gPad->SetGrid();
            auto* gZ = mkG(vZ, getP6(4));
            const double zMx = std::max(axMax(vZ), 1.10 * tgtZ);
            gZ->SetTitle(Form("spot width z;%s;RMS_{z} [mm]", xlab));
            gZ->GetYaxis()->SetRangeUser(0, zMx);
            gZ->Draw("AL"); hline(tgtZ, Form("ref %.0f mm", tgtZ));
            vline(markX, zMx);

            // [4] combined FOM
            c->cd(4); gPad->SetGrid();
            auto* gF = mkG(vF, kBlack);
            gF->SetTitle(Form("combined FOM (%s);%s;FOM [arb.]",
                              FomTitle(fom), xlab));
            gF->GetYaxis()->SetRangeUser(0, axMax(vF));
            gF->Draw("AL"); vline(markX, axMax(vF));
        };

        sweep(0);   // canvas C: low-edge marginal
        sweep(1);   // canvas D: high-edge marginal
    }
}
#endif







//  // optimize_gate.C
//  //   root -l 'optimize_gate.C("gamma","nominal")'
//  //   root -l 'optimize_gate.C("gamma","nominal",kCLYC,"SsqrtB",0,2,0.02,10, 1.0,1.0)'
//  //
//  // TWO-EDGE energy-gate optimizer. Scans EVERY window [lo,hi] and maximizes
//  //   FOM(lo,hi) = base(S,B) * (refPhi/sigmaPhi)^wPhi * (refZ/sigmaZ)^wZ
//  // base(S,B): S/sqrt(B) (default), S/sqrt(S+B), or S/B. sigmaPhi/sigmaZ are the
//  // FULL-RANGE RMS spreads of cask-0 WALL hits with edep in [lo,hi], read from the
//  // energy-resolved maps h2_wall_{phi,z}_vs_edep. refPhi/refZ are those RMS values
//  // at the CURRENT hard-coded gate, so the spatial factor == 1 there.
//  //
//  // PHYSICS (verified from the data): the (phi,z) spot is NARROWEST in the Cs-137
//  // photopeak band (~0.5-0.74 MeV). BOTH extremes broaden it: the <0.2 MeV tail is
//  // scattered/albedo (wide), and the >0.74 MeV region is penetrating Co-60/Eu-154
//  // (wide). Hence BOTH edges matter -> genuine two-edge optimization. The 2D FOM
//  // map (canvas B) is the two-way picture; canvases C and D are the low-edge and
//  // high-edge marginal sweeps through the optimum.
//  #ifndef OPTIMIZE_GATE
//  #define OPTIMIZE_GATE 1
//  
//  #include "dcs_setup.h"
//  #include "dcs_cuts.h"
//  #include "geometry_constraints.h"
//  #include "style.h"
//  
//  #include <TAxis.h>
//  #include <TBox.h>
//  #include <TCanvas.h>
//  #include <TFile.h>
//  #include <TGraph.h>
//  #include <TH1D.h>
//  #include <TH2.h>
//  #include <TH2D.h>
//  #include <TLatex.h>
//  #include <TLegend.h>
//  #include <TLine.h>
//  #include <TMarker.h>
//  #include <TStyle.h>
//  #include <cmath>
//  #include <cstdio>
//  #include <string>
//  #include <vector>
//  
//  using namespace PetroffPalette;
//  using namespace dcs;
//  
//  namespace {
//  
//  TString GateCacheName(const std::string& cfg, const std::string& part, int ic) {
//      return Form("root_output/analysis_%s_%s_cask%d.root",
//                  cfg.c_str(), part.c_str(), ic);
//  }
//  
//  enum FomKind { kSoverB, kSoverSqrtB, kSoverSqrtSpB };
//  FomKind ParseFom(const std::string& s) {
//      if (s == "SB"      || s == "S/B")         return kSoverB;
//      if (s == "SsqrtSB" || s == "S/sqrt(S+B)") return kSoverSqrtSpB;
//      return kSoverSqrtB;
//  }
//  const char* FomTitle(FomKind k) {
//      switch (k) {
//          case kSoverB:       return "S/B";
//          case kSoverSqrtSpB: return "S/#sqrt{S+B}";
//          default:            return "S/#sqrt{B}";
//      }
//  }
//  double BaseFom(FomKind k, double S, double B) {
//      switch (k) {
//          case kSoverB:       return (B > 0)     ? S / B              : 0.0;
//          case kSoverSqrtSpB: return (S + B > 0) ? S / std::sqrt(S+B) : 0.0;
//          default:            return (B > 0)     ? S / std::sqrt(B)   : 0.0;
//      }
//  }
//  
//  // Prefix sums of the 0th/1st/2nd y-moments per edep(x) bin of a TH2 (x=edep,
//  // y=spatial), summed over the FULL y-range so rms(lo,hi) == the value of
//  // ProjectionY(binLo,binHi)->GetRMS(). O(1) per window.
//  struct MomentScan {
//      const TAxis* ax = nullptr;
//      std::vector<double> cN, cS, cS2;
//      bool ok = false;
//      void build(TH2* h) {
//          if (!h) { ok = false; return; }
//          ax = h->GetXaxis();
//          const int nx = ax->GetNbins(), ny = h->GetYaxis()->GetNbins();
//          cN.assign(nx + 1, 0.0); cS.assign(nx + 1, 0.0); cS2.assign(nx + 1, 0.0);
//          for (int ix = 1; ix <= nx; ++ix) {
//              double n = 0, s = 0, s2 = 0;
//              for (int iy = 1; iy <= ny; ++iy) {         // ALL y bins
//                  const double c = h->GetBinContent(ix, iy);
//                  if (c == 0) continue;
//                  const double y = h->GetYaxis()->GetBinCenter(iy);
//                  n += c; s += c * y; s2 += c * y * y;
//              }
//              cN[ix]  = cN[ix - 1]  + n;
//              cS[ix]  = cS[ix - 1]  + s;
//              cS2[ix] = cS2[ix - 1] + s2;
//          }
//          ok = true;
//      }
//      double rms(double lo, double hi, double& nOut) const {
//          nOut = 0;
//          if (!ok) return 0;
//          const int i = std::max(1, ax->FindBin(lo + 1e-9));
//          const int j = std::min(ax->GetNbins(), ax->FindBin(hi - 1e-9));
//          if (j < i) return 0;
//          const double N  = cN[j]  - cN[i - 1];
//          const double S1 = cS[j]  - cS[i - 1];
//          const double S2 = cS2[j] - cS2[i - 1];
//          nOut = N;
//          if (N <= 0) return 0;
//          const double var = S2 / N - (S1 / N) * (S1 / N);
//          return var > 0 ? std::sqrt(var) : 0.0;
//      }
//  };
//  
//  } // anonymous
//  
//  void optimize_gate(std::string  part    = "gamma",
//                     std::string  cfgDir   = "nominal",
//                     DetectorType det      = kCLYC,
//                     std::string  fomName  = "SsqrtB", // "SB"|"SsqrtB"|"SsqrtSB"
//                     double       eLo      = 0.0,       // scan low edge  [MeV]
//                     double       eHi      = 2.0,       // scan high edge [MeV]
//                     double       minWidth = 0.02,      // min gate width [MeV]
//                     int          rebin    = 10,        // spectrum bin grouping
//                     double       wPhi     = 1.0,       // phi-narrowness weight
//                     double       wZ       = 1.0,       // z-narrowness weight
//                     double       minSpatN = 50.0)      // min wall cts to trust width
//  {
//      gStyle->SetOptStat(0);
//      gStyle->SetPalette(kBird);
//      gStyle->SetNumberContours(255);
//      const FomKind fom   = ParseFom(fomName);
//      const int     ncask = 4;
//  
//      // ---- build S (cask 0) & B (casks 1+2+3) + load cask-0 wall maps --------
//      TH1D *S = nullptr, *B = nullptr;
//      TH2  *phiMap = nullptr, *zMap = nullptr;
//      for (int ic = 0; ic < ncask; ++ic) {
//          TString fn = GateCacheName(cfgDir, part, ic);
//          auto* f = TFile::Open(fn);
//          if (!f || f->IsZombie()) {
//              std::cerr << "[err] missing " << fn
//                        << "\n       run plot_design_comp(\"" << part
//                        << "\", ...) first to (re)generate caches.\n";
//              if (f) f->Close();
//              return;
//          }
//          auto* h = (TH1*)f->Get("h1_edep_sum_all_assemblies");
//          if (!h) { std::cerr << "[err] no spectrum in " << fn << "\n"; f->Close(); return; }
//          if      (ic == 0) { S = (TH1D*)h->Clone("S_spec"); S->SetDirectory(nullptr); }
//          else if (!B)      { B = (TH1D*)h->Clone("B_spec"); B->SetDirectory(nullptr); }
//          else              { B->Add(h); }
//  
//          if (ic == 0) {   // spatial maps only for the signal cask
//              if (auto* p = (TH2*)f->Get("h2_wall_phi_vs_edep")) {
//                  phiMap = (TH2*)p->Clone("phiMap"); phiMap->SetDirectory(nullptr);
//              }
//              if (auto* z = (TH2*)f->Get("h2_wall_z_vs_edep")) {
//                  zMap = (TH2*)z->Clone("zMap"); zMap->SetDirectory(nullptr);
//              }
//          }
//          f->Close();
//      }
//      if (!S || !B) { std::cerr << "[err] failed to build S/B spectra\n"; return; }
//  
//      const bool spatialOn = (phiMap && zMap && (wPhi > 0 || wZ > 0));
//      if (!spatialOn && (wPhi > 0 || wZ > 0))
//          std::cerr << "[warn] h2_wall_{phi,z}_vs_edep not found -- regenerate "
//                       "caches. Falling back to pure-energy optimization.\n";
//  
//      MomentScan Pm, Zm;
//      if (spatialOn) { Pm.build(phiMap); Zm.build(zMap); }
//  
//      if (rebin > 1) { S->Rebin(rebin); B->Rebin(rebin); }
//      const double wkeV = 1000.0 * S->GetBinWidth(1);
//  
//      // ---- current (hard-coded) gate + its reference widths ------------------
//      CutConfig ref{ det, ParseParticle(part), kPidEnergy };
//      double curLo, curHi; EnergyWindow(ref, curLo, curHi);
//      auto integ = [](TH1D* h, double a, double b) {
//          return h->Integral(h->GetXaxis()->FindBin(a + 1e-9),
//                             h->GetXaxis()->FindBin(b - 1e-9));
//      };
//      const double curS = integ(S, curLo, curHi), curB = integ(B, curLo, curHi);
//      double nTmp = 0;
//      double refPhi = spatialOn ? Pm.rms(curLo, curHi, nTmp) : 1.0;
//      double refZ   = spatialOn ? Zm.rms(curLo, curHi, nTmp) : 1.0;
//      if (refPhi <= 0) refPhi = 1.0;
//      if (refZ   <= 0) refZ   = 1.0;
//  
//      // ---- prefix sums for S/B (O(1) window integrals) -----------------------
//      const int nb = S->GetNbinsX();
//      std::vector<double> csS(nb + 1, 0.0), csB(nb + 1, 0.0);
//      for (int b = 1; b <= nb; ++b) {
//          csS[b] = csS[b - 1] + S->GetBinContent(b);
//          csB[b] = csB[b - 1] + B->GetBinContent(b);
//      }
//      auto lowE = [&](int b) { return S->GetXaxis()->GetBinLowEdge(b); };
//      auto upE  = [&](int b) { return S->GetXaxis()->GetBinUpEdge(b);  };
//      const int bStart = std::max(1,  S->GetXaxis()->FindBin(eLo + 1e-9));
//      const int bStop  = std::min(nb, S->GetXaxis()->FindBin(eHi - 1e-9));
//  
//      // combined FOM + physical outputs for one window
//      auto evalWindow = [&](double lo, double hi, double& S_, double& B_,
//                            double& sPhi_, double& sZ_) {
//          const int i = S->GetXaxis()->FindBin(lo + 1e-9);
//          const int j = S->GetXaxis()->FindBin(hi - 1e-9);
//          S_ = csS[j] - csS[i - 1];
//          B_ = csB[j] - csB[i - 1];
//          const double base = BaseFom(fom, S_, B_);
//          sPhi_ = sZ_ = 0;
//          double Lphi = 1.0, Lz = 1.0;
//          if (spatialOn) {
//              double nP, nZ;
//              sPhi_ = Pm.rms(lo, hi, nP);
//              sZ_   = Zm.rms(lo, hi, nZ);
//              if (sPhi_ > 0 && nP >= minSpatN) Lphi = std::pow(refPhi / sPhi_, wPhi);
//              if (sZ_   > 0 && nZ >= minSpatN) Lz   = std::pow(refZ   / sZ_,   wZ);
//          }
//          return base * Lphi * Lz;
//      };
//  
//      // ---- FULL TWO-EDGE brute-force scan over (lo,hi) + FOM map -------------
//      const double axLo = lowE(bStart), axHi = upE(bStop);
//      const int    nax  = bStop - bStart + 1;
//      auto* hMap = new TH2D("h_fom_map",
//          Form("combined gate FOM  (%s, w_{#phi}=%.2g, w_{z}=%.2g);"
//               "gate LOW edge [MeV];gate HIGH edge [MeV]",
//               FomTitle(fom), wPhi, wZ), nax, axLo, axHi, nax, axLo, axHi);
//  
//      double bestFom = -1, bestLo = 0, bestHi = 0, bS = 0, bB = 0, bP = 0, bZ = 0;
//      for (int i = bStart; i <= bStop; ++i) {
//          const double lo = lowE(i);
//          for (int j = i; j <= bStop; ++j) {
//              const double hi = upE(j);
//              if (hi - lo < minWidth) continue;
//              double s_, b_, sp_, sz_;
//              const double fv = evalWindow(lo, hi, s_, b_, sp_, sz_);
//              hMap->SetBinContent(i - bStart + 1, j - bStart + 1, fv);
//              if (fv > bestFom) { bestFom = fv; bestLo = lo; bestHi = hi;
//                                  bS = s_; bB = b_; bP = sp_; bZ = sz_; }
//          }
//      }
//      if (bestFom < 0) { std::cerr << "[err] scan found no valid window\n"; return; }
//  
//      double curSp = 0, curZ = 0, dS, dB;
//      evalWindow(curLo, curHi, dS, dB, curSp, curZ);
//  
//      // ---- report ------------------------------------------------------------
//      printf("\n=== TWO-EDGE gate optimisation (part=%s, cfg=%s, FOM=%s) ===\n",
//             part.c_str(), cfgDir.c_str(), FomTitle(fom));
//      printf("  scan [%.3f,%.3f] MeV | minWidth=%.3f | %.1f keV/bin | "
//             "spatial=%s (w_phi=%.2g, w_z=%.2g)\n",
//             eLo, eHi, minWidth, wkeV, spatialOn ? "ON" : "OFF", wPhi, wZ);
//      printf("  %-8s %7s %7s %11s %11s %10s %11s %10s\n",
//             "gate", "lo", "hi", "S(c0)", "B(123)", "S/B", "RMSphi[rad]", "RMSz[mm]");
//      printf("  %-8s %7.3f %7.3f %11.4g %11.4g %10.3g %11.4g %10.4g\n",
//             "current", curLo, curHi, curS, curB, curB>0?curS/curB:0, curSp, curZ);
//      printf("  %-8s %7.3f %7.3f %11.4g %11.4g %10.3g %11.4g %10.4g\n",
//             "OPTIMAL", bestLo, bestHi, bS, bB, bB>0?bS/bB:0, bP, bZ);
//      printf("  base(%s): current=%.4g  optimal=%.4g\n",
//             FomTitle(fom), BaseFom(fom, curS, curB), BaseFom(fom, bS, bB));
//      if (spatialOn)
//          printf("  RMS vs current:  RMSphi x%.3f   RMSz x%.3f\n",
//                 curSp>0?bP/curSp:0, curZ>0?bZ/curZ:0);
//      printf("\n  paste into dcs_cuts.h EnergyWindow():\n");
//      printf("     if (c.part == kGamma) { lo = %.3f; hi = %.3f; }\n\n",
//             bestLo, bestHi);
//  
//      // ---- (A) S vs B spectra with the optimal band shaded -------------------
//      auto* c1 = new TCanvas("c_gate_spectra",
//          "signal vs background with optimal gate", 1100, 700);
//      c1->SetLogy(); c1->SetGrid();
//      const double ymax = 3.0 * std::max(1.0, std::max(S->GetMaximum(), B->GetMaximum()));
//      S->SetTitle(Form("%s gate optimisation (cfg=%s);edep [MeV];counts",
//                       part.c_str(), cfgDir.c_str()));
//      S->SetLineColor(getP6(0)); S->SetLineWidth(2);
//      B->SetLineColor(getP6(1)); B->SetLineWidth(2);
//      S->GetXaxis()->SetRangeUser(eLo, eHi);
//      S->SetMinimum(0.5); S->SetMaximum(ymax);
//      S->Draw("HIST"); B->Draw("HIST SAME");
//      auto* box = new TBox(bestLo, 0.5, bestHi, ymax);
//      box->SetFillColorAlpha(kGreen + 2, 0.15);
//      box->SetLineColor(kGreen + 2); box->SetLineStyle(2); box->Draw("l");
//      S->Draw("HIST SAME"); B->Draw("HIST SAME");
//      for (double x : { curLo, curHi }) {
//          auto* l = new TLine(x, 0.5, x, ymax);
//          l->SetLineColor(kGray + 2); l->SetLineStyle(3); l->Draw();
//      }
//      auto* leg = new TLegend(0.55, 0.68, 0.96, 0.90); leg->SetTextSize(0.026);
//      leg->AddEntry(S, "signal (cask 0)", "l");
//      leg->AddEntry(B, "background (casks 1+2+3)", "l");
//      leg->AddEntry(box, Form("optimal [%.3f, %.3f] MeV", bestLo, bestHi), "f");
//      leg->Draw();
//      auto* tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.026);
//      tx->DrawLatex(0.13, 0.16, Form("current [%.3f, %.3f] (grey dashed)",
//                                     curLo, curHi));
//  
//      // ---- (B) THE two-way view: combined-FOM map over (lo,hi) ---------------
//      auto* c2 = new TCanvas("c_gate_fom_map", "gate FOM map (two-edge)", 780, 700);
//      c2->SetRightMargin(0.15);
//      hMap->Draw("COLZ");
//      auto* diag = new TLine(axLo, axLo, axHi, axHi);          // hi=lo boundary
//      diag->SetLineStyle(2); diag->SetLineColor(kGray + 2); diag->Draw();
//      auto* mk = new TMarker(bestLo, bestHi, 29);
//      mk->SetMarkerColor(kRed + 1); mk->SetMarkerSize(2.4); mk->Draw();
//      auto* mc = new TMarker(curLo, curHi, 24);
//      mc->SetMarkerColor(kGray + 3); mc->SetMarkerSize(1.9); mc->Draw();
//      auto* legM = new TLegend(0.16, 0.80, 0.55, 0.90); legM->SetTextSize(0.026);
//      legM->AddEntry(mk, Form("optimum [%.3f,%.3f]", bestLo, bestHi), "p");
//      legM->AddEntry(mc, Form("current [%.3f,%.3f]", curLo, curHi), "p");
//      legM->Draw();
//  
//      // ---- marginal-sweep helper: hold one edge, sweep the other -------------
//      //  mode=0 -> fix hi=bestHi, sweep lo ;  mode=1 -> fix lo=bestLo, sweep hi
//      if (spatialOn) {
//          auto sweep = [&](int mode) {
//              std::vector<double> vx, vS, vP, vZ, vF;
//              if (mode == 0) {                       // sweep LOW edge
//                  const int jHi = S->GetXaxis()->FindBin(bestHi - 1e-9);
//                  for (int i = bStart; i <= jHi; ++i) {
//                      const double lo = lowE(i);
//                      if (bestHi - lo < minWidth) break;
//                      double s_, b_, sp_, sz_;
//                      const double fv = evalWindow(lo, bestHi, s_, b_, sp_, sz_);
//                      vx.push_back(lo); vS.push_back(s_);
//                      vP.push_back(sp_); vZ.push_back(sz_); vF.push_back(fv);
//                  }
//              } else {                               // sweep HIGH edge
//                  const int iLo = S->GetXaxis()->FindBin(bestLo + 1e-9);
//                  for (int j = iLo; j <= bStop; ++j) {
//                      const double hi = upE(j);
//                      if (hi - bestLo < minWidth) continue;
//                      double s_, b_, sp_, sz_;
//                      const double fv = evalWindow(bestLo, hi, s_, b_, sp_, sz_);
//                      vx.push_back(hi); vS.push_back(s_);
//                      vP.push_back(sp_); vZ.push_back(sz_); vF.push_back(fv);
//                  }
//              }
//              const bool low = (mode == 0);
//              const double fixedVal = low ? bestHi : bestLo;
//              const char* xlab = low ? "gate LOW edge [MeV]" : "gate HIGH edge [MeV]";
//              const double mark = low ? bestLo : bestHi;
//  
//              auto mkG = [&](std::vector<double>& y, Color_t col) {
//                  auto* g = new TGraph((int)vx.size());
//                  for (size_t k = 0; k < vx.size(); ++k) g->SetPoint((int)k, vx[k], y[k]);
//                  g->SetLineColor(col); g->SetLineWidth(3); g->SetMarkerColor(col);
//                  return g;
//              };
//              auto axMax = [&](std::vector<double>& v) {
//                  double m = 0; for (double x : v) m = std::max(m, x); return 1.15*std::max(m,1e-9);
//              };
//              auto vline = [&](double x, double y1) {
//                  auto* l = new TLine(x, 0, x, y1);
//                  l->SetLineColor(kRed + 1); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
//              };
//  
//              auto* c = new TCanvas(
//                  Form("c_gate_sweep_%s", low ? "low" : "high"),
//                  low ? "sweep low edge (hi fixed)" : "sweep high edge (lo fixed)",
//                  1200, 820);
//              c->Divide(2, 2, 0.003, 0.003);
//              if (vx.empty()) return;
//  
//              c->cd(1); gPad->SetGrid();
//              auto* gS = mkG(vS, getP6(0));
//              gS->SetTitle(Form("signal (%s = %.3f MeV);%s;S (cask 0) [counts]",
//                                low ? "hi" : "lo", fixedVal, xlab));
//              gS->GetYaxis()->SetRangeUser(0, axMax(vS));
//              gS->Draw("AL"); vline(mark, axMax(vS));
//  
//              c->cd(2); gPad->SetGrid();
//              auto* gP = mkG(vP, getP6(2));
//              gP->SetTitle(Form("spot width #phi;%s;RMS_{#phi} [rad]", xlab));
//              gP->GetYaxis()->SetRangeUser(0, axMax(vP));
//              gP->Draw("AL"); vline(mark, axMax(vP));
//  
//              c->cd(3); gPad->SetGrid();
//              auto* gZ = mkG(vZ, getP6(4));
//              gZ->SetTitle(Form("spot width z;%s;RMS_{z} [mm]", xlab));
//              gZ->GetYaxis()->SetRangeUser(0, axMax(vZ));
//              gZ->Draw("AL"); vline(mark, axMax(vZ));
//  
//              c->cd(4); gPad->SetGrid();
//              auto* gF = mkG(vF, kBlack);
//              gF->SetTitle(Form("combined FOM (%s);%s;FOM [arb.]", FomTitle(fom), xlab));
//              gF->GetYaxis()->SetRangeUser(0, axMax(vF));
//              gF->Draw("AL"); vline(mark, axMax(vF));
//          };
//  
//          sweep(0);   // canvas C: low-edge marginal
//          sweep(1);   // canvas D: high-edge marginal
//      }
//  }
//  #endif
//  
