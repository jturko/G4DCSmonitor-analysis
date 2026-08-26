// design_metrics.C
//   root -l 'design_metrics.C("gamma", kCLYC, 2, 0.30, 400.)'
//
// Cross-design metrics from the per-cask files in root_output/ produced by
// analyze_cask.C. Assumes equal activity per assembly (raw ratios valid).
//
// Config labels come from dcs_setup.h: configs[j].dir is the raw directory tag
// (used to locate files); configs[j].title is the short human label (used on
// every axis / legend). Nominal (index 0) is drawn solid black on the overlays.
#ifndef DESIGN_METRICS
#define DESIGN_METRICS 1

#include "dcs_setup.h"
#include "dcs_cuts.h"
#include "utils.h"
#include "style.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace PetroffPalette;
using namespace dcs;

namespace {

// Integral of a TH2Poly's bins (the per-cask total detected signal).
double PolyIntegral(TH2Poly* h, double* err = nullptr) {
    double s = 0, e2 = 0;
    for (int b = 1; b <= h->GetNumberOfBins(); ++b) {
        s += h->GetBinContent(b);
        double be = h->GetBinError(b); e2 += be * be;
    }
    if (err) *err = std::sqrt(e2);
    return s;
}

TString CacheName(const std::string& cfg, const std::string& part, int ic) {
    return Form("root_output/analysis_%s_%s_cask%d.root", cfg.c_str(),
                part.c_str(), ic);
}

// FWHM from a Gaussian fit over the histogram's core; falls back to the
// RMS-based estimate if the fit fails.
double FitFWHM(TH1* h, double& fwhmErr) {
    fwhmErr = 0;
    if (!h || h->GetEntries() < 20) return 0;
    double pk = h->GetBinCenter(h->GetMaximumBin());
    double rms = h->GetRMS();
    TF1 g("g", "gaus", pk - 1*rms, pk + 1*rms);
    g.SetParameters(h->GetMaximum(), pk, rms);
    g.SetNpx(400);
    //new TCanvas;
    if (h->Fit(&g, "QNR", "", pk-1*rms, pk+1*rms) == 0) {
        //h->DrawClone();
        //g.DrawClone("same");
        fwhmErr = 2.3548 * g.GetParError(2);
        return 2.3548 * g.GetParameter(2);
    }
    return 2.3548 * rms;
}
//double FitFWHM(TH1* h, double& fwhmErr) {
//    fwhmErr = abs(2.3548 * h->GetRMSError());
//    cout << "StdDev: "  << h->GetRMS() << ", Error: " << h->GetRMSError() << endl;
//    return abs(2.3548 *    h->GetRMS());
//}

// Per-config overlay line style: index 0 (nominal) -> solid black, heavier;
// others -> Petroff palette shifted so config 1 == getP6(0), reused dashed on
// wrap. Line width 1 (2 looked too heavy on dense overlays).
void StyleConfigLine(TH1* h, int j) {
    if (!h) return;
    if (j == 0) { h->SetLineColor(kBlack);              h->SetLineStyle(1);
                  h->SetLineWidth(2); }
    else        { h->SetLineColor(getP6((j - 1) % 6));  h->SetLineStyle(1 + (j - 1) / 6);
                  h->SetLineWidth(1); }
}

} // anonymous

void design_metrics(std::string part = "gamma",
                    DetectorType det  = kCLYC,
                    int level         = 2,
                    double phiGate    = 0.30,   // |phi| < phiGate  [rad]
                    double zGate      = 400.0)  // |z|   < zGate    [mm]
{
    //SetDivergingPalette();
    gStyle->SetOptStat(0);
    const auto& configs = Configs();
    const int nconfig = (int)configs.size(), ncask = kNCask;

    // convenience accessors for the two label flavours
    auto D = [&](int j) { return configs[j].dir.c_str();   };  // raw tag (files)
    auto T = [&](int j) { return configs[j].title.c_str(); };  // pretty (plots)

    // per-design accumulators
    std::vector<double> sig(nconfig), sigE(nconfig);       // cask 0 total
    std::vector<double> bkg(nconfig), bkgE(nconfig);       // casks 1-11
    std::vector<double> sbr(nconfig), impr(nconfig);       // S/B, x nominal
    std::vector<double> fwPhi(nconfig), fwPhiE(nconfig);
    std::vector<double> fwZ(nconfig),   fwZE(nconfig);
    std::vector<double> central(nconfig);                  // gated purity
    std::vector<double> meanCos(nconfig);
    std::vector<TH2F*>  emit(nconfig, nullptr);            // summed money plot
    std::vector<TH1*> sigSpec(nconfig, nullptr);   // cask-0 spectrum
    std::vector<TH1*> bkgSpec(nconfig, nullptr);   // casks 1-11 summed spectrum
    std::vector<TH1*> phiProj(nconfig, nullptr);   // cask-0 phi projection
    std::vector<TH1*> zProj(nconfig, nullptr);     // cask-0 z   projection

    for (int j = 0; j < nconfig; ++j) {
        double s = 0, sE2 = 0, b = 0, bE2 = 0;
        TH1* zphi_px = nullptr; TH1* zphi_py = nullptr;   // cask-0 projections
        double cosNum = 0, cosDen = 0;

        for (int ic = 0; ic < ncask; ++ic) {
            auto* f = TFile::Open(CacheName(configs[j].dir, part, ic));
            if (!f || f->IsZombie()) {
                std::cerr << "[warn] missing " << CacheName(configs[j].dir, part, ic)
                          << " -- run plot_design_comp first\n";
                if (f) f->Close();
                continue;
            }
            // total detected signal per cask (hex map)
            if (auto* hp = (TH2Poly*)f->Get("h2_counts_total")) {
                double e; double v = PolyIntegral(hp, &e);
                if (ic == 0) { s += v; sE2 += e*e; }
                else         { b += v; bE2 += e*e; }
            }
            // directionality (all casks contribute)
            if (auto* hc = (TH1*)f->Get("h1_cos_alpha")) {
                cosNum += hc->GetMean() * hc->GetEntries();
                cosDen += hc->GetEntries();
            }
            // cask-0 wall map -> projections + central-gate purity
            if (ic == 0) {
                if (auto* hw = (TH2F*)f->Get("h2_primary_zphi_wall")) {
                    hw->SetDirectory(nullptr);
                    const int xlo = hw->GetXaxis()->FindBin(-phiGate);
                    const int xhi = hw->GetXaxis()->FindBin(+phiGate);
                    const int ylo = hw->GetYaxis()->FindBin(-zGate);
                    const int yhi = hw->GetYaxis()->FindBin(+zGate);
                    double tot = hw->Integral();
                    double gate = hw->Integral(xlo, xhi, ylo, yhi);
                    central[j] = (tot > 0) ? gate / tot : 0;
                    zphi_px = hw->ProjectionX(Form("px_%d", j));
                    zphi_py = hw->ProjectionY(Form("py_%d", j));
                    zphi_px->SetDirectory(nullptr);
                    zphi_py->SetDirectory(nullptr);
                }
            }

            // full-range spectrum: cask 0 -> signal, casks 1-11 -> background
            if (auto* hs = (TH1*)f->Get("h1_edep_sum_all_assemblies")) {
                if (ic == 0) {
                    sigSpec[j] = (TH1*)hs->Clone(Form("sigspec_%d", j));
                    sigSpec[j]->SetDirectory(nullptr);
                } else if (!bkgSpec[j]) {
                    bkgSpec[j] = (TH1*)hs->Clone(Form("bkgspec_%d", j));
                    bkgSpec[j]->SetDirectory(nullptr);
                } else {
                    bkgSpec[j]->Add(hs);
                }
            }

            f->Close();
        }

        sig[j]  = s; sigE[j] = std::sqrt(sE2);
        bkg[j]  = b; bkgE[j] = std::sqrt(bE2);
        sbr[j]  = (b > 0) ? s / b : 0;
        meanCos[j] = (cosDen > 0) ? cosNum / cosDen : 0;
        if (zphi_px) fwPhi[j] = FitFWHM(zphi_px, fwPhiE[j]);
        if (zphi_py) fwZ[j]   = FitFWHM(zphi_py, fwZE[j]);
        phiProj[j] = zphi_px;   // already SetDirectory(nullptr) above -> survives
        zProj[j]   = zphi_py;
    }
    // improvement factor relative to nominal (index 0)
    for (int j = 0; j < nconfig; ++j)
        impr[j] = (sbr[0] > 0) ? sbr[j] / sbr[0] : 0;

    // ---- (1) scorecard table -----------------------------------------------
    printf("\n=== design scorecard  (part=%s, gate |phi|<%.2f rad, |z|<%.0f mm) ===\n",
           part.c_str(), phiGate, zGate);
    printf("%-22s %10s %10s %8s %8s %10s %8s %8s\n",
           "config", "sig(c0)", "bkg(1-11)", "S/B", "xNom",
           "central", "FWHMphi", "FWHMz");
    for (int j = 0; j < nconfig; ++j)
        printf("%-22s %10.4g %10.4g %8.3f %8.2f %10.3f %8.3f %8.1f\n",
               T(j), sig[j], bkg[j], sbr[j], impr[j],
               central[j], fwPhi[j], fwZ[j]);
    printf("\n");

    // ---- (2) S/B improvement bar chart -------------------------------------
    auto* cImp = new TCanvas("c_improvement", "S/B improvement x nominal", 900, 600);
    cImp->SetBottomMargin(0.28);
    auto* hImp = new TH1F("h_improvement",
        "signal-to-background improvement;;S/B  (#times nominal)", nconfig, 0, nconfig);
    for (int j = 0; j < nconfig; ++j) {
        hImp->SetBinContent(j + 1, impr[j]);
        hImp->GetXaxis()->SetBinLabel(j + 1, T(j));
    }
    hImp->SetFillColor(getP6(0)); hImp->GetXaxis()->LabelsOption("v");
    hImp->GetXaxis()->SetLabelSize(0.03); hImp->SetMinimum(0);
    hImp->Draw("hist");
    auto* one = new TLine(0, 1, nconfig, 1);
    one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw();

    // ---- (4) collimation FWHM per design (phi and z together) --------------
    // A FRAME HISTOGRAM carries the categorical x-axis (native SetBinLabel +
    // LabelsOption("v"), same as canvases 2 & 6) and the left FWHM_phi axis.
    // Native axis labels are redrawn on every zoom -- the old manual TLatex at
    // y=0 vanished the moment the frame rescaled, which caused the labels to
    // disappear on zoom. Short titles from dcs_setup now fit cleanly vertical.
    auto* cFW = new TCanvas("c_fwhm", "signal-spot width", 950, 650);
    cFW->SetBottomMargin(0.28); cFW->SetRightMargin(0.12);

    auto* gPhi = new TGraphErrors(nconfig);
    auto* gZ   = new TGraphErrors(nconfig);
    double phiMax = 0, zMax = 0;
    for (int j = 0; j < nconfig; ++j) {
        gPhi->SetPoint(j, j + 0.5, fwPhi[j]); gPhi->SetPointError(j, 0, fwPhiE[j]);
        gZ  ->SetPoint(j, j + 0.5, fwZ[j]);   gZ  ->SetPointError(j, 0, fwZE[j]);
        phiMax = std::max(phiMax, fwPhi[j] + fwPhiE[j]);
        zMax   = std::max(zMax,   fwZ[j]   + fwZE[j]);
    }
    if (phiMax <= 0) phiMax = 1;  if (zMax <= 0) zMax = 1;

    // frame: owns the x labels + left (phi) axis
    auto* hFrame = new TH1F("h_fwhm_frame",
        "cask-0 signal-spot FWHM;;FWHM_{#phi} [rad]", nconfig, 0, nconfig);
    for (int j = 0; j < nconfig; ++j)
        hFrame->GetXaxis()->SetBinLabel(j + 1, T(j));
    hFrame->SetStats(0);
    hFrame->SetMinimum(0); hFrame->SetMaximum(1.15 * phiMax);
    hFrame->GetXaxis()->LabelsOption("v"); hFrame->GetXaxis()->SetLabelSize(0.03);
    hFrame->GetYaxis()->SetTitleColor(getP6(2));
    hFrame->GetYaxis()->SetLabelColor(getP6(2));
    hFrame->GetYaxis()->SetAxisColor(getP6(2));
    hFrame->SetLineColor(kWhite);        // hide the empty flat line at y=0
    hFrame->Draw("HIST");

    // left series: FWHM_phi [rad]
    gPhi->SetMarkerStyle(20); gPhi->SetMarkerColor(getP6(2));
    gPhi->SetLineColor(getP6(2)); gPhi->SetLineWidth(2);
    gPhi->Draw("P SAME");

    // right series: FWHM_z [mm], rescaled onto the phi frame
    const double sc = (1.15 * phiMax) / (1.15 * zMax);
    auto* gZscaled = new TGraphErrors(nconfig);
    for (int j = 0; j < nconfig; ++j) {
        gZscaled->SetPoint(j, j + 0.5, fwZ[j] * sc);
        gZscaled->SetPointError(j, 0, fwZE[j] * sc);
    }
    gZscaled->SetMarkerStyle(21); gZscaled->SetMarkerColor(getP6(4));
    gZscaled->SetLineColor(getP6(4)); gZscaled->SetLineWidth(2);
    gZscaled->Draw("P SAME");

    auto* axZ = new TGaxis(nconfig, 0, nconfig, 1.15 * phiMax,
                           0, 1.15 * zMax, 510, "+L");
    axZ->SetTitle("FWHM_{z} [mm]");
    axZ->SetLineColor(getP6(4)); axZ->SetLabelColor(getP6(4));
    axZ->SetTitleColor(getP6(4)); axZ->Draw();

    auto* legFW = new TLegend(0.15, 0.80, 0.45, 0.90);
    legFW->AddEntry(gPhi,     "FWHM_{#phi} (left)",  "lp");
    legFW->AddEntry(gZscaled, "FWHM_{z} (right)",    "lp");
    legFW->Draw();

    // ---- (6) absolute signal vs background, grouped bars (log-y) -----------
    //  Exposes the mechanism the ratio chart (c_improvement) hides: does a
    //  design suppress background WITHOUT also killing the cask-0 signal?
    auto* cSvB = new TCanvas("c_sig_vs_bkg", "signal vs background", 950, 680);
    cSvB->SetBottomMargin(0.28); cSvB->SetLogy();
    auto* hSig = new TH1F("h_sig_c0",   "", nconfig, 0, nconfig);
    auto* hBkg = new TH1F("h_bkg_c123", "", nconfig, 0, nconfig);
    double ymaxSB = 1.0;
    for (int j = 0; j < nconfig; ++j) {
        hSig->SetBinContent(j + 1, sig[j]); hSig->SetBinError(j + 1, sigE[j]);
        hBkg->SetBinContent(j + 1, bkg[j]); hBkg->SetBinError(j + 1, bkgE[j]);
        hSig->GetXaxis()->SetBinLabel(j + 1, T(j));
        ymaxSB = std::max({ymaxSB, sig[j], bkg[j]});
    }
    hSig->SetFillColor(getP6(0)); hSig->SetBarWidth(0.40); hSig->SetBarOffset(0.08);
    hBkg->SetFillColor(getP6(1)); hBkg->SetBarWidth(0.40); hBkg->SetBarOffset(0.52);
    hSig->SetMinimum(0.5); hSig->SetMaximum(3.0 * ymaxSB);
    hSig->GetYaxis()->SetTitle("detected counts");
    hSig->GetXaxis()->LabelsOption("v"); hSig->GetXaxis()->SetLabelSize(0.03);
    hSig->Draw("bar"); hBkg->Draw("bar same");
    auto* legSvB = new TLegend(0.62, 0.80, 0.88, 0.90);
    legSvB->AddEntry(hSig, "signal (cask 0)", "f");
    legSvB->AddEntry(hBkg, "background (casks 1-11)", "f");
    legSvB->Draw();

    // ---- (7) retention vs suppression, ALL designs (ratio to nominal) ------
    //  Generalizes plot (3) beyond the Pb boron scan. Ideal design: retention
    //  (orange) pinned near 1, suppression (blue) pushed toward 0.
    auto* gRetAll = new TGraphErrors(nconfig);   // sig[j]/sig[0]
    auto* gSupAll = new TGraphErrors(nconfig);   // bkg[j]/bkg[0]
    for (int j = 0; j < nconfig; ++j) {
        double ret = (sig[0] > 0) ? sig[j] / sig[0] : 0.0;
        double sup = (bkg[0] > 0) ? bkg[j] / bkg[0] : 0.0;
        double retE = (sig[0] > 0 && sig[j] > 0)
            ? ret * std::sqrt((sigE[j]/sig[j])*(sigE[j]/sig[j]) +
                              (sigE[0]/sig[0])*(sigE[0]/sig[0])) : 0.0;
        double supE = (bkg[0] > 0 && bkg[j] > 0)
            ? sup * std::sqrt((bkgE[j]/bkg[j])*(bkgE[j]/bkg[j]) +
                              (bkgE[0]/bkg[0])*(bkgE[0]/bkg[0])) : 0.0;
        gRetAll->SetPoint(j, j + 0.5, ret); gRetAll->SetPointError(j, 0, retE);
        gSupAll->SetPoint(j, j + 0.5, sup); gSupAll->SetPointError(j, 0, supE);
    }
    auto* cRSall = new TCanvas("c_ret_sup_all",
        "retention vs suppression (all designs)", 950, 680);
    cRSall->SetBottomMargin(0.28);
    gRetAll->SetTitle("signal retention & background suppression;;ratio to nominal");
    gRetAll->SetMarkerStyle(20); gRetAll->SetMarkerColor(getP6(0));
    gRetAll->SetLineColor(getP6(0)); gRetAll->SetLineWidth(2);
    gSupAll->SetMarkerStyle(21); gSupAll->SetMarkerColor(getP6(1));
    gSupAll->SetLineColor(getP6(1)); gSupAll->SetLineWidth(2);
    gRetAll->GetXaxis()->SetLimits(0, nconfig);
    gRetAll->GetYaxis()->SetRangeUser(0, 1.3);
    gRetAll->Draw("AP"); gSupAll->Draw("P SAME");
    auto* lone = new TLine(0, 1, nconfig, 1);
    lone->SetLineStyle(2); lone->SetLineColor(kGray + 2); lone->Draw();
    { TLatex tx; tx.SetTextSize(0.026); tx.SetTextAngle(90); tx.SetTextAlign(32);
      for (int j = 0; j < nconfig; ++j)
          tx.DrawLatex(j + 0.5, -0.02, Form(" %s", T(j))); }
    auto* legRSall = new TLegend(0.50, 0.78, 0.88, 0.90); legRSall->SetTextSize(0.026);
    legRSall->AddEntry(gRetAll, "signal retention (cask 0 / nominal)", "lp");
    legRSall->AddEntry(gSupAll, "background suppression (casks 1-11 / nominal)", "lp");
    legRSall->Draw();

    // =======================================================================
    //  SIGNAL vs BACKGROUND spectra across designs -- (8)/(9)/(10).
    //  Nominal is solid black; other designs use the palette (dashed on wrap).
    // =======================================================================
    auto styleByConfig = [&](TH1* h, int j) {
        if (!h) return;
        StyleConfigLine(h, j);       // nominal (j==0) -> black; else palette
        h->SetDirectory(nullptr);
    };

    // shared y-range so signal and background panels are directly comparable
    double ymax = 1.0;
    for (int j = 0; j < nconfig; ++j) {
        if (sigSpec[j]) ymax = std::max(ymax, sigSpec[j]->GetMaximum());
        if (bkgSpec[j]) ymax = std::max(ymax, bkgSpec[j]->GetMaximum());
    }

    auto drawSpecOverlay = [&](std::vector<TH1*>& v, const char* title) {
        auto* leg = new TLegend(0.50, 0.55, 0.97, 0.90); leg->SetTextSize(0.022);
        bool first = true;
        for (int j = 0; j < nconfig; ++j) {
            if (!v[j]) continue;
            styleByConfig(v[j], j);
            v[j]->SetTitle(Form("%s;edep [MeV];counts", title));
            v[j]->SetMinimum(0.5); v[j]->SetMaximum(3.0 * ymax);
            v[j]->Draw(first ? "HIST" : "HIST SAME");
            first = false;
            leg->AddEntry(v[j], T(j), "l");
        }
        leg->Draw();
    };

    // // ---- (8) signal spectra, all configs -----------------------------------
    // auto* cSigSpec = new TCanvas("c_sigspec_allconfig",
    //     "cask-0 (signal) spectra, all configs", 1100, 700);
    // cSigSpec->SetLogy(); cSigSpec->SetGrid();
    // drawSpecOverlay(sigSpec, "Signal spectrum (cask 0)");
    // gPad->Modified(); gPad->Update();

    // // ---- (9) background spectra, all configs -------------------------------
    // auto* cBkgSpec = new TCanvas("c_bkgspec_allconfig",
    //     "casks 1+2+3 (background) spectra, all configs", 1100, 700);
    // cBkgSpec->SetLogy(); cBkgSpec->SetGrid();
    // drawSpecOverlay(bkgSpec, "Background spectrum (casks 1+2+3)");
    // gPad->Modified(); gPad->Update();

    // ---- (10) side-by-side signal | background, shared scale ---------------
    auto* cSB = new TCanvas("c_sig_bkg_spec",
        "signal vs background spectra", 1500, 700);
    cSB->Divide(2, 1, 0.001, 0.001);
    cSB->cd(1); gPad->SetLogy(); gPad->SetGrid();
    drawSpecOverlay(sigSpec, "Signal spectrum (cask 0)");
    cSB->cd(2); gPad->SetLogy(); gPad->SetGrid();
    drawSpecOverlay(bkgSpec, "Background spectrum (casks 1-11)");
    gPad->Modified(); gPad->Update();

    // =======================================================================
    //  (11)/(12)/(13) cask-0 signal-spot PROFILES: phi and z projections.
    //  Normalized to unit area so we compare NARROWNESS/shape, not raw flux.
    //  Gate lines mark |phi|<phiGate and |z|<zGate. Nominal is solid black.
    // =======================================================================
    auto styleProj = [&](TH1* h, int j) {
        StyleConfigLine(h, j);       // nominal (j==0) -> black; else palette
        h->SetStats(0);
        if (h->Integral() > 0) h->Scale(1.0 / h->Integral());  // unit area
    };
    auto vgateLines = [&](double g, double ymaxL) {
        //for (double x : { -g, +g }) {
        //    auto* l = new TLine(x, 0, x, ymaxL);
        //    l->SetLineColor(kGray + 2); l->SetLineStyle(3); l->SetLineWidth(1);
        //    l->Draw();
        //}
    };

    // // ---- (11) phi projection overlay ---------------------------------------
    // auto* cPhi = new TCanvas("c_proj_phi",
    //     "cask-0 signal profile in #phi (all designs)", 1000, 680);
    // cPhi->SetGrid();
    // auto* legP = new TLegend(0.62, 0.55, 0.97, 0.90); legP->SetTextSize(0.022);
     double yPhiMax = 0;
     for (int j = 0; j < nconfig; ++j) if (phiProj[j]) {
         styleProj(phiProj[j], j);
         yPhiMax = std::max(yPhiMax, phiProj[j]->GetMaximum());
     }
    // bool firstP = true;
    // for (int j = 0; j < nconfig; ++j) if (phiProj[j]) {
    //     phiProj[j]->SetTitle("cask-0 emission profile in #phi;"
    //                          "#phi [rad];fraction / bin");
    //     phiProj[j]->SetMaximum(1.25 * yPhiMax); phiProj[j]->SetMinimum(0);
    //     phiProj[j]->Draw(firstP ? "HIST" : "HIST SAME");
    //     firstP = false;
    //     legP->AddEntry(phiProj[j], T(j), "l");
    // }
    // vgateLines(phiGate, 1.25 * yPhiMax);
    // legP->Draw();
    // gPad->Modified(); gPad->Update();

    // // ---- (12) z projection overlay -----------------------------------------
    // auto* cZ = new TCanvas("c_proj_z",
    //     "cask-0 signal profile in z (all designs)", 1000, 680);
    // cZ->SetGrid();
    // auto* legZ = new TLegend(0.62, 0.55, 0.97, 0.90); legZ->SetTextSize(0.022);
    double yZMax = 0;
    for (int j = 0; j < nconfig; ++j) if (zProj[j]) {
        styleProj(zProj[j], j);
        yZMax = std::max(yZMax, zProj[j]->GetMaximum());
    }
    // bool firstZ = true;
    // for (int j = 0; j < nconfig; ++j) if (zProj[j]) {
    //     zProj[j]->SetTitle("cask-0 emission profile in z;"
    //                        "z [mm];fraction / bin");
    //     zProj[j]->SetMaximum(1.25 * yZMax); zProj[j]->SetMinimum(0);
    //     zProj[j]->Draw(firstZ ? "HIST" : "HIST SAME");
    //     firstZ = false;
    //     legZ->AddEntry(zProj[j], T(j), "l");
    // }
    // vgateLines(zGate, 1.25 * yZMax);
    // legZ->Draw();
    // gPad->Modified(); gPad->Update();

    // ---- (13) side-by-side phi | z -----------------------------------------
    auto* cPZ = new TCanvas("c_proj_phi_z",
        "cask-0 signal profiles: #phi | z", 1500, 680);
    cPZ->Divide(2, 1, 0.001, 0.001);
    cPZ->cd(1); gPad->SetGrid();
    { bool first = true;
      for (int j = 0; j < nconfig; ++j) if (phiProj[j]) {
          phiProj[j]->Draw(first ? "HIST" : "HIST SAME"); first = false; }
      vgateLines(phiGate, 1.25 * yPhiMax); }
    cPZ->cd(2); gPad->SetGrid();
    { bool first = true;
      for (int j = 0; j < nconfig; ++j) if (zProj[j]) {
          zProj[j]->Draw(first ? "HIST" : "HIST SAME"); first = false; }
      vgateLines(zGate, 1.25 * yZMax); }
    gPad->Modified(); gPad->Update();
}
#endif

