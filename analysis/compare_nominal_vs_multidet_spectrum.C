#ifndef COMPARE_NOMINAL_VS_MULTIDET_SPECTRUM_C
#define COMPARE_NOMINAL_VS_MULTIDET_SPECTRUM_C 1

// compare_nominal_vs_multidet_spectrum.C
//
// Purpose:
//   Compare the detector spectrum from:
//
//     1) data/nominal
//        - single physical detector at the nominal position
//
//     2) data/nominal-all-positions
//        - multidet run, but selecting only the detector whose position matches
//          the nominal detector.
//
//   The nominal-equivalent multidet detector is det=96 because
//     z order:   -1050,-750,-450,-150,150,450,750,1050  -> z=150 is iz=4
//     phi order: 0,15,...,345                            -> phi=0 is iphi=0
//     det = iz*24 + iphi = 4*24 + 0 = 96
//
//   The script sums over all 12 casks by default, because that is the measured
//   all-cask field. You can optionally restrict to one cask with caskSelect.
//
// Main entry point:
//   root -l 'compare_nominal_vs_multidet_spectrum.C("gamma")'
//
// Examples:
//   root -l 'compare_nominal_vs_multidet_spectrum.C("gamma")'
//   root -l 'compare_nominal_vs_multidet_spectrum.C("neutron")'
//
//   // compare cask 0 only:
//   root -l 'compare_nominal_vs_multidet_spectrum.C("gamma", kCLYC, 2, "data/nominal", "data/nominal-all-positions", 96, 10, false, true, 0)'
//
// Notes:
//   - By default, this matches the existing analysis convention: unweighted
//     histogram filling. The current analysis/analyze_cask.C fills spectra with
//     h->Fill(edep), not h->Fill(edep, weight).
//   - If you want to explicitly test weighted spectra, set useWeights=true.
//   - The residual convention is:
//         diff    = multidet_det96 - nominal
//         percent = 100 * (multidet_det96 - nominal) / nominal
//

#include "geometry_constraints.h"
#include "dcs_cuts.h"
#include "dcs_sources.h"
#include "dcs_setup.h"
#include "style.h"

#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPad.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

using namespace geo;
using namespace dcs;
using namespace PetroffPalette;

namespace compare_nominal_md_detail {

struct FileSlot {
    int         cask = -1;
    int         globalFuel = -1;
    int         sourceIndex = -1;
    std::string sourceTag;
    std::string path;
};

struct SpectrumResult {
    TH1D* h = nullptr;

    Long64_t nFilesFound = 0;
    Long64_t nFilesOpened = 0;
    Long64_t nFilesMissing = 0;
    Long64_t nHitsTotal = 0;
    Long64_t nHitsAcceptedPid = 0;
    Long64_t nHitsAcceptedDet = 0;
    Long64_t nHitsFilled = 0;

    Long64_t nFilesWithoutDetBranch = 0;
    Long64_t nFilesWithoutWeightBranch = 0;
};

static std::string ParticleLatex(const std::string& part)
{
    if (part == "gamma")   return "#gamma";
    if (part == "neutron") return "neutron";
    return "all";
}

static std::string CaskLabel(int caskSelect)
{
    if (caskSelect < 0) return "all casks 0-11";
    return Form("cask %d only", caskSelect);
}

static bool SourceWanted(const std::string& sourceTag,
                         const std::string& requestedSource)
{
    if (requestedSource.empty()) return true;
    if (requestedSource == "all") return true;
    return sourceTag == requestedSource;
}

static void DiscoverCaskFiles(const std::string& caskDir,
                              int caskIndex,
                              const std::string& requestedSource,
                              std::vector<FileSlot>& out)
{
    if (!fs::exists(caskDir)) {
        std::cerr << "[warn] input cask directory does not exist: "
                  << caskDir << "\n";
        return;
    }

    static const std::regex reFuel("globalFuel0*([0-9]+)");
    const auto& sources = Sources();

    for (const auto& e : fs::directory_iterator(caskDir)) {
        if (!e.is_regular_file()) continue;
        if (e.path().extension() != ".root") continue;

        const std::string name = e.path().filename().string();

        if (name.rfind("detector-response", 0) != 0) continue;

        std::smatch m;
        if (!std::regex_search(name, m, reFuel)) continue;

        const int g = std::stoi(m[1].str());
        if (g < 0 || g >= kNAssem) continue;

        int si = -1;
        for (size_t k = 0; k < sources.size(); ++k) {
            if (name.find(sources[k].tag) != std::string::npos) {
                si = static_cast<int>(k);
                break;
            }
        }

        if (si < 0) continue;

        const std::string tag = sources[si].tag;
        if (!SourceWanted(tag, requestedSource)) continue;

        FileSlot slot;
        slot.cask = caskIndex;
        slot.globalFuel = g;
        slot.sourceIndex = si;
        slot.sourceTag = tag;
        slot.path = e.path().string();

        out.push_back(slot);
    }
}

static std::vector<FileSlot> DiscoverDatasetFiles(const std::string& baseDir,
                                                  int caskSelect,
                                                  const std::string& requestedSource)
{
    std::vector<FileSlot> files;

    const int c0 = (caskSelect >= 0) ? caskSelect : 0;
    const int c1 = (caskSelect >= 0) ? caskSelect : (kNCask - 1);

    for (int ic = c0; ic <= c1; ++ic) {
        const std::string caskDir = Form("%s/cask%d", baseDir.c_str(), ic);
        DiscoverCaskFiles(caskDir, ic, requestedSource, files);
    }

    std::sort(files.begin(), files.end(),
              [](const FileSlot& a, const FileSlot& b) {
                  if (a.cask != b.cask) return a.cask < b.cask;
                  if (a.globalFuel != b.globalFuel) return a.globalFuel < b.globalFuel;
                  return a.sourceIndex < b.sourceIndex;
              });

    return files;
}

static void PrintOneLineProgress(const char* label,
                                 Long64_t i,
                                 Long64_t n,
                                 const FileSlot& slot)
{
    std::cout << "\33[2K\r"
              << "[" << label << "] file " << i << "/" << n
              << " | cask=" << slot.cask
              << " | globalFuel=" << Form("%02d", slot.globalFuel)
              << " | source=" << slot.sourceTag
              << std::flush;
}

static bool OpenHitsTree(const std::string& filename,
                         TFile*& f,
                         TTree*& hits)
{
    f = nullptr;
    hits = nullptr;

    const Int_t prev = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;
    f = TFile::Open(filename.c_str(), "READ");
    gErrorIgnoreLevel = prev;

    if (!f || f->IsZombie() || f->TestBit(TFile::kRecovered)) {
        if (f) {
            f->Close();
            delete f;
            f = nullptr;
        }
        return false;
    }

    hits = dynamic_cast<TTree*>(f->Get("hits"));
    if (!hits) {
        f->Close();
        delete f;
        f = nullptr;
        return false;
    }

    return true;
}

static void ConfigureStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetNumberContours(255);
}

static SpectrumResult BuildSpectrumFromRawFiles(const std::string& label,
                                                const std::string& baseDir,
                                                const std::string& part,
                                                DetectorType detectorType,
                                                int level,
                                                int detFilter,
                                                bool requireDetBranch,
                                                bool useWeights,
                                                int caskSelect,
                                                const std::string& sourceTag,
                                                int nBins,
                                                double eMin,
                                                double eMax)
{
    SpectrumResult R;

    const std::vector<FileSlot> files =
        DiscoverDatasetFiles(baseDir, caskSelect, sourceTag);

    R.nFilesFound = files.size();

    R.h = new TH1D(Form("h_%s_%s", label.c_str(), part.c_str()),
                   Form("%s spectrum;E_{dep} [MeV];counts / bin", label.c_str()),
                   nBins, eMin, eMax);
    R.h->Sumw2();
    R.h->SetDirectory(nullptr);

    CutConfig cfg;
    cfg.det = detectorType;
    cfg.part = ParseParticle(part);
    cfg.level = static_cast<CutLevel>(level);

    std::cout << "[" << label << "] discovered " << files.size()
              << " files under " << baseDir
              << " (" << CaskLabel(caskSelect) << ", source="
              << (sourceTag.empty() ? "all" : sourceTag) << ")\n";

    Long64_t idx = 0;
    for (const auto& slot : files) {
        ++idx;
        PrintOneLineProgress(label.c_str(), idx, files.size(), slot);

        TFile* f = nullptr;
        TTree* hits = nullptr;

        if (!OpenHitsTree(slot.path, f, hits)) {
            ++R.nFilesMissing;
            continue;
        }

        ++R.nFilesOpened;

        double h_edep = 0.0;
        double h_pid  = 0.0;
        double h_w    = 1.0;
        Int_t  h_det  = -999999;

        hits->SetBranchStatus("*", 0);

        if (!hits->GetBranch("edep") || !hits->GetBranch("pid")) {
            std::cerr << "\n[warn] required branch edep/pid missing in "
                      << slot.path << "; skipping\n";
            f->Close();
            delete f;
            continue;
        }

        hits->SetBranchStatus("edep", 1);
        hits->SetBranchStatus("pid",  1);

        const bool haveDet = (hits->GetBranch("det") != nullptr);
        const bool haveW   = (hits->GetBranch("weight") != nullptr);

        if (haveDet) {
            hits->SetBranchStatus("det", 1);
        } else {
            ++R.nFilesWithoutDetBranch;
        }

        if (useWeights && haveW) {
            hits->SetBranchStatus("weight", 1);
        } else if (useWeights && !haveW) {
            ++R.nFilesWithoutWeightBranch;
        }

        hits->SetBranchAddress("edep", &h_edep);
        hits->SetBranchAddress("pid",  &h_pid);

        if (haveDet) {
            hits->SetBranchAddress("det", &h_det);
        }

        if (useWeights && haveW) {
            hits->SetBranchAddress("weight", &h_w);
        }

        const Long64_t n = hits->GetEntries();

        for (Long64_t i = 0; i < n; ++i) {
            hits->GetEntry(i);

            ++R.nHitsTotal;

            const long long pid = std::llround(h_pid);

            if (!PassSpectrum(cfg, pid)) continue;
            ++R.nHitsAcceptedPid;

            bool passDet = true;

            if (detFilter >= 0) {
                if (!haveDet) {
                    passDet = !requireDetBranch;
                } else {
                    passDet = (h_det == detFilter);
                }
            }

            if (!passDet) continue;
            ++R.nHitsAcceptedDet;

            const double w = (useWeights && haveW) ? h_w : 1.0;
            R.h->Fill(h_edep, w);
            ++R.nHitsFilled;
        }

        f->Close();
        delete f;
    }

    std::cout << "\33[2K\r"
              << "[" << label << "] completed: opened "
              << R.nFilesOpened << "/" << R.nFilesFound
              << " files, filled " << R.nHitsFilled << " hits"
              << "\n";

    if (requireDetBranch && R.nFilesWithoutDetBranch > 0) {
        std::cerr << "[warn] " << label << ": "
                  << R.nFilesWithoutDetBranch
                  << " files had no hits.det branch while det filtering was required\n";
    }

    if (useWeights && R.nFilesWithoutWeightBranch > 0) {
        std::cerr << "[warn] " << label << ": "
                  << R.nFilesWithoutWeightBranch
                  << " files had no hits.weight branch; those files were filled unweighted\n";
    }

    return R;
}

static TH1D* MakeDifferenceHist(const TH1D* hNom,
                                const TH1D* hMD,
                                const char* name)
{
    auto* h = dynamic_cast<TH1D*>(hMD->Clone(name));
    h->SetDirectory(nullptr);
    h->Reset("ICES");
    h->SetTitle("absolute residual;E_{dep} [MeV];multidet det96 - nominal");

    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double m  = hMD ->GetBinContent(b);
        const double n  = hNom->GetBinContent(b);
        const double em = hMD ->GetBinError(b);
        const double en = hNom->GetBinError(b);

        h->SetBinContent(b, m - n);
        h->SetBinError  (b, std::sqrt(em*em + en*en));
    }

    return h;
}

static TH1D* MakePercentResidualHist(const TH1D* hNom,
                                     const TH1D* hMD,
                                     const char* name)
{
    auto* h = dynamic_cast<TH1D*>(hMD->Clone(name));
    h->SetDirectory(nullptr);
    h->Reset("ICES");
    h->SetTitle("percent residual;E_{dep} [MeV];100 #times (multidet det96 - nominal) / nominal [%]");

    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double m  = hMD ->GetBinContent(b);
        const double n  = hNom->GetBinContent(b);
        const double em = hMD ->GetBinError(b);
        const double en = hNom->GetBinError(b);

        if (n <= 0.0) {
            h->SetBinContent(b, 0.0);
            h->SetBinError  (b, 0.0);
            continue;
        }

        const double pct = 100.0 * (m - n) / n;

        // f = 100 * (M/N - 1)
        // df/dM = 100/N
        // df/dN = -100*M/N^2
        const double epct = 100.0 * std::sqrt(
            (em / n) * (em / n) +
            (m * en / (n * n)) * (m * en / (n * n))
        );

        h->SetBinContent(b, pct);
        h->SetBinError  (b, epct);
    }

    return h;
}

static TH1D* MakeRatioHist(const TH1D* hNom,
                           const TH1D* hMD,
                           const char* name)
{
    auto* h = dynamic_cast<TH1D*>(hMD->Clone(name));
    h->SetDirectory(nullptr);
    h->Reset("ICES");
    h->SetTitle("ratio;E_{dep} [MeV];multidet det96 / nominal");

    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double m  = hMD ->GetBinContent(b);
        const double n  = hNom->GetBinContent(b);
        const double em = hMD ->GetBinError(b);
        const double en = hNom->GetBinError(b);

        if (n <= 0.0) {
            h->SetBinContent(b, 0.0);
            h->SetBinError  (b, 0.0);
            continue;
        }

        const double r = m / n;
        double er = 0.0;

        if (m > 0.0) {
            er = r * std::sqrt((em/m)*(em/m) + (en/n)*(en/n));
        } else {
            er = em / n;
        }

        h->SetBinContent(b, r);
        h->SetBinError  (b, er);
    }

    return h;
}

static double IntegralError(const TH1D* h, int b1, int b2)
{
    double e2 = 0.0;
    for (int b = b1; b <= b2; ++b) {
        const double e = h->GetBinError(b);
        e2 += e * e;
    }
    return std::sqrt(e2);
}

static void PrintSummary(const TH1D* hNom,
                         const TH1D* hMD,
                         const TH1D* hPct,
                         const std::string& part,
                         DetectorType detectorType,
                         int level)
{
    const int b1 = 1;
    const int b2 = hNom->GetNbinsX();

    const double IN = hNom->Integral(b1, b2);
    const double IM = hMD ->Integral(b1, b2);
    const double EN = IntegralError(hNom, b1, b2);
    const double EM = IntegralError(hMD,  b1, b2);

    const double diff = IM - IN;
    const double ediff = std::sqrt(EN*EN + EM*EM);

    const double pct = (IN > 0.0) ? 100.0 * diff / IN : 0.0;
    const double epct = (IN > 0.0)
        ? 100.0 * std::sqrt((EM/IN)*(EM/IN) +
                            (IM*EN/(IN*IN))*(IM*EN/(IN*IN)))
        : 0.0;

    double chi2 = 0.0;
    int ndf = 0;

    for (int b = 1; b <= hNom->GetNbinsX(); ++b) {
        const double n  = hNom->GetBinContent(b);
        const double m  = hMD ->GetBinContent(b);
        const double en = hNom->GetBinError(b);
        const double em = hMD ->GetBinError(b);
        const double e2 = en*en + em*em;

        if (e2 <= 0.0) continue;

        const double d = m - n;
        chi2 += d*d / e2;
        ++ndf;
    }

    double eLo = 0.0;
    double eHi = 0.0;
    CutConfig cfg{detectorType, ParseParticle(part), static_cast<CutLevel>(level)};
    EnergyWindow(cfg, eLo, eHi);

    const int r1 = hNom->GetXaxis()->FindBin(eLo + 1e-9);
    const int r2 = hNom->GetXaxis()->FindBin(eHi - 1e-9);

    const double RIN = hNom->Integral(r1, r2);
    const double RIM = hMD ->Integral(r1, r2);
    const double REN = IntegralError(hNom, r1, r2);
    const double REM = IntegralError(hMD,  r1, r2);

    const double rdiff = RIM - RIN;
    const double rediff = std::sqrt(REN*REN + REM*REM);
    const double rpct = (RIN > 0.0) ? 100.0 * rdiff / RIN : 0.0;
    const double repct = (RIN > 0.0)
        ? 100.0 * std::sqrt((REM/RIN)*(REM/RIN) +
                            (RIM*REN/(RIN*RIN))*(RIM*REN/(RIN*RIN)))
        : 0.0;

    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << " Nominal single-detector vs multidet det96 spectrum summary\n";
    std::cout << "============================================================\n";
    std::cout << " particle: " << part << "\n";
    std::cout << " detector: " << (detectorType == kCLYC ? "CLYC" : "plastic") << "\n";
    std::cout << " cut level for spectrum PID selection: " << level << "\n";
    std::cout << "\n";
    std::cout << " Full plotted range:\n";
    std::cout << "   nominal integral        = " << IN << " +/- " << EN << "\n";
    std::cout << "   multidet det96 integral = " << IM << " +/- " << EM << "\n";
    std::cout << "   difference MD - nominal = " << diff << " +/- " << ediff << "\n";
    std::cout << "   percent difference      = " << pct << " +/- " << epct << " %\n";
    std::cout << "   chi2/ndf                = "
              << chi2 << " / " << ndf
              << " = " << (ndf > 0 ? chi2 / ndf : 0.0) << "\n";
    std::cout << "\n";
    std::cout << " Default ROI from dcs_cuts.h EnergyWindow(): ["
              << eLo << ", " << eHi << "] MeV\n";
    std::cout << "   nominal ROI integral        = " << RIN << " +/- " << REN << "\n";
    std::cout << "   multidet det96 ROI integral = " << RIM << " +/- " << REM << "\n";
    std::cout << "   ROI difference MD - nominal = " << rdiff << " +/- " << rediff << "\n";
    std::cout << "   ROI percent difference      = " << rpct << " +/- " << repct << " %\n";
    std::cout << "============================================================\n";
    std::cout << "\n";
}

static void StyleSpectra(TH1D* hNom, TH1D* hMD)
{
    hNom->SetLineColor(kBlack);
    hNom->SetMarkerColor(kBlack);
    hNom->SetMarkerStyle(20);
    hNom->SetMarkerSize(0.55);
    hNom->SetLineWidth(2);

    hMD->SetLineColor(getP6(1));
    hMD->SetMarkerColor(getP6(1));
    hMD->SetMarkerStyle(24);
    hMD->SetMarkerSize(0.55);
    hMD->SetLineWidth(2);
}

static void StyleResiduals(TH1D* hDiff, TH1D* hPct, TH1D* hRatio)
{
    hDiff->SetLineColor(getP6(2));
    hDiff->SetMarkerColor(getP6(2));
    hDiff->SetMarkerStyle(20);
    hDiff->SetMarkerSize(0.45);

    hPct->SetLineColor(getP6(4));
    hPct->SetMarkerColor(getP6(4));
    hPct->SetMarkerStyle(20);
    hPct->SetMarkerSize(0.45);

    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.45);
}

static TCanvas* DrawComparisonCanvas(TH1D* hNom,
                                     TH1D* hMD,
                                     TH1D* hDiff,
                                     TH1D* hPct,
                                     const std::string& canvasName,
                                     const std::string& part,
                                     int mdDetId,
                                     bool logY,
                                     int caskSelect,
                                     const std::string& sourceTag)
{
    auto* c = new TCanvas(canvasName.c_str(),
                          "nominal vs multidet nominal-position spectrum",
                          1200, 950);

    auto* p1 = new TPad("p_spectra", "spectra", 0.0, 0.46, 1.0, 1.0);
    auto* p2 = new TPad("p_diff",    "difference", 0.0, 0.23, 1.0, 0.46);
    auto* p3 = new TPad("p_pct",     "percent residual", 0.0, 0.0, 1.0, 0.23);

    p1->SetBottomMargin(0.015);
    p2->SetTopMargin(0.04);
    p2->SetBottomMargin(0.02);
    p3->SetTopMargin(0.04);
    p3->SetBottomMargin(0.32);

    p1->Draw();
    p2->Draw();
    p3->Draw();

    p1->cd();
    if (logY) gPad->SetLogy();
    gPad->SetGrid();

    const double maxY = std::max(hNom->GetMaximum(), hMD->GetMaximum());
    const double minY = logY ? 0.5 : 0.0;

    hNom->SetTitle(Form("%s spectrum comparison: nominal single detector vs multidet det%d;"
                        "E_{dep} [MeV];counts / bin",
                        ParticleLatex(part).c_str(), mdDetId));
    hNom->SetMinimum(minY);
    hNom->SetMaximum((logY ? 20.0 : 1.35) * std::max(1.0, maxY));

    hNom->GetXaxis()->SetLabelSize(0.0);
    hNom->GetYaxis()->SetTitleSize(0.055);
    hNom->GetYaxis()->SetLabelSize(0.045);
    hNom->GetYaxis()->SetTitleOffset(0.90);

    hNom->Draw("E1");
    hMD->Draw("E1 SAME");

    auto* leg = new TLegend(0.50, 0.68, 0.92, 0.90);
    leg->SetTextSize(0.035);
    leg->AddEntry(hNom, "data/nominal: single detector", "lep");
    leg->AddEntry(hMD,  Form("data/nominal-all-positions: det %d only", mdDetId), "lep");
    leg->Draw();

    TLatex tx;
    tx.SetNDC(true);
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.13, 0.86, Form("%s", CaskLabel(caskSelect).c_str()));
    tx.DrawLatex(0.13, 0.80, Form("source selection: %s",
                                  sourceTag.empty() ? "all" : sourceTag.c_str()));

    p2->cd();
    gPad->SetGrid();

    double dmax = 0.0;
    for (int b = 1; b <= hDiff->GetNbinsX(); ++b) {
        dmax = std::max(dmax,
                        std::abs(hDiff->GetBinContent(b)) + hDiff->GetBinError(b));
    }
    if (dmax <= 0.0) dmax = 1.0;

    hDiff->SetMinimum(-1.25 * dmax);
    hDiff->SetMaximum(+1.25 * dmax);
    hDiff->GetXaxis()->SetLabelSize(0.0);
    hDiff->GetYaxis()->SetTitleSize(0.10);
    hDiff->GetYaxis()->SetLabelSize(0.085);
    hDiff->GetYaxis()->SetTitleOffset(0.43);
    hDiff->GetYaxis()->SetNdivisions(505);
    hDiff->Draw("E1");

    auto* l0a = new TLine(hDiff->GetXaxis()->GetXmin(), 0.0,
                          hDiff->GetXaxis()->GetXmax(), 0.0);
    l0a->SetLineColor(kGray + 2);
    l0a->SetLineStyle(2);
    l0a->Draw();

    p3->cd();
    gPad->SetGrid();

    double pmax = 0.0;
    for (int b = 1; b <= hPct->GetNbinsX(); ++b) {
        const double y = hPct->GetBinContent(b);
        const double e = hPct->GetBinError(b);
        if (e <= 0.0 && y == 0.0) continue;
        pmax = std::max(pmax, std::abs(y) + e);
    }
    if (pmax <= 0.0) pmax = 1.0;

    // Avoid one huge empty/noisy bin setting an unreadably large scale.
    // This only affects drawing range, not the histogram contents.
    pmax = std::min(pmax, 100.0);

    hPct->SetMinimum(-1.25 * pmax);
    hPct->SetMaximum(+1.25 * pmax);
    hPct->GetXaxis()->SetTitleSize(0.13);
    hPct->GetXaxis()->SetLabelSize(0.11);
    hPct->GetYaxis()->SetTitleSize(0.10);
    hPct->GetYaxis()->SetLabelSize(0.085);
    hPct->GetYaxis()->SetTitleOffset(0.43);
    hPct->GetYaxis()->SetNdivisions(505);
    hPct->Draw("E1");

    auto* l0b = new TLine(hPct->GetXaxis()->GetXmin(), 0.0,
                          hPct->GetXaxis()->GetXmax(), 0.0);
    l0b->SetLineColor(kGray + 2);
    l0b->SetLineStyle(2);
    l0b->Draw();

    c->cd();
    c->Modified();
    c->Update();

    return c;
}

static TCanvas* DrawRatioCanvas(TH1D* hRatio,
                                const std::string& canvasName,
                                int mdDetId)
{
    auto* c = new TCanvas(canvasName.c_str(),
                          "multidet det96 / nominal ratio",
                          1100, 550);
    c->SetGrid();

    hRatio->SetTitle(Form("spectrum ratio: multidet det%d / nominal;"
                          "E_{dep} [MeV];ratio", mdDetId));
    hRatio->SetMinimum(0.5);
    hRatio->SetMaximum(1.5);
    hRatio->Draw("E1");

    auto* one = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0,
                          hRatio->GetXaxis()->GetXmax(), 1.0);
    one->SetLineColor(kGray + 2);
    one->SetLineStyle(2);
    one->SetLineWidth(2);
    one->Draw();

    c->Modified();
    c->Update();
    return c;
}

} // namespace compare_nominal_md_detail



void compare_nominal_vs_multidet_spectrum(
    std::string part = "gamma",
    DetectorType detectorType = kCLYC,
    int level = 2,
    std::string nominalDir = "data/nominal",
    std::string multidetDir = "data/nominal-all-positions",
    int multidetNominalDetId = 96,
    int rebin = 10,
    bool useWeights = false,
    bool logY = true,
    int caskSelect = -1,
    std::string sourceTag = "all",
    std::string outFile = "")
{
    using namespace compare_nominal_md_detail;

    ConfigureStyle();

    if (part != "gamma" && part != "neutron" && part != "all") {
        std::cerr << "[err] part must be one of: gamma, neutron, all\n";
        return;
    }

    if (sourceTag == "all") sourceTag.clear();

    if (outFile.empty()) {
        outFile = Form("root_output/compare_nominal_vs_multidet_det%03d_%s.root",
                       multidetNominalDetId, part.c_str());
    }

    fs::create_directories(fs::path(outFile).parent_path());

    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << " Comparing nominal single-detector spectrum to multidet\n";
    std::cout << "============================================================\n";
    std::cout << " nominal dir:       " << nominalDir << "\n";
    std::cout << " multidet dir:      " << multidetDir << "\n";
    std::cout << " multidet det id:   " << multidetNominalDetId << "\n";
    std::cout << " particle:          " << part << "\n";
    std::cout << " detector type:     " << (detectorType == kCLYC ? "CLYC" : "plastic") << "\n";
    std::cout << " cut level:         " << level << "\n";
    std::cout << " cask selection:    " << CaskLabel(caskSelect) << "\n";
    std::cout << " source selection:  " << (sourceTag.empty() ? "all" : sourceTag) << "\n";
    std::cout << " use weights:       " << (useWeights ? "true" : "false") << "\n";
    std::cout << " rebin:             " << rebin << "\n";
    std::cout << " output file:       " << outFile << "\n";
    std::cout << "============================================================\n\n";

    const int nBins = kNbE;
    const double eMin = kEmin;
    const double eMax = kEmax;

    // Nominal single-detector run:
    //   no detector filter by default, because there is only one detector.
    SpectrumResult nominal = BuildSpectrumFromRawFiles(
        "nominal",
        nominalDir,
        part,
        detectorType,
        level,
        -1,                // no det filtering
        false,             // do not require det branch
        useWeights,
        caskSelect,
        sourceTag,
        nBins,
        eMin,
        eMax);

    // Multidet run:
    //   require hits.det and select only det=96 by default.
    SpectrumResult multidet = BuildSpectrumFromRawFiles(
        "multidet_det96",
        multidetDir,
        part,
        detectorType,
        level,
        multidetNominalDetId,
        true,              // require det branch for meaningful comparison
        useWeights,
        caskSelect,
        sourceTag,
        nBins,
        eMin,
        eMax);

    if (!nominal.h || !multidet.h) {
        std::cerr << "[err] failed to build one or both spectra\n";
        return;
    }

    if (rebin > 1) {
        nominal.h->Rebin(rebin);
        multidet.h->Rebin(rebin);
    }

    StyleSpectra(nominal.h, multidet.h);

    auto* hDiff = MakeDifferenceHist(nominal.h, multidet.h, "h_residual_difference");
    auto* hPct  = MakePercentResidualHist(nominal.h, multidet.h, "h_residual_percent");
    auto* hRat  = MakeRatioHist(nominal.h, multidet.h, "h_ratio_multidet_over_nominal");

    StyleResiduals(hDiff, hPct, hRat);

    PrintSummary(nominal.h, multidet.h, hPct,
                 part, detectorType, level);

    auto* cMain = DrawComparisonCanvas(
        nominal.h,
        multidet.h,
        hDiff,
        hPct,
        Form("c_compare_nominal_vs_multidet_det%03d_%s",
             multidetNominalDetId, part.c_str()),
        part,
        multidetNominalDetId,
        logY,
        caskSelect,
        sourceTag);

    auto* cRatio = DrawRatioCanvas(
        hRat,
        Form("c_ratio_nominal_vs_multidet_det%03d_%s",
             multidetNominalDetId, part.c_str()),
        multidetNominalDetId);

    TFile* fout = TFile::Open(outFile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[err] could not create output file " << outFile << "\n";
        return;
    }

    nominal.h->Write("h_spectrum_nominal_single_detector");
    multidet.h->Write(Form("h_spectrum_multidet_det%03d", multidetNominalDetId));
    hDiff->Write();
    hPct->Write();
    hRat->Write();

    cMain->Write();
    cRatio->Write();

    TParameter<Int_t>("multidet_nominal_detector_id", multidetNominalDetId).Write();
    TParameter<Int_t>("rebin", rebin).Write();
    TParameter<Int_t>("use_weights", useWeights ? 1 : 0).Write();
    TParameter<Int_t>("cask_select", caskSelect).Write();

    TParameter<Long64_t>("nominal_files_found", nominal.nFilesFound).Write();
    TParameter<Long64_t>("nominal_files_opened", nominal.nFilesOpened).Write();
    TParameter<Long64_t>("nominal_hits_total", nominal.nHitsTotal).Write();
    TParameter<Long64_t>("nominal_hits_filled", nominal.nHitsFilled).Write();

    TParameter<Long64_t>("multidet_files_found", multidet.nFilesFound).Write();
    TParameter<Long64_t>("multidet_files_opened", multidet.nFilesOpened).Write();
    TParameter<Long64_t>("multidet_hits_total", multidet.nHitsTotal).Write();
    TParameter<Long64_t>("multidet_hits_filled", multidet.nHitsFilled).Write();

    fout->Close();

    const std::string pngMain =
        Form("root_output/compare_nominal_vs_multidet_det%03d_%s.png",
             multidetNominalDetId, part.c_str());
    const std::string pdfMain =
        Form("root_output/compare_nominal_vs_multidet_det%03d_%s.pdf",
             multidetNominalDetId, part.c_str());

    cMain->SaveAs(pngMain.c_str());
    cMain->SaveAs(pdfMain.c_str());

    std::cout << "[ok] wrote comparison ROOT file: " << outFile << "\n";
    std::cout << "[ok] wrote comparison PNG:       " << pngMain << "\n";
    std::cout << "[ok] wrote comparison PDF:       " << pdfMain << "\n";
}

#endif

