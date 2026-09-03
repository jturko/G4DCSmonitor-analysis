#ifndef PLOT_DESIGN_COMP
#define PLOT_DESIGN_COMP 1

#include "analyze_cask.C"
#include "dcs_geometry.h"
#include "dcs_cuts.h"
#include "dcs_multidet.h"
#include "utils.h"
#include "style.h"

#include <THStack.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <filesystem>
#include <algorithm>
#include <vector>

namespace fs = std::filesystem;
using namespace PetroffPalette;
using namespace dcs;

// New hall layout is axis-aligned -> no per-map rotation.
namespace { constexpr double kRotDeg = 0.0; }

static long CountInputs(const std::string& caskDir, fs::file_time_type& newest)
{
    std::vector<std::vector<std::string>> fp;
    DiscoverFiles(caskDir.c_str(), Sources(), fp);

    long n = 0;
    newest = fs::file_time_type::min();

    for (auto& row : fp) {
        for (auto& s : row) {
            if (!s.empty()) {
                ++n;
                std::error_code ec;
                auto t = fs::last_write_time(s, ec);
                if (!ec && t > newest) newest = t;
            }
        }
    }

    return n;
}

static bool HaveGoodFile(const std::string& fn)
{
    const Int_t pv = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;

    auto* t = TFile::Open(fn.c_str());
    const bool ok = t && !t->IsZombie() && !t->TestBit(TFile::kRecovered);

    if (t) {
        t->Close();
        delete t;
    }

    gErrorIgnoreLevel = pv;
    return ok;
}

static long CachedInputCount(const std::string& fn)
{
    const Int_t pv = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;

    auto* f = TFile::Open(fn.c_str());
    long v = -1;

    if (f && !f->IsZombie() && !f->TestBit(TFile::kRecovered)) {
        if (auto* p = (TParameter<Long64_t>*)f->Get("n_input_files")) {
            v = (long)p->GetVal();
        }
    }

    if (f) {
        f->Close();
        delete f;
    }

    gErrorIgnoreLevel = pv;
    return v;
}

static bool CacheHasObject(const std::string& fn, const char* obj)
{
    const Int_t pv = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;

    auto* f = TFile::Open(fn.c_str());
    bool ok = false;

    if (f && !f->IsZombie() && !f->TestBit(TFile::kRecovered)) {
        ok = (f->Get(obj) != nullptr);
    }

    if (f) {
        f->Close();
        delete f;
    }

    gErrorIgnoreLevel = pv;
    return ok;
}

struct DetectorPhiZSummary {
    TH2F* measuredTotal = nullptr;  // all casks, observable T = S+B
    TH2F* targetCask0   = nullptr;  // cask 0 only, S
    TH2F* background    = nullptr;  // casks 1-11, B
    TH2F* targetFrac    = nullptr;  // S/(S+B)
    TH2F* targetOverBkg = nullptr;  // S/B
    bool  ok            = false;
};

static DetectorPhiZSummary BuildDetectorPhiZSummary(const std::string& pathTag,
                                                    const std::string& part,
                                                    const std::string& label)
{
    DetectorPhiZSummary R;

    R.measuredTotal = MakeDetectorPhiZHist(
        Form("h2_detector_phiz_measured_total_%s", label.c_str()),
        Form("MEASURED total detector-position ROI counts, all casks, %s;"
             "detector #phi [deg];detector z [mm];T = cask0 + casks1-11",
             label.c_str()));

    R.targetCask0 = MakeDetectorPhiZHist(
        Form("h2_detector_phiz_target_cask0_%s", label.c_str()),
        Form("target-cask contribution, cask 0 only, %s;"
             "detector #phi [deg];detector z [mm];S = cask0 counts",
             label.c_str()));

    R.background = MakeDetectorPhiZHist(
        Form("h2_detector_phiz_background_casks1to11_%s", label.c_str()),
        Form("background-cask contribution, casks 1-11, %s;"
             "detector #phi [deg];detector z [mm];B = casks1-11 counts",
             label.c_str()));

    R.targetFrac = MakeDetectorPhiZHist(
        Form("h2_detector_phiz_target_fraction_%s", label.c_str()),
        Form("target-cask fraction, %s;"
             "detector #phi [deg];detector z [mm];S/(S+B)",
             label.c_str()));

    R.targetOverBkg = MakeDetectorPhiZHist(
        Form("h2_detector_phiz_target_over_background_%s", label.c_str()),
        Form("target/background ratio, %s;"
             "detector #phi [deg];detector z [mm];S/B",
             label.c_str()));

    for (int ic = 0; ic < kNCask; ++ic) {
        const std::string fn =
            Form("root_output/analysis_%s_%s_cask%d.root",
                 pathTag.c_str(), part.c_str(), ic);

        const Int_t pv = gErrorIgnoreLevel;
        gErrorIgnoreLevel = kFatal;
        auto* f = TFile::Open(fn.c_str());
        gErrorIgnoreLevel = pv;

        if (!f || f->IsZombie()) {
            if (f) {
                f->Close();
                delete f;
            }
            continue;
        }

        auto* h = dynamic_cast<TH2F*>(f->Get("h2_detector_phiz_counts_total"));
        if (!h) {
            f->Close();
            delete f;
            continue;
        }

        R.ok = true;

        for (int d = 0; d < kNMultidet; ++d) {
            const double v = GetDetectorPhiZBin(h, d);

            AddDetectorPhiZBin(R.measuredTotal, d, v);

            if (IsSignalCask(ic)) {
                AddDetectorPhiZBin(R.targetCask0, d, v);
            } else {
                AddDetectorPhiZBin(R.background, d, v);
            }
        }

        f->Close();
        delete f;
    }

    if (R.ok) {
        DivideDetectorPhiZ(R.targetCask0, R.measuredTotal, R.targetFrac,    0.0);
        DivideDetectorPhiZ(R.targetCask0, R.background,    R.targetOverBkg, 0.0);
    }

    return R;
}

static void DrawDetectorPhiZSummary(const DetectorPhiZSummary& R,
                                    const std::string& canvasTag)
{
    if (!R.ok) return;

    auto* c = new TCanvas(Form("c_detector_phiz_summary_%s", canvasTag.c_str()),
                          "detector-position all-cask/target/background summary",
                          1650, 950);

    c->Divide(3, 2, 0.001, 0.001);

    c->cd(1);
    gPad->SetRightMargin(0.15);
    R.measuredTotal->Draw("COLZ TEXT");

    c->cd(2);
    gPad->SetRightMargin(0.15);
    R.targetCask0->Draw("COLZ TEXT");

    c->cd(3);
    gPad->SetRightMargin(0.15);
    R.background->Draw("COLZ TEXT");

    c->cd(4);
    gPad->SetRightMargin(0.15);
    R.targetFrac->SetMinimum(0.0);
    R.targetFrac->SetMaximum(1.0);
    R.targetFrac->Draw("COLZ TEXT");

    c->cd(5);
    gPad->SetRightMargin(0.15);
    R.targetOverBkg->Draw("COLZ TEXT");

    c->cd(6);
    TLatex tx;
    tx.SetNDC(true);
    tx.SetTextSize(0.035);
    tx.DrawLatex(0.08, 0.82, "Detector-position summary");
    tx.DrawLatex(0.08, 0.68, "Pad 1: measured total T = all casks");
    tx.DrawLatex(0.08, 0.56, "Pad 2: target contribution S = cask 0");
    tx.DrawLatex(0.08, 0.44, "Pad 3: background B = casks 1-11");
    tx.DrawLatex(0.08, 0.32, "Pad 4: target fraction S/(S+B)");
    tx.DrawLatex(0.08, 0.20, "Pad 5: target/background S/B");

    gPad->Modified();
    gPad->Update();
}

// Produce/refresh the per-cask analysis files and return the grand-sum spectrum
// all casks for one config, as an independent clone.
TH1* plot_heatmap(fs::path p = "data/roomReturn/hall_3cluster/2cm_Pb",
                  std::string typePart = "gamma",
                  DetectorType typeDet = kCLYC,
                  int level = 2,
                  bool slimAnalysis = false)
{
    std::string pathTag = p.lexically_relative("data").generic_string();
    if (pathTag.empty() || pathTag == ".") pathTag = p.filename().string();

    std::string dir = pathTag;
    std::replace(dir.begin(), dir.end(), '/', '_');

    const int ncask = kNCask;

    std::vector<TH2Poly*> maps(ncask, nullptr);
    std::vector<TH1*>     spec(ncask, nullptr);
    TH1* specSum = nullptr;

    std::cout << "[plot_heatmap] config=" << pathTag
              << " particle=" << typePart
              << " mode=" << (slimAnalysis ? "slim" : "full")
              << "\n";

    for (int ic = 0; ic < ncask; ++ic) {
        const std::string caskDir = Form("%s/cask%d", p.string().c_str(), ic);

        const std::string fnameG =
            Form("root_output/analysis_%s_gamma_cask%d.root",   pathTag.c_str(), ic);
        const std::string fnameN =
            Form("root_output/analysis_%s_neutron_cask%d.root", pathTag.c_str(), ic);

        fs::file_time_type inNewest;
        const long nNow = CountInputs(caskDir, inNewest);
        const long nExp = (long)kNAssem * (long)Sources().size();

        bool rebuild = false;

        if (!HaveGoodFile(fnameG) || !HaveGoodFile(fnameN)) {
            rebuild = true;
            std::cerr << "[rebuild] cask " << ic
                      << ": cache missing / zombie / recovered\n";
        } else {
            const long nCached = CachedInputCount(fnameG);
            const auto outMt   = fs::last_write_time(fnameG);

            if (nCached < nNow) {
                rebuild = true;
                std::cerr << "[rebuild] cask " << ic << ": inputs grew "
                          << nCached << " -> " << nNow
                          << " since cache was built\n";
            } else if (inNewest > outMt) {
                rebuild = true;
                std::cerr << "[rebuild] cask " << ic
                          << ": an input file is newer than the cache\n";
            }
        }

        if (LooksLikeMultidetPath(pathTag)) {
            const char* mdObj = "h2_detector_phiz_counts_total";
            if (!CacheHasObject(fnameG, mdObj) || !CacheHasObject(fnameN, mdObj)) {
                rebuild = true;
                std::cerr << "[rebuild] cask " << ic
                          << ": multidet cache lacks " << mdObj << "\n";
            }
        }

        if (nNow < nExp) {
            std::cerr << "\033[1;31m[WARN] cask " << ic << ": only "
                      << nNow << "/" << nExp
                      << " input files present (" << (nExp - nNow)
                      << " missing) -- results reflect a PARTIAL dataset\033[0m\n";
        } else {
            std::cout << "[ok] cask " << ic << ": all "
                      << nExp << " inputs present\n";
        }

        if (rebuild) {
            fs::create_directories(fs::path(fnameG).parent_path());

            std::cout << "[plot_heatmap] rebuilding cask "
                      << (ic + 1) << "/" << ncask
                      << " from " << caskDir
                      << " using " << (slimAnalysis ? "slim" : "full")
                      << " analysis\n";

            analyze_cask(caskDir.c_str(),
                         { {"gamma", fnameG}, {"neutron", fnameN} },
                         typeDet, level,
                         kCaskPos[ic].x, kCaskPos[ic].y,
                         false,
                         slimAnalysis);
        } else {
            std::cout << "[plot_heatmap] using existing cache for cask "
                      << (ic + 1) << "/" << ncask << "\n";
        }

        const std::string fname =
            Form("root_output/analysis_%s_%s_cask%d.root",
                 pathTag.c_str(), typePart.c_str(), ic);

        const Int_t prev = gErrorIgnoreLevel;
        gErrorIgnoreLevel = kFatal;
        auto* f = TFile::Open(fname.c_str());
        gErrorIgnoreLevel = prev;

        if (!f || f->IsZombie()) {
            std::cerr << "[err] could not obtain " << fname << std::endl;
            if (f) {
                f->Close();
                delete f;
            }
            continue;
        }

        auto* raw = (TH2Poly*)f->Get("h2_counts_total");
        if (raw) {
            maps[ic] = RotateTH2Poly(raw, kRotDeg,
                          Form("h2_counts_total_rot_%s_cask%d", dir.c_str(), ic));
        }

        auto* s = (TH1*)f->Get("h1_edep_sum_all_assemblies");
        if (s) {
            s = (TH1*)s->Clone(Form("spec_%s_cask%d", dir.c_str(), ic));
            s->SetDirectory(nullptr);
            spec[ic] = s;

            if (!specSum) {
                specSum = (TH1*)s->Clone(Form("h1_spectra_sum_%s", dir.c_str()));
                specSum->SetDirectory(nullptr);
            } else {
                specSum->Add(s);
            }
        }

        f->Close();
        delete f;
    }

    std::vector<TH2*> arr(maps.begin(), maps.end());
    ScaleHistogramsToGlobalMax(arr.data(), (int)arr.size());

    auto* c = new TCanvas(Form("c_grid_%s_%s", dir.c_str(), typePart.c_str()),
                          Form("%s, %s", dir.c_str(), typePart.c_str()),
                          1900, 680);
    c->Divide(kPadCols, kPadRows, 0.001, 0.001);

    for (int ic = 0; ic < ncask; ++ic) {
        if (!maps[ic]) continue;

        c->cd(kCaskPad[ic]);
        gPad->SetRightMargin(0.14);
        gPad->SetLeftMargin(0.10);
        maps[ic]->SetTitle(Form("cask %d", ic));
        maps[ic]->Draw("COLZ TEXT L");
    }

    auto* t = new TLatex();
    t->SetNDC(true);
    t->SetTextSize(0.03);
    c->cd(0);
    t->DrawLatex(0.35, 0.975, Form("%s  (%s)", dir.c_str(), typePart.c_str()));

    auto* c2 = new TCanvas(Form("c_spec_%s_%s", dir.c_str(), typePart.c_str()),
                           "spectra by cask", 1100, 750);
    c2->SetLogy();
    c2->SetGrid();

    auto* leg = new TLegend(0.78, 0.45, 0.97, 0.90);
    leg->SetTextSize(0.022);

    if (specSum) {
        specSum->SetLineColor(kBlack);
        specSum->SetLineStyle(2);
        specSum->SetLineWidth(1);
        specSum->Draw("HIST");
        leg->AddEntry(specSum, "sum (all casks)", "l");
    }

    for (int ic = 0; ic < ncask; ++ic) {
        if (!spec[ic]) continue;

        spec[ic]->SetLineColor(CaskColor(ic));
        spec[ic]->SetLineWidth(CaskLineWidth(ic));
        spec[ic]->Draw("HIST SAME");
        leg->AddEntry(spec[ic],
                      Form("cask %d%s", ic, ic == 0 ? " (target)" : ""),
                      "l");
    }

    if (specSum) specSum->Draw("HIST SAME");
    leg->Draw();

    if (LooksLikeMultidetPath(pathTag)) {
        auto R = BuildDetectorPhiZSummary(pathTag, typePart, dir);
        DrawDetectorPhiZSummary(R, Form("%s_%s", dir.c_str(), typePart.c_str()));
    }

    TH1* out = nullptr;
    if (specSum) {
        out = (TH1*)specSum->Clone(Form("h1_spectra_sum_overlay_%s", dir.c_str()));
        out->SetDirectory(nullptr);
    }

    gSystem->ProcessEvents();

    return out;
}

// Cross-config comparison.
void plot_design_comp(std::string part = "gamma",
                      DetectorType det = kCLYC,
                      int level = 2,
                      bool slimAnalysis = false)
{
    const auto& configs = Configs();
    const int nconfig = (int)configs.size();
    const int ncask   = kNCask;

    std::cout << "[plot_design_comp] particle=" << part
              << " mode=" << (slimAnalysis ? "slim" : "full")
              << "\n";

    TH1* hcounts[kNCask];
    TH1* hnorm[kNCask];

    auto* stack = new THStack();
    auto* stackNorm = new THStack();

    for (int ic = 0; ic < ncask; ++ic) {
        hcounts[ic] = new TH1F(Form("hcounts_%d", ic), "", nconfig, 0, nconfig);
        hcounts[ic]->SetFillColor(CaskColor(ic));
        stack->Add(hcounts[ic]);

        hnorm[ic] = new TH1F(Form("hnorm_%d", ic), "", nconfig, 0, nconfig);
        hnorm[ic]->SetFillColor(CaskColor(ic));
        stackNorm->Add(hnorm[ic]);
    }

    std::vector<TH1*> configSum(nconfig, nullptr);

    for (int j = 0; j < nconfig; ++j) {
        std::cout << "[plot_design_comp] config "
                  << (j + 1) << "/" << nconfig
                  << ": " << configs[j].dir << "\n";

        configSum[j] = plot_heatmap("data/" + configs[j].dir,
                                    part, det, level, slimAnalysis);

        for (int ic = 0; ic < ncask; ++ic) {
            std::string fname = Form("root_output/analysis_%s_%s_cask%d.root",
                                     configs[j].dir.c_str(), part.c_str(), ic);

            auto* f = TFile::Open(fname.c_str());
            if (!f || f->IsZombie()) {
                if (f) {
                    f->Close();
                    delete f;
                }
                continue;
            }

            auto* htmp = (TH2Poly*)f->Get("h2_counts_total");

            double integral = 0.0;
            double err2 = 0.0;

            if (htmp) {
                for (int b = 1; b <= htmp->GetNumberOfBins(); ++b) {
                    integral += htmp->GetBinContent(b);
                    double be = htmp->GetBinError(b);
                    err2 += be * be;
                }
            }

            hcounts[ic]->SetBinContent(j + 1, integral);
            hcounts[ic]->SetBinError  (j + 1, std::sqrt(err2));
            hcounts[ic]->GetXaxis()->SetBinLabel(j + 1, configs[j].title.c_str());
            hnorm[ic]  ->GetXaxis()->SetBinLabel(j + 1, configs[j].title.c_str());

            f->Close();
            delete f;
        }

        double colT = 0.0;
        double colE2 = 0.0;

        for (int ic = 0; ic < ncask; ++ic) {
            colT += hcounts[ic]->GetBinContent(j + 1);
            double e = hcounts[ic]->GetBinError(j + 1);
            colE2 += e * e;
        }

        if (colT > 0) {
            for (int ic = 0; ic < ncask; ++ic) {
                double a = hcounts[ic]->GetBinContent(j + 1);
                double sa2 = hcounts[ic]->GetBinError(j + 1);
                sa2 *= sa2;

                double fr = a / colT;
                double restE2 = colE2 - sa2;
                double var = ((1 - fr) * (1 - fr) * sa2 + fr * fr * restE2)
                           / (colT * colT);

                hnorm[ic]->SetBinContent(j + 1, 100.0 * fr);
                hnorm[ic]->SetBinError  (j + 1, 100.0 * std::sqrt(var));
            }
        }
    }

    auto drawLabels = [&](THStack*, TH1** h, const char* fmt) {
        TLatex tex;
        tex.SetTextAlign(22);
        tex.SetTextFont(42);
        tex.SetTextSize(0.018);
        tex.SetTextColor(kWhite);

        for (int j = 0; j < nconfig; ++j) {
            double cum = 0.0;
            double x = h[0]->GetXaxis()->GetBinCenter(j + 1);

            for (int ic = 0; ic < ncask; ++ic) {
                double v = h[ic]->GetBinContent(j + 1);
                double e = h[ic]->GetBinError(j + 1);
                double y = cum + v / 2.0;
                cum += v;

                if (v > 0) tex.DrawLatex(x, y, Form(fmt, v, e));
            }
        }
    };

    auto makeLegend = [&](TH1** h) {
        auto* leg = new TLegend(0.80, 0.35, 0.93, 0.92);
        leg->SetHeader("Cask");
        leg->SetTextSize(0.024);

        for (int ic = 0; ic < ncask; ++ic) {
            leg->AddEntry(h[ic],
                          Form("Cask %d%s", ic, ic == 0 ? " (target)" : ""),
                          "f");
        }

        return leg;
    };

    auto* c1 = new TCanvas("c_counts", "raw counts", 1000, 720);
    c1->SetBottomMargin(0.25);
    stack->Draw("hist");
    stack->GetXaxis()->LabelsOption("v");
    stack->GetXaxis()->SetLabelSize(0.03);
    gPad->Update();
    drawLabels(stack, hcounts, "%.0f#pm%.0f");
    makeLegend(hcounts)->Draw();
    gPad->Modified();
    gPad->Update();

    auto* c2 = new TCanvas("c_norm", "normalized", 1000, 720);
    c2->SetBottomMargin(0.25);
    stackNorm->SetMaximum(110.0);
    stackNorm->Draw("hist");
    stackNorm->GetXaxis()->LabelsOption("v");
    stackNorm->GetXaxis()->SetLabelSize(0.03);
    stackNorm->GetYaxis()->SetTitle("Fraction of total per config [%]");
    gPad->Update();
    drawLabels(stackNorm, hnorm, "%.1f#pm%.1f");
    makeLegend(hnorm)->Draw();
    gPad->Modified();
    gPad->Update();

    auto* c3 = new TCanvas("c_specsum_allconfig",
                           "summed spectra, all configs", 1100, 750);
    c3->SetLogy();
    c3->SetGrid();

    auto* leg3 = new TLegend(0.50, 0.75, 0.96, 0.90);
    leg3->SetTextSize(0.022);

    bool first = true;

    for (int j = 0; j < nconfig; ++j) {
        if (!configSum[j]) continue;

        configSum[j]->SetLineColor(j == 0 ? kBlack : getP6((j - 1) % 6));
        configSum[j]->SetLineStyle(j == 0 ? 1 : 1 + (j - 1) / 6);
        configSum[j]->SetLineWidth(j == 0 ? 2 : 1);
        configSum[j]->SetTitle("Summed cask spectrum by configuration;"
                               "edep [MeV];counts");

        configSum[j]->Draw(first ? "HIST" : "HIST SAME");
        first = false;

        leg3->AddEntry(configSum[j], configs[j].title.c_str(), "l");
    }

    leg3->Draw();
    gPad->Modified();
    gPad->Update();
}

#endif

