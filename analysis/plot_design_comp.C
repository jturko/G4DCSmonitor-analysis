#ifndef PLOT_DESIGN_COMP
#define PLOT_DESIGN_COMP 1

#include "analyze_cask.C"
#include "dcs_geometry.h"
#include "dcs_cuts.h"
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

// New hall layout is axis-aligned -> no per-map rotation (old diamond used -45).
namespace { constexpr double kRotDeg = 0.0; }

// Count currently-present (globalFuel, source) inputs under caskDir and the
// newest input mtime. Reuses analyze_cask's own discovery so the count matches.
static long CountInputs(const std::string& caskDir, fs::file_time_type& newest)
{
    std::vector<std::vector<std::string>> fp;
    DiscoverFiles(caskDir.c_str(), Sources(), fp);      // fills [kNAssem][nsrc]
    long n = 0; newest = fs::file_time_type::min();
    for (auto& row : fp) for (auto& s : row) if (!s.empty()) {
        ++n;
        std::error_code ec; auto t = fs::last_write_time(s, ec);
        if (!ec && t > newest) newest = t;
    }
    return n;
}

static long CachedInputCount(const std::string& fn)     // -1 if unreadable/unstamped
{
    const Int_t pv = gErrorIgnoreLevel; gErrorIgnoreLevel = kFatal;
    auto* f = TFile::Open(fn.c_str());
    long v = -1;
    if (f && !f->IsZombie() && !f->TestBit(TFile::kRecovered)) {
        if (auto* p = (TParameter<Long64_t>*)f->Get("n_input_files")) v = (long)p->GetVal();
    }
    if (f) { f->Close(); delete f; }
    gErrorIgnoreLevel = pv;
    return v;
}

// Produce/refresh the per-cask analysis files and return the grand-sum spectrum
// (all casks) for one config, as an INDEPENDENT clone.
TH1* plot_heatmap(fs::path p = "data/roomReturn/hall_3cluster/2cm_Pb",
                  std::string typePart = "gamma",
                  DetectorType typeDet = kCLYC,
                  int level = 2)
{
    std::string pathTag = p.lexically_relative("data").generic_string();
    if (pathTag.empty() || pathTag == ".") pathTag = p.filename().string();
    std::string dir = pathTag;
    std::replace(dir.begin(), dir.end(), '/', '_');   // slash-free ROOT names

    const int ncask = kNCask;

    std::vector<TH2Poly*> maps(ncask, nullptr);   // (un)rotated total maps
    std::vector<TH1*>     spec(ncask, nullptr);
    TH1* specSum = nullptr;                        // owned by canvas c2 (black)

    for (int ic = 0; ic < ncask; ++ic) {
        //
        const std::string caskDir = Form("%s/cask%d", p.string().c_str(), ic);
        const std::string fnameG =
            Form("root_output/analysis_%s_gamma_cask%d.root",   pathTag.c_str(), ic);
        const std::string fnameN =
            Form("root_output/analysis_%s_neutron_cask%d.root", pathTag.c_str(), ic);
        
        // reject absent / zombie / *recovered* (in-progress) caches
        auto haveFile = [](const std::string& fn) {
            const Int_t pv = gErrorIgnoreLevel; gErrorIgnoreLevel = kFatal;
            auto* t = TFile::Open(fn.c_str());
            const bool ok = t && !t->IsZombie() && !t->TestBit(TFile::kRecovered);
            if (t) { t->Close(); delete t; }
            gErrorIgnoreLevel = pv;
            return ok;
        };
        
        fs::file_time_type inNewest;
        const long nNow = CountInputs(caskDir, inNewest);
        const long nExp = (long)kNAssem * (long)Sources().size();   // 84 * 5 = 420
        
        bool rebuild = false;
        if (!haveFile(fnameG) || !haveFile(fnameN)) {
            rebuild = true;
            std::cerr << "[rebuild] cask " << ic << ": cache missing / zombie / recovered\n";
        } else {
            const long nCached = CachedInputCount(fnameG);
            const auto outMt    = fs::last_write_time(fnameG);
            if (nCached < nNow) {
                rebuild = true;
                std::cerr << "[rebuild] cask " << ic << ": inputs grew "
                          << nCached << " -> " << nNow << " since cache was built\n";
            } else if (inNewest > outMt) {
                rebuild = true;
                std::cerr << "[rebuild] cask " << ic << ": an input file is newer than the cache\n";
            }
        }
        
        // LOUD completeness report, every run, rebuild or not
        if (nNow < nExp) {
            std::cerr << "\033[1;31m[WARN] cask " << ic << ": only " << nNow << "/" << nExp
                      << " input files present (" << (nExp - nNow)
                      << " missing) -- results reflect a PARTIAL dataset\033[0m\n";
        } else {
            std::cout << "[ok] cask " << ic << ": all " << nExp << " inputs present\n";
        }
        
        if (rebuild) {
            fs::create_directories(fs::path(fnameG).parent_path());
            analyze_cask(caskDir.c_str(),
                         { {"gamma", fnameG}, {"neutron", fnameN} },
                         typeDet, level, kCaskPos[ic].x, kCaskPos[ic].y, false);
        }
        //


        //// One analyze_cask pass produces BOTH particle files; regenerate only
        //// if a file is missing, then display only typePart.
        //const std::string fnameG =
        //    Form("root_output/analysis_%s_gamma_cask%d.root",   pathTag.c_str(), ic);
        //const std::string fnameN =
        //    Form("root_output/analysis_%s_neutron_cask%d.root", pathTag.c_str(), ic);

        //auto haveFile = [](const std::string& fn) {
        //    const Int_t pv = gErrorIgnoreLevel; gErrorIgnoreLevel = kFatal;
        //    auto* t = TFile::Open(fn.c_str());
        //    const bool ok = t && !t->IsZombie();
        //    if (t) { t->Close(); delete t; }
        //    gErrorIgnoreLevel = pv;
        //    return ok;
        //};
        //if (!haveFile(fnameG) || !haveFile(fnameN)) {
        //    fs::create_directories(fs::path(fnameG).parent_path());
        //    analyze_cask(Form("%s/cask%d", p.string().c_str(), ic),
        //                 { {"gamma", fnameG}, {"neutron", fnameN} },
        //                 typeDet, level, kCaskPos[ic].x, kCaskPos[ic].y, false);
        //}

        const std::string fname =
            Form("root_output/analysis_%s_%s_cask%d.root",
                 pathTag.c_str(), typePart.c_str(), ic);
        const Int_t prev = gErrorIgnoreLevel; gErrorIgnoreLevel = kFatal;
        auto* f = TFile::Open(fname.c_str());
        gErrorIgnoreLevel = prev;
        if (!f || f->IsZombie()) {
            std::cerr << "[err] could not obtain " << fname << std::endl;
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
            } else specSum->Add(s);
        }
        f->Close();
    }

    // common colour scale across all cask maps
    std::vector<TH2*> arr(maps.begin(), maps.end());
    ScaleHistogramsToGlobalMax(arr.data(), (int)arr.size());

    // spatially-faithful 6x2 hall grid (columns increase with x; top row = +y)
    auto* c = new TCanvas(Form("c_grid_%s_%s", dir.c_str(), typePart.c_str()),
                          Form("%s, %s", dir.c_str(), typePart.c_str()), 1900, 680);
    c->Divide(kPadCols, kPadRows, 0.001, 0.001);
    for (int ic = 0; ic < ncask; ++ic) {
        if (!maps[ic]) continue;
        c->cd(kCaskPad[ic]);
        gPad->SetRightMargin(0.14); gPad->SetLeftMargin(0.10);
        maps[ic]->SetTitle(Form("cask %d", ic));
        maps[ic]->Draw("COLZ TEXT L");
    }
    auto* t = new TLatex(); t->SetNDC(true); t->SetTextSize(0.03);
    c->cd(0);
    t->DrawLatex(0.35, 0.975, Form("%s  (%s)", dir.c_str(), typePart.c_str()));

    // per-cask spectrum overlay (this canvas OWNS 'specSum' and keeps it black)
    auto* c2 = new TCanvas(Form("c_spec_%s_%s", dir.c_str(), typePart.c_str()),
                           "spectra by cask", 1100, 750);
    c2->SetLogy(); c2->SetGrid();
    auto* leg = new TLegend(0.78, 0.45, 0.97, 0.90); leg->SetTextSize(0.022);
    if (specSum) {
        specSum->SetLineColor(kBlack); specSum->SetLineStyle(2);
        specSum->SetLineWidth(1);      specSum->Draw("HIST");
        leg->AddEntry(specSum, "sum (all casks)", "l");
    }
    for (int ic = 0; ic < ncask; ++ic) if (spec[ic]) {
        spec[ic]->SetLineColor(CaskColor(ic));
        spec[ic]->SetLineWidth(CaskLineWidth(ic));
        spec[ic]->Draw("HIST SAME");
        leg->AddEntry(spec[ic], Form("cask %d%s", ic, ic == 0 ? " (signal)" : ""), "l");
    }
    if (specSum) specSum->Draw("HIST SAME");
    leg->Draw();

    TH1* out = nullptr;
    if (specSum) {
        out = (TH1*)specSum->Clone(Form("h1_spectra_sum_overlay_%s", dir.c_str()));
        out->SetDirectory(nullptr);
    }

    gSystem->ProcessEvents();

    return out;
}

// Cross-config comparison (stacked per-cask counts + per-config spectra).
void plot_design_comp(std::string part = "gamma",
                      DetectorType det = kCLYC, int level = 2)
{
    const auto& configs = Configs();
    const int nconfig = (int)configs.size(), ncask = kNCask;

    TH1* hcounts[kNCask]; TH1* hnorm[kNCask];
    auto* stack = new THStack(); auto* stackNorm = new THStack();
    for (int ic = 0; ic < ncask; ++ic) {
        hcounts[ic] = new TH1F(Form("hcounts_%d", ic), "", nconfig, 0, nconfig);
        hcounts[ic]->SetFillColor(CaskColor(ic)); stack->Add(hcounts[ic]);
        hnorm[ic]   = new TH1F(Form("hnorm_%d", ic),   "", nconfig, 0, nconfig);
        hnorm[ic]->SetFillColor(CaskColor(ic));   stackNorm->Add(hnorm[ic]);
    }

    std::vector<TH1*> configSum(nconfig, nullptr);

    for (int j = 0; j < nconfig; ++j) {
        configSum[j] = plot_heatmap("data/" + configs[j].dir, part, det, level);

        for (int ic = 0; ic < ncask; ++ic) {
            std::string fname = Form("root_output/analysis_%s_%s_cask%d.root",
                                     configs[j].dir.c_str(), part.c_str(), ic);
            auto* f = TFile::Open(fname.c_str());
            if (!f || f->IsZombie()) continue;
            auto* htmp = (TH2Poly*)f->Get("h2_counts_total");
            double integral = 0, err2 = 0;
            if (htmp) for (int b = 1; b <= htmp->GetNumberOfBins(); ++b) {
                integral += htmp->GetBinContent(b);
                double be = htmp->GetBinError(b); err2 += be * be;
            }
            hcounts[ic]->SetBinContent(j + 1, integral);
            hcounts[ic]->SetBinError  (j + 1, std::sqrt(err2));
            hcounts[ic]->GetXaxis()->SetBinLabel(j + 1, configs[j].title.c_str());
            hnorm[ic]  ->GetXaxis()->SetBinLabel(j + 1, configs[j].title.c_str());
            f->Close();
        }
        double colT = 0, colE2 = 0;
        for (int ic = 0; ic < ncask; ++ic) {
            colT += hcounts[ic]->GetBinContent(j + 1);
            double e = hcounts[ic]->GetBinError(j + 1); colE2 += e * e;
        }
        if (colT > 0) for (int ic = 0; ic < ncask; ++ic) {
            double a = hcounts[ic]->GetBinContent(j + 1);
            double sa2 = hcounts[ic]->GetBinError(j + 1); sa2 *= sa2;
            double fr = a / colT, restE2 = colE2 - sa2;
            double var = ((1 - fr) * (1 - fr) * sa2 + fr * fr * restE2) / (colT * colT);
            hnorm[ic]->SetBinContent(j + 1, 100.0 * fr);
            hnorm[ic]->SetBinError  (j + 1, 100.0 * std::sqrt(var));
        }
    }

    auto drawLabels = [&](THStack*, TH1** h, const char* fmt) {
        TLatex tex; tex.SetTextAlign(22); tex.SetTextFont(42);
        tex.SetTextSize(0.018); tex.SetTextColor(kWhite);
        for (int j = 0; j < nconfig; ++j) {
            double cum = 0, x = h[0]->GetXaxis()->GetBinCenter(j + 1);
            for (int ic = 0; ic < ncask; ++ic) {
                double v = h[ic]->GetBinContent(j + 1), e = h[ic]->GetBinError(j + 1);
                double y = cum + v / 2.0; cum += v;
                if (v > 0) tex.DrawLatex(x, y, Form(fmt, v, e));
            }
        }
    };
    auto makeLegend = [&](TH1** h) {
        auto* leg = new TLegend(0.80, 0.35, 0.93, 0.92);
        leg->SetHeader("Cask"); leg->SetTextSize(0.024);
        for (int ic = 0; ic < ncask; ++ic)
            leg->AddEntry(h[ic], Form("Cask %d%s", ic, ic == 0 ? " (sig)" : ""), "f");
        return leg;
    };

    auto* c1 = new TCanvas("c_counts", "raw counts", 1000, 720);
    c1->SetBottomMargin(0.25);
    stack->Draw("hist");
    stack->GetXaxis()->LabelsOption("v"); stack->GetXaxis()->SetLabelSize(0.03);
    gPad->Update(); drawLabels(stack, hcounts, "%.0f#pm%.0f");
    makeLegend(hcounts)->Draw(); gPad->Modified(); gPad->Update();

    auto* c2 = new TCanvas("c_norm", "normalized", 1000, 720);
    c2->SetBottomMargin(0.25); stackNorm->SetMaximum(110.0);
    stackNorm->Draw("hist");
    stackNorm->GetXaxis()->LabelsOption("v"); stackNorm->GetXaxis()->SetLabelSize(0.03);
    stackNorm->GetYaxis()->SetTitle("Fraction of total per config [%]");
    gPad->Update(); drawLabels(stackNorm, hnorm, "%.1f#pm%.1f");
    makeLegend(hnorm)->Draw(); gPad->Modified(); gPad->Update();

    // per-config summed spectrum (nominal j==0 -> black; palette otherwise)
    auto* c3 = new TCanvas("c_specsum_allconfig",
                           "summed spectra, all configs", 1100, 750);
    c3->SetLogy(); c3->SetGrid();
    auto* leg3 = new TLegend(0.50, 0.75, 0.96, 0.90); leg3->SetTextSize(0.022);
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
    gPad->Modified(); gPad->Update();
}
#endif

