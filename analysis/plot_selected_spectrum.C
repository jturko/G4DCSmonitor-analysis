#ifndef PLOT_SELECTED_SPECTRUM_C
#define PLOT_SELECTED_SPECTRUM_C 1

// plot_selected_spectrum.C
//   General recreator of plot_design_comp's overhead 2x6 per-assembly TH2Poly
//   grid, restricted by ANY combination of {cask, detector, globalFuel, source}.
//   Any selector = -1 (or "all" for source) means "sum over that axis".
//   Every one of the 12 cask pads is always drawn (kCaskPad layout); pads/bins
//   outside the selection are simply left at zero. The matching edep spectrum
//   for the same selection is produced alongside. Everything is saved to a
//   ROOT file.
//
//   Reuses: DiscoverDatasetFiles / OpenHitsTree / FileSlot
//           (compare_nominal_vs_multidet_spectrum.C),
//           MakeAssemblyHexMap / BuildFuelPositions (dcs_geometry.h),
//           kCaskPad / kPadCols / kPadRows (dcs_setup.h),
//           ScaleHistogramsToGlobalMax (utils.h),
//           PassSpectrum / PassSignal / EnergyWindow (dcs_cuts.h).
//
//   root -l 'plot_selected_spectrum.C()'
//     -> all casks, all detectors, all assemblies, all sources
//        (reproduces the plot_design_comp c_grid_..._gamma grid)
//   root -l 'plot_selected_spectrum.C("data/nominal-all-positions","gamma",0,96)'
//     -> cask 0, detector 96, all assemblies: only the cask-0 pad populated
//   root -l 'plot_selected_spectrum.C("data/nominal-all-positions","gamma",1,-1,61)'
//     -> cask 1, all detectors, globalFuel 61: only bin 61 of cask-1 pad lit
//   root -l 'plot_selected_spectrum.C("data/nominal-all-positions","gamma",-1,96,-1,"gamma_662keV")'
//     -> all casks, det 96, 662 keV source only

#include "compare_nominal_vs_multidet_spectrum.C"   // detail ns + core headers
#include "dcs_geometry.h"                           // MakeAssemblyHexMap, BuildFuelPositions
#include "dcs_multidet.h"
#include "utils.h"                                   // ScaleHistogramsToGlobalMax

#include <TCanvas.h>
#include <TH2Poly.h>
#include <TLatex.h>
#include <TLine.h>

#include <vector>

void plot_selected_spectrum(int          caskSelect = -1,   // -1 = all casks
                            int          detSelect  = -1,   // -1 = all detectors
                            int          globalFuel = -1,   // -1 = all assemblies
                            std::string  dataDir    = "data/nominal-all-positions",
                            std::string  part       = "gamma",
                            std::string  sourceTag  = "all",// "all" = all sources
                            DetectorType det        = kCLYC,
                            int          level      = 2,
                            int          rebin      = 10,
                            bool         logY       = true,
                            std::string  outFile    = "")
{
    using namespace compare_nominal_md_detail;   // DiscoverDatasetFiles, OpenHitsTree, FileSlot
    using namespace dcs;                          // MakeAssemblyHexMap, kCaskPad, ...

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(255);
    if (sourceTag == "all") sourceTag.clear();

    CutConfig cfg{ det, ParseParticle(part), (CutLevel)level };
    double eLo = 0, eHi = 0;
    EnergyWindow(cfg, eLo, eHi);

    const std::string caskLbl = (caskSelect < 0) ? "allcask" : Form("cask%d",  caskSelect);
    const std::string detLbl  = (detSelect  < 0) ? "alldet"  : Form("det%03d", detSelect);
    const std::string fuelLbl  = (globalFuel < 0) ? "allfuel" : Form("gf%02d",   globalFuel);
    const std::string srcLbl  = sourceTag.empty() ? "allsrc" : sourceTag;

    if (outFile.empty())
        outFile = Form("root_output/spectrum_%s_%s_%s_%s_%s_%s.root",
                       fs::path(dataDir).filename().string().c_str(),
                       part.c_str(), caskLbl.c_str(), detLbl.c_str(),
                       fuelLbl.c_str(), srcLbl.c_str());
    fs::create_directories(fs::path(outFile).parent_path());

    std::cout << "\n[plot_selected_spectrum]"
              << " cask="   << (caskSelect < 0 ? "all" : std::to_string(caskSelect))
              << " det="    << (detSelect  < 0 ? "all" : std::to_string(detSelect))
              << " gFuel="  << (globalFuel < 0 ? "all" : std::to_string(globalFuel))
              << " source=" << srcLbl
              << " | spectrum=PassSpectrum ; grid=PassSignal ROI ["
              << eLo << "," << eHi << "] MeV\n";

    // Discover files for ALL casks; the {cask,fuel,source} cuts are applied
    // here, the detector cut is applied per-hit.
    auto files = DiscoverDatasetFiles(dataDir, -1, sourceTag);

    // ---- spectrum (PassSpectrum, PID only) for the SELECTED slice ----
    auto* h = new TH1D(
        Form("h_edep_%s_%s_%s_%s",
             caskLbl.c_str(), detLbl.c_str(), fuelLbl.c_str(), srcLbl.c_str()),
        Form("edep  (%s, %s, %s, %s, %s);E_{dep} [MeV];counts / bin",
             part.c_str(), caskLbl.c_str(), detLbl.c_str(),
             fuelLbl.c_str(), srcLbl.c_str()),
        kNbE, kEmin, kEmax);
    h->Sumw2();
    h->SetDirectory(nullptr);

    // ---- 12 overhead per-assembly maps (PassSignal, ROI), one per cask ----
    const auto positions = BuildFuelPositions();
    std::vector<TH2Poly*> maps(kNCask, nullptr);
    std::vector<std::vector<double>> assemCounts(
        kNCask, std::vector<double>(kNAssem, 0.0));   // [cask][globalFuel] ROI signal

    for (int ic = 0; ic < kNCask; ++ic) {
        maps[ic] = MakeAssemblyHexMap(
            Form("h2_counts_total_%s_%s_%s_%s_cask%d",
                 caskLbl.c_str(), detLbl.c_str(), fuelLbl.c_str(),
                 srcLbl.c_str(), ic),
            Form("cask %d", ic),
            positions);
    }

    long long nFiles = 0, nHits = 0, nFilled = 0, nSigTot = 0, nNoDet = 0;

    for (const auto& s : files) {
        if (s.cask < 0 || s.cask >= kNCask)                 continue;   // safety
        if (caskSelect >= 0 && s.cask != caskSelect)        continue;   // cask cut
        if (globalFuel >= 0 && s.globalFuel != globalFuel)  continue;   // fuel cut
        if (s.globalFuel < 0 || s.globalFuel >= kNAssem)    continue;   // safety
        // (source cut already applied by DiscoverDatasetFiles via sourceTag)

        TFile* f = nullptr; TTree* t = nullptr;
        if (!OpenHitsTree(s.path, f, t)) continue;
        ++nFiles;

        if (!t->GetBranch("edep") || !t->GetBranch("pid")) { f->Close(); delete f; continue; }

        double h_edep = 0, h_pid = 0; Int_t h_det = -1;
        t->SetBranchStatus("*", 0);
        t->SetBranchStatus("edep", 1);
        t->SetBranchStatus("pid",  1);
        const bool haveDet = (t->GetBranch("det") != nullptr);
        if (haveDet) t->SetBranchStatus("det", 1); else ++nNoDet;
        t->SetBranchAddress("edep", &h_edep);
        t->SetBranchAddress("pid",  &h_pid);
        if (haveDet) t->SetBranchAddress("det", &h_det);

        long long nSigThisFile = 0;
        const Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            ++nHits;

            // detector selection applies to BOTH products
            if (detSelect >= 0) { if (!haveDet || h_det != detSelect) continue; }

            const long long pid = std::llround(h_pid);

            // spectrum: PID-only (matches analyze_cask's spectrum filling)
            if (PassSpectrum(cfg, pid)) { h->Fill(h_edep); ++nFilled; }

            // grid: PID + ROI window (matches analyze_cask's h2_counts_total)
            if (PassSignal(cfg, pid, h_edep)) ++nSigThisFile;
        }

        assemCounts[s.cask][s.globalFuel] += (double)nSigThisFile;
        nSigTot                           += nSigThisFile;

        f->Close(); delete f;
    }

    // push per-assembly ROI counts into each cask's TH2Poly (bin g+1 <-> gFuel g)
    for (int ic = 0; ic < kNCask; ++ic)
        for (int g = 0; g < kNAssem; ++g)
            if (assemCounts[ic][g] != 0.0)
                maps[ic]->SetBinContent(g + 1, assemCounts[ic][g]);

    if ((detSelect >= 0) && nNoDet > 0)
        std::cerr << "[warn] " << nNoDet
                  << " files lacked hits.det while det filtering was on\n";

    const double roiInt = h->Integral(h->GetXaxis()->FindBin(eLo + 1e-9),
                                      h->GetXaxis()->FindBin(eHi - 1e-9));
    std::cout << "[plot_selected_spectrum] files=" << nFiles
              << " hits=" << nHits << " spectrum-filled=" << nFilled
              << " | ROI(spectrum PID+window) integral=" << roiInt
              << " | grid ROI signal total=" << nSigTot << "\n";

    if (rebin > 1) h->Rebin(rebin);

    // ---- spectrum canvas ----
    auto* cSpec = new TCanvas(
        Form("c_spec_%s_%s_%s_%s",
             caskLbl.c_str(), detLbl.c_str(), fuelLbl.c_str(), srcLbl.c_str()),
        "selected edep spectrum", 1100, 720);
    if (logY) cSpec->SetLogy();
    cSpec->SetGrid();
    h->SetLineColor(getP6(0));
    h->SetLineWidth(2);
    h->SetMinimum(logY ? 0.5 : 0.0);
    h->Draw("HIST");
    if (cfg.level == kPidEnergy) {
        const double y1 = (logY ? 20.0 : 1.35) * std::max(1.0, h->GetMaximum());
        for (double x : {eLo, eHi}) {
            auto* l = new TLine(x, logY ? 0.5 : 0.0, x, y1);
            l->SetLineColor(kGray + 2); l->SetLineStyle(2); l->Draw();
        }
    }

    // ---- overhead 2x6 assembly grid (same layout/scale as plot_heatmap) ----
    std::vector<TH2*> arr(maps.begin(), maps.end());
    ScaleHistogramsToGlobalMax(arr.data(), (int)arr.size());   // common z-range

    auto* cGrid = new TCanvas(
        Form("c_grid_%s_%s_%s_%s",
             caskLbl.c_str(), detLbl.c_str(), fuelLbl.c_str(), srcLbl.c_str()),
        Form("assembly ROI signal grid (%s, %s, %s, %s, %s)",
             part.c_str(), caskLbl.c_str(), detLbl.c_str(),
             fuelLbl.c_str(), srcLbl.c_str()),
        1900, 680);
    cGrid->Divide(kPadCols, kPadRows, 0.001, 0.001);

    for (int ic = 0; ic < kNCask; ++ic) {
        if (!maps[ic]) continue;
        cGrid->cd(kCaskPad[ic]);
        gPad->SetRightMargin(0.14);
        gPad->SetLeftMargin(0.10);
        maps[ic]->SetTitle(Form("cask %d%s", ic, ic == 0 ? " (target)" : ""));
        maps[ic]->Draw("COLZ TEXT L");
    }

    cGrid->cd(0);
    {
        TLatex tx; tx.SetNDC(true); tx.SetTextSize(0.03);
        tx.DrawLatex(0.22, 0.975,
            Form("%s | %s | %s | %s   (ROI [%.3g,%.3g] MeV, total=%.4g cts)",
                 caskLbl.c_str(), detLbl.c_str(), fuelLbl.c_str(), srcLbl.c_str(),
                 eLo, eHi, (double)nSigTot));
    }

    // ---- write everything ----
    TFile* fout = TFile::Open(outFile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[err] cannot open output " << outFile << "\n";
        return;
    }
    h->Write();
    for (int ic = 0; ic < kNCask; ++ic) if (maps[ic]) maps[ic]->Write();
    cSpec->Write();
    cGrid->Write();
    TParameter<Int_t>("cask_select", caskSelect).Write();
    TParameter<Int_t>("det_select",  detSelect).Write();
    TParameter<Int_t>("global_fuel", globalFuel).Write();
    TNamed("source_select", srcLbl.c_str()).Write();
    TParameter<Long64_t>("grid_roi_signal_total", (Long64_t)nSigTot).Write();
    fout->Close();

    std::cout << "[ok] wrote " << outFile << "\n";
}

#endif

