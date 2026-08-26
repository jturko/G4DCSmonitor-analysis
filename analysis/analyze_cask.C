// analyze_cask.C
//   root -l -b -q 'analyze_cask.C("data/nominal/cask0","out.root","gamma",kCLYC,2)'
//
// One pass per file over the 'hits' tree (with the 'primary' join available),
// producing EVERYTHING for one cask into a single ROOT file:
//   edep_per_assembly/    h1_edep_globalFuelGG_<src>
//   edep_sum_per_assembly/ h1_edep_sum_globalFuelGG
//   h1_edep_sum_all_assemblies[_<src>]
//   h2_counts_<src>, h2_counts_total            (TH2Poly hex maps)
//   h2_edep_vs_ekin
//   h2_primary_zphi_wall, h2_primary_xy_top, h2_primary_xy_bottom
//   c_assembly_heatmaps, c_edep_sum_all_assemblies, (optional per-assembly)

#include "geometry_constraints.h"
#include "dcs_geometry.h"
#include "dcs_cuts.h"
#include "dcs_sources.h"
#include "dcs_setup.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>   // std::pair
#include <TParameter.h>

#include <filesystem>
#include <regex>
namespace fs = std::filesystem;

using namespace geo;
using namespace dcs;

namespace {

#include <cerrno>
#include <cstring>

static int CountOpenFDs() {
    int n = 0; std::error_code ec;
    for (auto it = fs::directory_iterator("/proc/self/fd", ec);
         !ec && it != fs::directory_iterator(); it.increment(ec)) ++n;
    return n;
}

bool OpenGood(const char* fname, TFile*& fin, TTree*& hits, TTree*& prim) {
    fin = nullptr; hits = nullptr; prim = nullptr;
    errno = 0;
    fin = TFile::Open(fname, "READ");
    const int  e     = errno;
    const bool isnull= (fin == nullptr);
    const bool zomb  = (fin && fin->IsZombie());
    const bool recov = (fin && fin->TestBit(TFile::kRecovered));
    if (fin && !zomb && !recov) hits = (TTree*)fin->Get("hits");
    if (hits) { prim = (TTree*)fin->Get("primary"); return true; }

    std::cerr << "[OpenGood FAIL] " << fname
              << "  null="  << isnull << " zombie=" << zomb
              << " recov="  << recov  << " hasHits=" << (hits!=nullptr)
              << " errno="  << e << "(" << std::strerror(e) << ")"
              << " openFDs="<< CountOpenFDs() << "\n";
    if (fin) { fin->Close(); delete fin; fin = nullptr; }
    hits = nullptr;
    return false;
}




void DrawWindow(double lo, double hi, double ymin, double ymax) {
    for (double x : {lo, hi}) {
        auto* l = new TLine(x, ymin, x, ymax);
        l->SetLineColor(kGray + 2); l->SetLineStyle(2); l->SetLineWidth(2);
        l->Draw();
    }
    auto* t = new TLatex();
    t->SetTextSize(0.028); t->SetTextColor(kGray + 2);
    t->DrawLatex(hi, ymax * 0.4, Form("  gate [%.3g, %.3g] MeV", lo, hi));
}

// Discover detector-response files in 'dir' and map them onto (globalFuel, src).
// Globs + parses instead of constructing exact names, so it is immune to
// naming variants (optional "cask<N>_" prefix, zero-padding of the index).
// Mirrors how compare_edep_ekin.C finds its files [[17]].
void DiscoverFiles(const char* dir,
                   const std::vector<SrcStyle>& sources,
                   std::vector<std::vector<std::string>>& out)
{
    out.assign(kNAssem, std::vector<std::string>(sources.size(), std::string()));
    if (!fs::exists(dir)) {
        std::cerr << "[warn] input dir does not exist: " << dir << std::endl;
        return;
    }
    static const std::regex reFuel("globalFuel0*([0-9]+)");
    for (const auto& e : fs::directory_iterator(dir)) {
        if (!e.is_regular_file()) continue;
        const std::string name = e.path().filename().string();
        if (name.rfind("detector-response", 0) != 0) continue;   // prefix
        if (e.path().extension() != ".root")        continue;

        std::smatch m;
        if (!std::regex_search(name, m, reFuel)) continue;        // globalFuel<idx>
        const int g = std::stoi(m[1].str());
        if (g < 0 || g >= kNAssem) continue;

        int si = -1;                                              // which source tag?
        for (size_t k = 0; k < sources.size(); ++k)
            if (name.find(sources[k].tag) != std::string::npos) { si = (int)k; break; }
        if (si < 0) continue;

        out[g][si] = e.path().string();
    }
}


} // anonymous



// One pass over the raw hits/primary trees, filling an INDEPENDENT set of
// histograms for every requested particle simultaneously, then writing each
// particle's set to its own file. 'outputs' = list of {particle, outFile}.
void analyze_cask(const char* inputDir = "data/nominal/cask0",
                  const std::vector<std::pair<std::string,std::string>>& outputs
                      = { {"gamma",   "analysis_gamma.root"},
                          {"neutron", "analysis_neutron.root"} },
                  DetectorType det     = kCLYC,
                  int level            = 2,
                  double caskX         = 0.0,   // global cask centre [mm]
                  double caskY         = 0.0,
                  bool makePerAssemblyCanvases = false)
{
    gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(256);
    gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);

    const auto& sources    = Sources();
    const size_t nsrc      = sources.size();
    const auto   positions = BuildFuelPositions();
    if ((int)positions.size() != kNAssem) {
        std::cerr << "[err] expected " << kNAssem << " positions, got "
                  << positions.size() << std::endl; return;
    }

    // ---- one fully independent output bundle per requested particle --------
    struct PartOut {
        std::string particle;
        CutConfig   cfg;
        double      eLo = 0, eHi = 0;
        TFile*      fout = nullptr;
        TDirectory *edepDir = nullptr, *sumDir = nullptr;
        std::vector<TH2Poly*> hMaps;
        TH2Poly*    hSum = nullptr;
        TH1D*       hAll = nullptr;
        std::vector<TH1D*> hAllSrc;
        TH2F *hEk=nullptr, *hWall=nullptr, *hTop=nullptr, *hBot=nullptr;
        TH2F *hPhiE=nullptr, *hZE=nullptr;   // energy-RESOLVED wall maps (PID, ungated)
        TH1D *hCos=nullptr;
        std::vector<std::vector<TH1D*>> hcell;
        std::vector<std::vector<bool>>  missing;
        long long nEmpty=0, nMissing=0, nMatched=0, nOrphan=0;
    };

    std::vector<PartOut> po;
    po.reserve(outputs.size());

    for (const auto& pr : outputs) {
        const std::string& particle = pr.first;
        const std::string& outFile  = pr.second;
        if (particle != "gamma" && particle != "neutron" && particle != "all") {
            std::cerr << "[err] unknown particle: " << particle << std::endl; continue;
        }
        fs::create_directories(fs::path(outFile).parent_path());
        TFile* fout = TFile::Open(outFile.c_str(), "RECREATE");
        if (!fout || fout->IsZombie()) {
            std::cerr << "[err] cannot open output: " << outFile << std::endl; continue;
        }

        po.emplace_back();
        PartOut& P   = po.back();
        P.particle   = particle;
        P.cfg        = CutConfig{ det, ParseParticle(particle), (CutLevel)level };
        EnergyWindow(P.cfg, P.eLo, P.eHi);
        P.fout       = fout;
        P.edepDir    = fout->mkdir("edep_per_assembly");
        P.sumDir     = fout->mkdir("edep_sum_per_assembly");

        fout->cd();
        for (const auto& s : sources)
            P.hMaps.push_back(MakeAssemblyHexMap(
                ("h2_counts_" + std::string(s.tag)).c_str(),
                ("Hits counts / assembly (" + std::string(s.tag) + ")").c_str(),
                positions));
        P.hSum = MakeAssemblyHexMap("h2_counts_total",
            "Hits counts / assembly (sum over sources)", positions);

        P.hAll = new TH1D("h1_edep_sum_all_assemblies",
            "edep, sum over all (assemblies, sources);edep [MeV];counts",
            kNbE, kEmin, kEmax);
        P.hAll->SetLineColor(kBlack);
        P.hAllSrc.assign(nsrc, nullptr);
        for (size_t si = 0; si < nsrc; ++si) {
            P.hAllSrc[si] = new TH1D(
                Form("h1_edep_sum_all_assemblies_%s", sources[si].tag),
                Form("edep, all assemblies, %s;edep [MeV];counts", sources[si].tag),
                kNbE, kEmin, kEmax);
            P.hAllSrc[si]->SetLineColor(sources[si].color);
        }
        P.hEk   = new TH2F("h2_edep_vs_ekin", "Edep vs Ekin;Ekin [MeV];Edep [MeV]",
                           1000, 0, 5, 1000, 0, 5);
        P.hWall = new TH2F("h2_primary_zphi_wall",
                           "primary z vs #phi (wall);#phi [rad];z [mm]",
                           720, -M_PI, M_PI, 816, -2040, 2040);
        P.hTop  = new TH2F("h2_primary_xy_top",
                           "primary y vs x (top surface, z=+2040);x [mm];y [mm]",
                           264, -1320, 1320, 264, -1320, 1320);
        P.hBot  = new TH2F("h2_primary_xy_bottom",
                           "primary y vs x (bottom surface, z=-2040);x [mm];y [mm]",
                           264, -1320, 1320, 264, -1320, 1320);
        P.hCos  = new TH1D("h1_cos_alpha",
            "emission-direction alignment;cos#alpha (toward detector);detected primaries",
            220, -1.1, 1.1);
        // Wall spatial coordinate vs deposited energy, filled for EVERY
        // PID-passing wall hit with NO energy gate. This is what lets the gate
        // optimizer measure the (phi,z) spread as a function of a trial window.
        // 10 keV energy bins (coarser than the 1 keV spectrum) keep the file small.
        P.hPhiE = new TH2F("h2_wall_phi_vs_edep",
            "wall #phi vs edep (PID, ungated);edep [MeV];#phi [rad]",
            1000, kEmin, kEmax, 180, -M_PI, M_PI);
        P.hZE   = new TH2F("h2_wall_z_vs_edep",
            "wall z vs edep (PID, ungated);edep [MeV];z [mm]",
            1000, kEmin, kEmax, 204, -2040, 2040);


        P.hcell.assign(kNAssem, std::vector<TH1D*>(nsrc, nullptr));
        P.missing.assign(nsrc, std::vector<bool>(kNAssem, false));
    }
    if (po.empty()) { std::cerr << "[err] no valid outputs to write\n"; return; }

    // ---- discover input files ----------------------------------------------
    std::vector<std::vector<std::string>> filepath;
    DiscoverFiles(inputDir, sources, filepath);
    long long nFound = 0;
    for (auto& row : filepath) for (auto& s : row) if (!s.empty()) ++nFound;
    std::cout << "[info] discovered " << nFound << " files under " << inputDir << "\n";

    const double surfTol = 1.0; // mm
    long long nProcessed = 0;

    // ---- SINGLE pass over the raw data, filling every particle in lockstep --
    for (int g = 0; g < kNAssem; ++g) {
        for (size_t si = 0; si < nsrc; ++si) {
            const std::string& fname = filepath[g][si];

            TFile* fin; TTree* hits; TTree* prim;
            if (fname.empty() || !OpenGood(fname.c_str(), fin, hits, prim)) {
                for (auto& P : po) { P.missing[si][g] = true; ++P.nMissing; }
                std::cout << "[missing] file: \"" << fname.c_str()
                          << "\" not found" << std::endl;
                continue;
            }

            // one spectrum histogram per particle for THIS (g,si) file
            std::vector<TH1D*> h1(po.size(), nullptr);
            for (size_t pi = 0; pi < po.size(); ++pi) {
                po[pi].edepDir->cd();
                TH1D* h = new TH1D(
                    Form("h1_edep_globalFuel%02d_%s", g, sources[si].tag),
                    Form("edep, globalFuel=%d, %s;edep [MeV];counts", g, sources[si].tag),
                    kNbE, kEmin, kEmax);
                h->SetDirectory(po[pi].edepDir);
                h->SetLineColor(sources[si].color);
                h->SetLineStyle(sources[si].line);
                po[pi].hcell[g][si] = h;
                h1[pi] = h;
            }

            // hits branch addresses
            double h_edep, h_evtNb, h_pid;
            hits->SetBranchAddress("edep",  &h_edep);
            hits->SetBranchAddress("evtNb", &h_evtNb);
            hits->SetBranchAddress("pid",   &h_pid);

            // primary branch addresses (shared by all particles)
            double p_ekin=0, p_x=0, p_y=0, p_z=0, p_evtNb=0, p_pid=0;
            double p_px=0, p_py=0, p_pz=0;
            // GUARD: never BuildIndex on an empty tree -> kills the
            // "Cannot build a TreeIndex with a Tree having no entries" spam.
            const bool havePrim = (prim != nullptr && prim->GetEntries() > 0);
            if (havePrim) {
                prim->BuildIndex("evtNb");            // evtNb unique WITHIN file
                prim->SetBranchAddress("ekin",  &p_ekin);
                prim->SetBranchAddress("evtNb", &p_evtNb);
                prim->SetBranchAddress("pid",   &p_pid);
                prim->SetBranchAddress("x",     &p_x);
                prim->SetBranchAddress("y",     &p_y);
                prim->SetBranchAddress("z",     &p_z);
                prim->SetBranchAddress("px",    &p_px);
                prim->SetBranchAddress("py",    &p_py);
                prim->SetBranchAddress("pz",    &p_pz);
            }

            std::vector<long long> nSig(po.size(), 0);
            std::vector<char>      sig(po.size(), 0);
            std::vector<char>      pidPass(po.size(), 0);

            const Long64_t n = hits->GetEntries();
            for (Long64_t i = 0; i < n; ++i) {
                hits->GetEntry(i);
                const long long pid = std::llround(h_pid);
                // full-range spectrum, PID pass, and signal decision, per particle
                bool anySig = false, anyPid = false;
                for (size_t pi = 0; pi < po.size(); ++pi) {
                    pidPass[pi] = PassSpectrum(po[pi].cfg, pid) ? 1 : 0;
                    if (pidPass[pi]) { h1[pi]->Fill(h_edep); anyPid = true; }
                    sig[pi] = PassSignal(po[pi].cfg, pid, h_edep) ? 1 : 0;
                    if (sig[pi]) { ++nSig[pi]; anySig = true; }
                }
                // Gate the primary join on anyPid (NOT anySig): the energy-resolved
                // wall maps must see every PID-passing hit at ALL energies so the
                // optimizer can measure spatial width vs. the trial window.
                if (!anyPid) continue;
                if (havePrim) {
                    const Long64_t pe =
                        prim->GetEntryWithIndex((Long64_t)std::llround(h_evtNb));
                    if (pe < 0) {
                        for (size_t pi = 0; pi < po.size(); ++pi)
                            if (sig[pi]) ++po[pi].nOrphan;
                        continue;
                    }
                    // geometry computed ONCE, shared across particles
                    const double lx  = p_x - caskX,  ly = p_y - caskY;
                    const double r   = std::hypot(lx, ly);
                    const double phi = std::atan2(ly, lx);
                    const double dx = kDetX - p_x, dy = kDetY - p_y, dz = kDetZ - p_z;
                    const double dn = std::sqrt(dx*dx + dy*dy + dz*dz);
                    const double pn = std::sqrt(p_px*p_px + p_py*p_py + p_pz*p_pz);
                    const bool   haveCos  = (dn > 0 && pn > 0);
                    const double cosAlpha = haveCos
                        ? (p_px*dx + p_py*dy + p_pz*dz) / (dn * pn) : 0.0;
                    const bool   onWall   = (std::abs(r - 1270.) < surfTol);
                    for (size_t pi = 0; pi < po.size(); ++pi) {
                        auto& P = po[pi];

                        // energy-RESOLVED wall spatial for every PID-passing hit
                        // (ungated) -> feeds the optimizer's width term.
                        if (pidPass[pi] && onWall) {
                            P.hPhiE->Fill(h_edep, phi);
                            P.hZE  ->Fill(h_edep, p_z);
                        }
                        if (!sig[pi]) continue;   // below: energy-GATED signal fills
                        ++P.nMatched;
                        P.hEk->Fill(p_ekin, h_edep);
                        if      (onWall)                          P.hWall->Fill(phi, p_z);
                        else if (std::abs(p_z - 2040.) < surfTol) P.hTop->Fill(lx, ly);
                        else if (std::abs(p_z + 2040.) < surfTol) P.hBot->Fill(lx, ly);
                        if (haveCos) P.hCos->Fill(cosAlpha);
                    }
                }
            }



            // per-file accumulation into each particle's maps/spectra
            for (size_t pi = 0; pi < po.size(); ++pi) {
                auto& P = po[pi];
                P.hMaps[si]->SetBinContent(g + 1, nSig[pi]);
                P.hSum->SetBinContent(g + 1, P.hSum->GetBinContent(g + 1) + nSig[pi]);
                if (nSig[pi] == 0) ++P.nEmpty;
                P.hAll->Add(h1[pi]);
                P.hAllSrc[si]->Add(h1[pi]);
            }

            ++nProcessed;
            fin->Close(); delete fin;
        }
    }

    // ---- per-particle finalize + write (identical layout to the old file) --
    for (auto& P : po) {
        P.fout->cd();

        std::vector<bool> missingSum(kNAssem, false);
        for (int g = 0; g < kNAssem; ++g) {
            bool all = true;
            for (size_t si = 0; si < nsrc; ++si)
                if (!P.missing[si][g]) { all = false; break; }
            missingSum[g] = all;
        }

        // shared colour scale across the 6 hex maps
        double zmin = 0, zmax = -1e30;
        auto scan = [&](TH2Poly* h) {
            for (int b = 1; b <= h->GetNumberOfBins(); ++b) {
                double v = h->GetBinContent(b); if (v > zmax) zmax = v;
            }
        };
        for (auto* h : P.hMaps) scan(h); scan(P.hSum);
        if (!(zmax > zmin)) zmax = 1.0;
        for (auto* h : P.hMaps) { h->SetMinimum(zmin); h->SetMaximum(zmax); }
        P.hSum->SetMinimum(zmin); P.hSum->SetMaximum(zmax);

        // 2x3 hex-map canvas
        auto* cH = new TCanvas(Form("c_assembly_heatmaps_%s", P.particle.c_str()),
                               "Per-assembly hits counts", 1400, 800);
        cH->Divide(3, 2);
        auto drawOne = [&](TVirtualPad* pad, TH2Poly* h, const std::vector<bool>& miss) {
            pad->cd(); pad->SetRightMargin(0.15); pad->SetLeftMargin(0.12);
            h->Draw("TEXT COLZ L");
            DrawMissingOverlay(positions, miss);
            h->Write();
        };
        for (size_t i = 0; i < nsrc; ++i) drawOne(cH->cd(i + 1), P.hMaps[i], P.missing[i]);
        drawOne(cH->cd(6), P.hSum, missingSum);
        cH->Write();

        // per-assembly summed spectra (+ optional canvases)
        auto seedSum = [&](int g) -> TH1D* {
            TH1D* s = nullptr;
            for (size_t si = 0; si < nsrc; ++si) if (P.hcell[g][si]) {
                s = (TH1D*)P.hcell[g][si]->Clone(Form("h1_edep_sum_globalFuel%02d", g));
                s->Reset(); s->SetDirectory(P.sumDir);
                s->SetTitle(Form("edep, globalFuel=%d;edep [MeV];counts", g));
                s->SetLineColor(kBlack); s->SetLineWidth(1);
                for (size_t sj = 0; sj < nsrc; ++sj)
                    if (P.hcell[g][sj]) s->Add(P.hcell[g][sj]);
                break;
            }
            return s;
        };
        for (int g = 0; g < kNAssem; ++g) {
            TH1D* hs = seedSum(g);
            if (!hs) continue;
            if (!makePerAssemblyCanvases) continue;

            auto* c = new TCanvas(Form("c_spectra_globalFuel%02d_%s", g, P.particle.c_str()),
                                  Form("edep spectra, globalFuel=%d", g), 1100, 750);
            c->SetLogy(); c->SetGrid();
            double ymin = 0.5, ymax = 5.0 * std::max(1.0, hs->GetMaximum());
            hs->SetMinimum(ymin); hs->SetMaximum(ymax);
            hs->Draw("HIST");
            for (size_t si = 0; si < nsrc; ++si)
                if (P.hcell[g][si]) P.hcell[g][si]->Draw("HIST SAME");
            hs->Draw("HIST SAME");
            auto* leg = new TLegend(0.62, 0.62, 0.96, 0.90); leg->SetTextSize(0.03);
            leg->AddEntry(hs, "sum of sources", "l");
            for (size_t si = 0; si < nsrc; ++si) if (P.hcell[g][si])
                leg->AddEntry(P.hcell[g][si],
                    Form("%s (%.3g cts)", PrettyLabel(sources[si].tag),
                         P.hcell[g][si]->Integral()), "l");
            for (size_t si = 0; si < nsrc; ++si) if (!P.hcell[g][si]) {
                auto* d = new TH1D(Form("h1_stub_%02d_%zu_%s", g, si, P.particle.c_str()),
                                   "", 1, 0, 1);
                d->SetLineColor(sources[si].color); d->SetLineStyle(2);
                leg->AddEntry(d, Form("%s [missing file]",
                                      PrettyLabel(sources[si].tag)), "l");
            }
            leg->Draw();
            if (P.cfg.level == kPidEnergy) DrawWindow(P.eLo, P.eHi, ymin, ymax);
            auto* t = new TLatex(); t->SetNDC(true); t->SetTextSize(0.035);
            t->DrawLatex(0.10, 0.945,
                Form("globalFuel = %d   (total = %.3g cts)", g, hs->Integral()));
            c->Write();
            delete c;
        }

        // grand-sum spectrum canvas
        P.fout->cd();
        auto* cAll = new TCanvas(Form("c_edep_sum_all_assemblies_%s", P.particle.c_str()),
                                 "edep, summed over all assemblies", 1100, 750);
        cAll->SetLogy(); cAll->SetGrid();
        if (P.hAll->GetEntries() > 0) {
            double ymin = 0.5, ymax = 5.0 * std::max(1.0, P.hAll->GetMaximum());
            P.hAll->SetMinimum(ymin); P.hAll->SetMaximum(ymax);
            P.hAll->Draw("HIST");
            for (size_t si = 0; si < nsrc; ++si) P.hAllSrc[si]->Draw("HIST SAME");
            P.hAll->Draw("HIST SAME");
            auto* leg = new TLegend(0.62, 0.62, 0.96, 0.90); leg->SetTextSize(0.03);
            leg->AddEntry(P.hAll, Form("sum (%.3g cts)", P.hAll->Integral()), "l");
            for (size_t si = 0; si < nsrc; ++si)
                leg->AddEntry(P.hAllSrc[si], Form("%s (%.3g cts)",
                    PrettyLabel(sources[si].tag), P.hAllSrc[si]->Integral()), "l");
            leg->Draw();
            if (P.cfg.level == kPidEnergy) DrawWindow(P.eLo, P.eHi, ymin, ymax);
        } else {
            std::cerr << "[warn] no spectrum entries (" << P.particle
                      << "); grand-sum canvas empty\n";
        }
        cAll->Write();

        // write remaining objects (hCos/hEmit ride along via fout->Write)
        P.hAll->Write();
        for (auto* h : P.hAllSrc) h->Write();
        P.hEk->Write(); P.hWall->Write(); P.hTop->Write(); P.hBot->Write();
        P.hCos->Write(); 
        P.hPhiE->Write(); P.hZE->Write();

        TNamed cfgStamp("cut_config",
            Form("detector=%s particle=%s level=%d window=[%.4g,%.4g]MeV",
                 det == kCLYC ? "CLYC" : "plastic",
                 ParticleName(P.cfg.part), (int)P.cfg.level, P.eLo, P.eHi));
        cfgStamp.Write();
        
        TParameter<Long64_t> pIn ("n_input_files",   nFound);      pIn.Write();
        TParameter<Long64_t> pMis("n_missing_files", P.nMissing);  pMis.Write();

        P.fout->Write();
        P.fout->Close();

        std::cout << kAnsiGreen << kAnsiBold << "[ok] " << kAnsiReset
                  << P.particle << ": processed " << kAnsiBold << nProcessed << kAnsiReset
                  << " (fuel,source) pairs; "
                  << kAnsiYellow << P.nEmpty  << kAnsiReset << " zero-signal, "
                  << kAnsiRed    << P.nMissing << kAnsiReset << " missing files; "
                  << "primary match=" << P.nMatched << " orphan=" << P.nOrphan << ".\n";
    }
}

