// analyze_cask.C
//
// Full legacy/default mode:
//   analyze_cask(..., makePerAssemblyCanvases=false)
//
// Optional fast/slim mode:
//   analyze_cask(..., makePerAssemblyCanvases=false, slimOutput=true)
//
// Full mode preserves the legacy output layout.
// Slim mode writes only the objects needed by plot_design_comp's count/spectrum
// summaries, plus compact multidet detector-position heatmaps when applicable.

#include "geometry_constraints.h"
#include "dcs_geometry.h"
#include "dcs_cuts.h"
#include "dcs_sources.h"
#include "dcs_setup.h"
#include "dcs_multidet.h"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <utility>
#include <filesystem>
#include <regex>
#include <cerrno>
#include <cstring>
#include <algorithm>

namespace fs = std::filesystem;

using namespace geo;
using namespace dcs;

namespace {

static int CountOpenFDs()
{
    int n = 0;
    std::error_code ec;

    for (auto it = fs::directory_iterator("/proc/self/fd", ec);
         !ec && it != fs::directory_iterator();
         it.increment(ec)) {
        ++n;
    }

    return n;
}

bool OpenGood(const char* fname, TFile*& fin, TTree*& hits, TTree*& prim)
{
    fin  = nullptr;
    hits = nullptr;
    prim = nullptr;

    errno = 0;
    fin = TFile::Open(fname, "READ");

    const int  e      = errno;
    const bool isnull = (fin == nullptr);
    const bool zomb   = (fin && fin->IsZombie());
    const bool recov  = (fin && fin->TestBit(TFile::kRecovered));

    if (fin && !zomb && !recov) {
        hits = (TTree*)fin->Get("hits");
    }

    if (hits) {
        prim = (TTree*)fin->Get("primary");
        return true;
    }

    std::cerr << "[OpenGood FAIL] " << fname
              << "  null="   << isnull
              << " zombie="  << zomb
              << " recov="   << recov
              << " hasHits=" << (hits != nullptr)
              << " errno="   << e << "(" << std::strerror(e) << ")"
              << " openFDs=" << CountOpenFDs()
              << "\n";

    if (fin) {
        fin->Close();
        delete fin;
        fin = nullptr;
    }

    hits = nullptr;
    prim = nullptr;
    return false;
}

void DrawWindow(double lo, double hi, double ymin, double ymax)
{
    for (double x : {lo, hi}) {
        auto* l = new TLine(x, ymin, x, ymax);
        l->SetLineColor(kGray + 2);
        l->SetLineStyle(2);
        l->SetLineWidth(2);
        l->Draw();
    }

    auto* t = new TLatex();
    t->SetTextSize(0.028);
    t->SetTextColor(kGray + 2);
    t->DrawLatex(hi, ymax * 0.4, Form("  gate [%.3g, %.3g] MeV", lo, hi));
}

// Discover detector-response files in 'dir' and map them onto (globalFuel, src).
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
        if (name.rfind("detector-response", 0) != 0) continue;
        if (e.path().extension() != ".root") continue;

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

        out[g][si] = e.path().string();
    }
}

void PrintFileProgress(const char* inputDir,
                       long long slot,
                       long long total,
                       int g,
                       const char* sourceTag,
                       bool slimOutput)
{
    std::cout << "\33[2K\r"
              << "[analyze] " << inputDir
              << " | file " << slot << "/" << total
              << " | globalFuel=" << Form("%02d", g)
              << " | source=" << sourceTag
              << " | mode=" << (slimOutput ? "slim" : "full")
              << std::flush;
}

} // anonymous namespace



void analyze_cask(const char* inputDir = "data/nominal/cask0",
                  const std::vector<std::pair<std::string,std::string>>& outputs
                      = { {"gamma",   "analysis_gamma.root"},
                          {"neutron", "analysis_neutron.root"} },
                  DetectorType det     = kCLYC,
                  int level            = 2,
                  double caskX         = 0.0,
                  double caskY         = 0.0,
                  bool makePerAssemblyCanvases = false,
                  bool slimOutput = false)
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(256);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);

    const auto& sources    = Sources();
    const size_t nsrc      = sources.size();
    const auto   positions = BuildFuelPositions();

    const bool multiDetMode = LooksLikeMultidetPath(inputDir);

    if ((int)positions.size() != kNAssem) {
        std::cerr << "[err] expected " << kNAssem << " positions, got "
                  << positions.size() << std::endl;
        return;
    }

    if (slimOutput) {
        std::cout << "[info] analyze_cask slim mode enabled for " << inputDir
                  << " -- skipping primary-join diagnostics and per-assembly spectra\n";
    }

    if (multiDetMode) {
        std::cout << "[info] multidet compact detector-position maps enabled for "
                  << inputDir << "\n";
    }

    struct PartOut {
        std::string particle;
        CutConfig   cfg;
        double      eLo = 0.0;
        double      eHi = 0.0;

        TFile*      fout = nullptr;
        TDirectory* edepDir = nullptr;
        TDirectory* sumDir  = nullptr;

        std::vector<TH2Poly*> hMaps;
        TH2Poly*    hSum = nullptr;

        TH1D*       hAll = nullptr;
        std::vector<TH1D*> hAllSrc;

        // Full-mode diagnostics.
        TH2F *hEk=nullptr, *hWall=nullptr, *hTop=nullptr, *hBot=nullptr;
        TH2F *hPhiE=nullptr, *hZE=nullptr;
        TH1D *hCos=nullptr;

        std::vector<std::vector<TH1D*>> hcell;
        std::vector<std::vector<bool>>  missing;

        // Compact multidet outputs.
        bool multiDet = false;
        std::vector<TH2F*> hDetPhiZSrc;
        TH2F* hDetPhiZTotal = nullptr;

        long long nEmpty=0, nMissing=0, nMatched=0, nOrphan=0;
    };

    std::vector<PartOut> po;
    po.reserve(outputs.size());

    for (const auto& pr : outputs) {
        const std::string& particle = pr.first;
        const std::string& outFile  = pr.second;

        if (particle != "gamma" && particle != "neutron" && particle != "all") {
            std::cerr << "[err] unknown particle: " << particle << std::endl;
            continue;
        }

        fs::create_directories(fs::path(outFile).parent_path());

        TFile* fout = TFile::Open(outFile.c_str(), "RECREATE");
        if (!fout || fout->IsZombie()) {
            std::cerr << "[err] cannot open output: " << outFile << std::endl;
            continue;
        }

        po.emplace_back();
        PartOut& P = po.back();

        P.particle = particle;
        P.cfg      = CutConfig{ det, ParseParticle(particle), (CutLevel)level };
        EnergyWindow(P.cfg, P.eLo, P.eHi);
        P.fout     = fout;

        if (!slimOutput) {
            P.edepDir = fout->mkdir("edep_per_assembly");
            P.sumDir  = fout->mkdir("edep_sum_per_assembly");
        }

        fout->cd();

        for (const auto& s : sources) {
            P.hMaps.push_back(MakeAssemblyHexMap(
                ("h2_counts_" + std::string(s.tag)).c_str(),
                ("Hits counts / assembly (" + std::string(s.tag) + ")").c_str(),
                positions));
        }

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

        if (!slimOutput) {
            P.hEk = new TH2F("h2_edep_vs_ekin",
                             "Edep vs Ekin;Ekin [MeV];Edep [MeV]",
                             1000, 0, 5, 1000, 0, 5);

            P.hWall = new TH2F("h2_primary_zphi_wall",
                               "primary z vs #phi (wall);#phi [rad];z [mm]",
                               720, -M_PI, M_PI, 816, -2040, 2040);

            P.hTop = new TH2F("h2_primary_xy_top",
                              "primary y vs x (top surface, z=+2040);x [mm];y [mm]",
                              264, -1320, 1320, 264, -1320, 1320);

            P.hBot = new TH2F("h2_primary_xy_bottom",
                              "primary y vs x (bottom surface, z=-2040);x [mm];y [mm]",
                              264, -1320, 1320, 264, -1320, 1320);

            P.hCos = new TH1D("h1_cos_alpha",
                "emission-direction alignment;cos#alpha (toward detector);detected primaries",
                220, -1.1, 1.1);

            P.hPhiE = new TH2F("h2_wall_phi_vs_edep",
                "wall #phi vs edep (PID, ungated);edep [MeV];#phi [rad]",
                1000, kEmin, kEmax, 180, -M_PI, M_PI);

            P.hZE = new TH2F("h2_wall_z_vs_edep",
                "wall z vs edep (PID, ungated);edep [MeV];z [mm]",
                1000, kEmin, kEmax, 204, -2040, 2040);

            P.hcell.assign(kNAssem, std::vector<TH1D*>(nsrc, nullptr));
        }

        P.missing.assign(nsrc, std::vector<bool>(kNAssem, false));

        P.multiDet = multiDetMode;

        if (P.multiDet) {
            P.hDetPhiZSrc.assign(nsrc, nullptr);

            for (size_t si = 0; si < nsrc; ++si) {
                P.hDetPhiZSrc[si] = MakeDetectorPhiZHist(
                    Form("h2_detector_phiz_counts_%s", sources[si].tag),
                    Form("per-cask detector-position ROI counts, %s;"
                         "detector #phi [deg];detector z [mm];counts",
                         sources[si].tag));
            }

            P.hDetPhiZTotal = MakeDetectorPhiZHist(
                "h2_detector_phiz_counts_total",
                "per-cask detector-position ROI counts, sum over sources;"
                "detector #phi [deg];detector z [mm];counts");
        }
    }

    if (po.empty()) {
        std::cerr << "[err] no valid outputs to write\n";
        return;
    }

    std::vector<std::vector<std::string>> filepath;
    DiscoverFiles(inputDir, sources, filepath);

    long long nFound = 0;
    for (auto& row : filepath) {
        for (auto& s : row) {
            if (!s.empty()) ++nFound;
        }
    }

    const long long nExpectedSlots =
        static_cast<long long>(kNAssem) * static_cast<long long>(nsrc);

    std::cout << "[info] discovered " << nFound
              << " files under " << inputDir << "\n";

    const double surfTol = 1.0;
    long long nProcessed = 0;
    long long slot = 0;

    for (int g = 0; g < kNAssem; ++g) {
        for (size_t si = 0; si < nsrc; ++si) {
            ++slot;

            PrintFileProgress(inputDir, slot, nExpectedSlots,
                              g, sources[si].tag, slimOutput);

            const std::string& fname = filepath[g][si];

            TFile* fin = nullptr;
            TTree* hits = nullptr;
            TTree* prim = nullptr;

            if (fname.empty() || !OpenGood(fname.c_str(), fin, hits, prim)) {
                for (auto& P : po) {
                    P.missing[si][g] = true;
                    ++P.nMissing;
                }
                continue;
            }

            std::vector<TH1D*> h1(po.size(), nullptr);

            if (!slimOutput) {
                for (size_t pi = 0; pi < po.size(); ++pi) {
                    po[pi].edepDir->cd();

                    TH1D* h = new TH1D(
                        Form("h1_edep_globalFuel%02d_%s", g, sources[si].tag),
                        Form("edep, globalFuel=%d, %s;edep [MeV];counts",
                             g, sources[si].tag),
                        kNbE, kEmin, kEmax);

                    h->SetDirectory(po[pi].edepDir);
                    h->SetLineColor(sources[si].color);
                    h->SetLineStyle(sources[si].line);

                    po[pi].hcell[g][si] = h;
                    h1[pi] = h;
                }
            }

            double h_edep  = 0.0;
            double h_evtNb = 0.0;
            double h_pid   = 0.0;
            Int_t  h_det   = 0;

            hits->SetBranchStatus("*", 0);
            hits->SetBranchStatus("edep",  1);
            hits->SetBranchStatus("evtNb", 1);
            hits->SetBranchStatus("pid",   1);

            const bool haveDetBranch = (hits->GetBranch("det") != nullptr);
            if (haveDetBranch) hits->SetBranchStatus("det", 1);

            hits->SetBranchAddress("edep",  &h_edep);
            hits->SetBranchAddress("evtNb", &h_evtNb);
            hits->SetBranchAddress("pid",   &h_pid);
            if (haveDetBranch) hits->SetBranchAddress("det", &h_det);

            double p_ekin=0, p_x=0, p_y=0, p_z=0, p_evtNb=0, p_pid=0;
            double p_px=0, p_py=0, p_pz=0;

            const bool havePrim =
                (!slimOutput && prim != nullptr && prim->GetEntries() > 0);

            if (havePrim) {
                prim->SetBranchStatus("*", 0);
                prim->SetBranchStatus("ekin",  1);
                prim->SetBranchStatus("evtNb", 1);
                prim->SetBranchStatus("pid",   1);
                prim->SetBranchStatus("x",     1);
                prim->SetBranchStatus("y",     1);
                prim->SetBranchStatus("z",     1);
                prim->SetBranchStatus("px",    1);
                prim->SetBranchStatus("py",    1);
                prim->SetBranchStatus("pz",    1);

                prim->BuildIndex("evtNb");

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

            std::vector<std::vector<long long>> nSigDet(
                po.size(), std::vector<long long>(kNMultidet, 0));

            std::vector<char> sig(po.size(), 0);
            std::vector<char> pidPass(po.size(), 0);

            const Long64_t n = hits->GetEntries();

            for (Long64_t i = 0; i < n; ++i) {
                hits->GetEntry(i);

                const long long pid = std::llround(h_pid);

                int detIdx = 0;
                bool validMultidetHit = false;

                if (haveDetBranch) {
                    detIdx = static_cast<int>(h_det);
                    validMultidetHit = IsValidMultidetId(detIdx);
                }

                if (!validMultidetHit) detIdx = 0;

                bool anyPid = false;

                for (size_t pi = 0; pi < po.size(); ++pi) {
                    auto& P = po[pi];

                    pidPass[pi] = PassSpectrum(P.cfg, pid) ? 1 : 0;

                    if (pidPass[pi]) {
                        if (slimOutput) {
                            P.hAll->Fill(h_edep);
                            P.hAllSrc[si]->Fill(h_edep);
                        } else {
                            h1[pi]->Fill(h_edep);
                        }

                        anyPid = true;
                    }

                    sig[pi] = PassSignal(P.cfg, pid, h_edep) ? 1 : 0;

                    if (sig[pi]) {
                        ++nSig[pi];

                        if (P.multiDet && validMultidetHit) {
                            ++nSigDet[pi][detIdx];
                        }
                    }
                }

                if (slimOutput) continue;

                if (!anyPid) continue;

                if (havePrim) {
                    const Long64_t pe =
                        prim->GetEntryWithIndex((Long64_t)std::llround(h_evtNb));

                    if (pe < 0) {
                        for (size_t pi = 0; pi < po.size(); ++pi) {
                            if (sig[pi]) ++po[pi].nOrphan;
                        }
                        continue;
                    }

                    const double lx  = p_x - caskX;
                    const double ly  = p_y - caskY;
                    const double r   = std::hypot(lx, ly);
                    const double phi = std::atan2(ly, lx);

                    const double dx = kDetX - p_x;
                    const double dy = kDetY - p_y;
                    const double dz = kDetZ - p_z;

                    const double dn = std::sqrt(dx*dx + dy*dy + dz*dz);
                    const double pn = std::sqrt(p_px*p_px + p_py*p_py + p_pz*p_pz);

                    const bool haveCos = (dn > 0 && pn > 0);
                    const double cosAlpha = haveCos
                        ? (p_px*dx + p_py*dy + p_pz*dz) / (dn * pn)
                        : 0.0;

                    const bool onWall = (std::abs(r - 1270.) < surfTol);

                    for (size_t pi = 0; pi < po.size(); ++pi) {
                        auto& P = po[pi];

                        if (pidPass[pi] && onWall) {
                            P.hPhiE->Fill(h_edep, phi);
                            P.hZE  ->Fill(h_edep, p_z);
                        }

                        if (!sig[pi]) continue;

                        ++P.nMatched;

                        P.hEk->Fill(p_ekin, h_edep);

                        if      (onWall)                          P.hWall->Fill(phi, p_z);
                        else if (std::abs(p_z - 2040.) < surfTol) P.hTop->Fill(lx, ly);
                        else if (std::abs(p_z + 2040.) < surfTol) P.hBot->Fill(lx, ly);

                        if (haveCos) P.hCos->Fill(cosAlpha);
                    }
                }
            }

            for (size_t pi = 0; pi < po.size(); ++pi) {
                auto& P = po[pi];

                P.hMaps[si]->SetBinContent(g + 1, nSig[pi]);
                P.hSum->SetBinContent(g + 1,
                    P.hSum->GetBinContent(g + 1) + nSig[pi]);

                if (nSig[pi] == 0) ++P.nEmpty;

                if (!slimOutput) {
                    P.hAll->Add(h1[pi]);
                    P.hAllSrc[si]->Add(h1[pi]);
                }

                if (P.multiDet) {
                    for (int d = 0; d < kNMultidet; ++d) {
                        const long long v = nSigDet[pi][d];
                        if (v == 0) continue;

                        AddDetectorPhiZBin(P.hDetPhiZSrc[si], d, v);
                        AddDetectorPhiZBin(P.hDetPhiZTotal,  d, v);
                    }
                }
            }

            ++nProcessed;

            fin->Close();
            delete fin;
        }
    }

    std::cout << "\33[2K\r"
              << "[analyze] completed " << inputDir
              << " | processed " << nProcessed << "/" << nFound
              << " discovered files"
              << " | mode=" << (slimOutput ? "slim" : "full")
              << "\n";

    for (auto& P : po) {
        P.fout->cd();

        std::vector<bool> missingSum(kNAssem, false);

        for (int g = 0; g < kNAssem; ++g) {
            bool all = true;
            for (size_t si = 0; si < nsrc; ++si) {
                if (!P.missing[si][g]) {
                    all = false;
                    break;
                }
            }
            missingSum[g] = all;
        }

        double zmin = 0.0;
        double zmax = -1e30;

        auto scan = [&](TH2Poly* h) {
            for (int b = 1; b <= h->GetNumberOfBins(); ++b) {
                double v = h->GetBinContent(b);
                if (v > zmax) zmax = v;
            }
        };

        for (auto* h : P.hMaps) scan(h);
        scan(P.hSum);

        if (!(zmax > zmin)) zmax = 1.0;

        for (auto* h : P.hMaps) {
            h->SetMinimum(zmin);
            h->SetMaximum(zmax);
        }

        P.hSum->SetMinimum(zmin);
        P.hSum->SetMaximum(zmax);

        if (!slimOutput) {
            auto* cH = new TCanvas(Form("c_assembly_heatmaps_%s", P.particle.c_str()),
                                   "Per-assembly hits counts", 1400, 800);
            cH->Divide(3, 2);

            auto drawOne = [&](TVirtualPad* pad,
                               TH2Poly* h,
                               const std::vector<bool>& miss) {
                pad->cd();
                pad->SetRightMargin(0.15);
                pad->SetLeftMargin(0.12);
                h->Draw("TEXT COLZ L");
                DrawMissingOverlay(positions, miss);
                h->Write();
            };

            for (size_t i = 0; i < nsrc; ++i) {
                drawOne(cH->cd(i + 1), P.hMaps[i], P.missing[i]);
            }

            drawOne(cH->cd(6), P.hSum, missingSum);
            cH->Write();

            auto seedSum = [&](int g) -> TH1D* {
                TH1D* s = nullptr;

                for (size_t si = 0; si < nsrc; ++si) {
                    if (P.hcell[g][si]) {
                        s = (TH1D*)P.hcell[g][si]->Clone(
                            Form("h1_edep_sum_globalFuel%02d", g));
                        s->Reset();
                        s->SetDirectory(P.sumDir);
                        s->SetTitle(Form("edep, globalFuel=%d;edep [MeV];counts", g));
                        s->SetLineColor(kBlack);
                        s->SetLineWidth(1);

                        for (size_t sj = 0; sj < nsrc; ++sj) {
                            if (P.hcell[g][sj]) s->Add(P.hcell[g][sj]);
                        }

                        break;
                    }
                }

                return s;
            };

            for (int g = 0; g < kNAssem; ++g) {
                TH1D* hs = seedSum(g);
                if (!hs) continue;
                if (!makePerAssemblyCanvases) continue;

                auto* c = new TCanvas(
                    Form("c_spectra_globalFuel%02d_%s", g, P.particle.c_str()),
                    Form("edep spectra, globalFuel=%d", g), 1100, 750);

                c->SetLogy();
                c->SetGrid();

                double ymin = 0.5;
                double ymax = 5.0 * std::max(1.0, hs->GetMaximum());

                hs->SetMinimum(ymin);
                hs->SetMaximum(ymax);
                hs->Draw("HIST");

                for (size_t si = 0; si < nsrc; ++si) {
                    if (P.hcell[g][si]) P.hcell[g][si]->Draw("HIST SAME");
                }

                hs->Draw("HIST SAME");

                auto* leg = new TLegend(0.62, 0.62, 0.96, 0.90);
                leg->SetTextSize(0.03);
                leg->AddEntry(hs, "sum of sources", "l");

                for (size_t si = 0; si < nsrc; ++si) {
                    if (P.hcell[g][si]) {
                        leg->AddEntry(P.hcell[g][si],
                            Form("%s (%.3g cts)", PrettyLabel(sources[si].tag),
                                 P.hcell[g][si]->Integral()), "l");
                    }
                }

                for (size_t si = 0; si < nsrc; ++si) {
                    if (!P.hcell[g][si]) {
                        auto* d = new TH1D(
                            Form("h1_stub_%02d_%zu_%s", g, si, P.particle.c_str()),
                            "", 1, 0, 1);
                        d->SetLineColor(sources[si].color);
                        d->SetLineStyle(2);
                        leg->AddEntry(d, Form("%s [missing file]",
                                              PrettyLabel(sources[si].tag)), "l");
                    }
                }

                leg->Draw();

                if (P.cfg.level == kPidEnergy) DrawWindow(P.eLo, P.eHi, ymin, ymax);

                auto* t = new TLatex();
                t->SetNDC(true);
                t->SetTextSize(0.035);
                t->DrawLatex(0.10, 0.945,
                    Form("globalFuel = %d   (total = %.3g cts)", g, hs->Integral()));

                c->Write();
                delete c;
            }

            P.fout->cd();

            auto* cAll = new TCanvas(
                Form("c_edep_sum_all_assemblies_%s", P.particle.c_str()),
                "edep, summed over all assemblies", 1100, 750);

            cAll->SetLogy();
            cAll->SetGrid();

            if (P.hAll->GetEntries() > 0) {
                double ymin = 0.5;
                double ymax = 5.0 * std::max(1.0, P.hAll->GetMaximum());

                P.hAll->SetMinimum(ymin);
                P.hAll->SetMaximum(ymax);
                P.hAll->Draw("HIST");

                for (size_t si = 0; si < nsrc; ++si) {
                    P.hAllSrc[si]->Draw("HIST SAME");
                }

                P.hAll->Draw("HIST SAME");

                auto* leg = new TLegend(0.62, 0.62, 0.96, 0.90);
                leg->SetTextSize(0.03);
                leg->AddEntry(P.hAll, Form("sum (%.3g cts)", P.hAll->Integral()), "l");

                for (size_t si = 0; si < nsrc; ++si) {
                    leg->AddEntry(P.hAllSrc[si], Form("%s (%.3g cts)",
                        PrettyLabel(sources[si].tag), P.hAllSrc[si]->Integral()), "l");
                }

                leg->Draw();

                if (P.cfg.level == kPidEnergy) DrawWindow(P.eLo, P.eHi, ymin, ymax);
            } else {
                std::cerr << "[warn] no spectrum entries (" << P.particle
                          << "); grand-sum canvas empty\n";
            }

            cAll->Write();
        }

        if (slimOutput) {
            for (auto* h : P.hMaps) h->Write();
            P.hSum->Write();
        }

        P.hAll->Write();
        for (auto* h : P.hAllSrc) h->Write();

        if (!slimOutput) {
            P.hEk->Write();
            P.hWall->Write();
            P.hTop->Write();
            P.hBot->Write();
            P.hCos->Write();
            P.hPhiE->Write();
            P.hZE->Write();
        }

        if (P.multiDet) {
            for (auto* h : P.hDetPhiZSrc) {
                if (h) h->Write();
            }

            if (P.hDetPhiZTotal) P.hDetPhiZTotal->Write();
        }

        TNamed cfgStamp("cut_config",
            Form("detector=%s particle=%s level=%d window=[%.4g,%.4g]MeV",
                 det == kCLYC ? "CLYC" : "plastic",
                 ParticleName(P.cfg.part), (int)P.cfg.level, P.eLo, P.eHi));
        cfgStamp.Write();

        TParameter<Long64_t> pIn ("n_input_files",   nFound);
        TParameter<Long64_t> pMis("n_missing_files", P.nMissing);
        pIn.Write();
        pMis.Write();

        TParameter<Int_t> pSlim("analysis_slim_output", slimOutput ? 1 : 0);
        pSlim.Write();

        TParameter<Int_t> pMD("is_multidet_analysis", P.multiDet ? 1 : 0);
        pMD.Write();

        TParameter<Int_t> pND("n_multidet_positions",
                              P.multiDet ? kNMultidet : 1);
        pND.Write();

        if (!slimOutput) {
            P.fout->Write();
        }

        P.fout->Close();

        std::cout << kAnsiGreen << kAnsiBold << "[ok] " << kAnsiReset
                  << P.particle << ": processed " << kAnsiBold << nProcessed << kAnsiReset
                  << " (fuel,source) pairs; "
                  << kAnsiYellow << P.nEmpty  << kAnsiReset << " zero-signal, "
                  << kAnsiRed    << P.nMissing << kAnsiReset << " missing files; "
                  << "primary match=" << P.nMatched
                  << " orphan=" << P.nOrphan
                  << "; mode=" << (slimOutput ? "slim" : "full")
                  << ".\n";
    }
}

