#ifndef DCS_CUTS_H
#define DCS_CUTS_H
#include "geometry_constraints.h"
#include <cmath>
#include <string>

namespace dcs {
using geo::DetectorType; using geo::kCLYC; using geo::kPlastic;

enum Particle { kGamma, kNeutron, kAll };
enum CutLevel { kNone = 0, kPidOnly = 1, kPidEnergy = 2 };

struct CutConfig {
    DetectorType det   = kCLYC;
    Particle     part  = kGamma;
    CutLevel     level = kPidEnergy;
};

inline Particle ParseParticle(const std::string& s) {
    if (s == "gamma")   return kGamma;
    if (s == "neutron") return kNeutron;
    return kAll;
}
inline const char* ParticleName(Particle p) {
    return p == kGamma ? "gamma" : p == kNeutron ? "neutron" : "all";
}

// Energy window [MeV], keyed on (detector, particle). Single definition,
// consolidating the values previously duplicated across the three macros.
inline void EnergyWindow(const CutConfig& c, double& lo, double& hi) {
    if (c.part == kGamma)        { lo = 0.500; hi = 0.740; }
    else if (c.part == kNeutron) {
        if (c.det == kCLYC)      { lo = 4.6;   hi = 5.2;   }  // Li-6(n,a) peak
        else                     { lo = 0.500; hi = 5.000; }  // plastic
    } else                       { lo = geo::kEmin; hi = geo::kEmax; }
}

// PID selection on the HITS-tree energy-depositing particle (this is how the
// real detector discriminates gamma vs. thermal-n capture -> e-/e+ vs. alpha).
inline bool PassPid(const CutConfig& c, long long pid) {
    if (c.part == kAll)   return true;
    if (c.part == kGamma) return (pid == 11 || pid == 22 || pid == -11);
    if (c.det == kCLYC)   return (pid == 1000020040 || pid == 2112);
    return (pid == 1000060120 || pid == 2112 || pid == 2212);   // plastic
}

// Entry contributes to the (full-range) SPECTRUM histogram.
inline bool PassSpectrum(const CutConfig& c, long long pid) {
    return (c.level == kNone) ? true : PassPid(c, pid);
}

// Entry counts as a detector SIGNAL (heatmap bin + primary maps).
inline bool PassSignal(const CutConfig& c, long long pid, double edep) {
    if (c.level == kNone) return true;
    if (!PassPid(c, pid)) return false;
    if (c.level == kPidEnergy) {
        double lo, hi; EnergyWindow(c, lo, hi);
        if (edep < lo || edep > hi) return false;
    }
    return true;
}
} // namespace dcs
#endif

