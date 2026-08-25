#ifndef GEOM_CONSTRAINTS
#define GEOM_CONSTRAINTS

namespace geo {
    inline constexpr int    kNAssem   = 84;     // C++17 inline → safe in many TUs
    inline constexpr double kPitch    = 147.0;
    inline constexpr int    kNRings   = 5;
    inline constexpr double kPhiStart = 30.0 * M_PI/180.;

    inline constexpr int    kNbE  = 10000;
    inline constexpr double kEmin = 0.0;
    inline constexpr double kEmax = 10.0;

    enum DetectorType { kCLYC, kPlastic };
}

#endif
