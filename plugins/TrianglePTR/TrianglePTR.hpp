// PluginTrianglePTR.hpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"

namespace TrianglePTR {

class TrianglePTR : public SCUnit {
public:
    TrianglePTR();

    // Destructor
    // ~TrianglePTR();

private:
    // Calc functions
    void next_aaaa(int nSamples);
    void next_aaak(int nSamples);
    void next_aaka(int nSamples);
    void next_aakk(int nSamples);
    void next_akaa(int nSamples);
    void next_akak(int nSamples);
    void next_akka(int nSamples);
    void next_akkk(int nSamples);
    void next_kaaa(int nSamples);
    void next_kaak(int nSamples);
    void next_kaka(int nSamples);
    void next_kakk(int nSamples);
    void next_kkaa(int nSamples);
    void next_kkak(int nSamples);
    void next_kkka(int nSamples);
    void next_kkkk(int nSamples);
    double algorithm(double p, double t0, double t2, double t3,
            double w, double a, double b, double c,
            double dc, double p1, double p2, double p3);

    // Member variables
    double mPhase   = 0;
    double mPhaseIn = 0;
    double mSync    = 0;
};

} // namespace TrianglePTR
