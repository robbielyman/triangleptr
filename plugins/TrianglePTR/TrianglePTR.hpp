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
    void next_aa(int nSamples);
    void next_ak(int nSamples);
    void next_ka(int nSamples);
    void next_kk(int nSamples);

    // Member variables
    double mPhase;
    double mSync;
};

} // namespace TrianglePTR
