// PluginTrianglePTR.hpp
// Rylee Lyman (ryleealanza@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"

namespace TrianglePTR {

class TrianglePTR : public SCUnit {
public:
    TrianglePTR();

    // Destructor
    // ~TrianglePTR();

private:
    // Calc function
    void next(int nSamples);

    // Member variables
};

} // namespace TrianglePTR
