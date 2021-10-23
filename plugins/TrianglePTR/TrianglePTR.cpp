// PluginTrianglePTR.cpp
// Rylee Lyman (ryleealanza@gmail.com)

#include "SC_PlugIn.hpp"
#include "TrianglePTR.hpp"

static InterfaceTable* ft;

namespace TrianglePTR {

TrianglePTR::TrianglePTR() {
    mCalcFunc = make_calc_function<TrianglePTR, &TrianglePTR::next>();
    next(1);
}

void TrianglePTR::next(int nSamples) {
    const float* input = in(0);
    const float* gain = in(1);
    float* outbuf = out(0);

    // simple gain function
    for (int i = 0; i < nSamples; ++i) {
        outbuf[i] = input[i] * gain[i];
    }
}

} // namespace TrianglePTR

PluginLoad(TrianglePTRUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<TrianglePTR::TrianglePTR>(ft, "TrianglePTR", false);
}
