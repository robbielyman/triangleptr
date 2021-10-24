// PluginTrianglePTR.cpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#include "SC_PlugIn.hpp"
#include "TrianglePTR.hpp"

static InterfaceTable* ft;

namespace TrianglePTR {

TrianglePTR::TrianglePTR() {
    if (isAudioRateIn(0)) {
        if (isAudioRateIn(1)) {
            mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_aa>();
        } else {
            mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_ak>();
        }
    } else {
        if (isAudioRateIn(1)) {
            mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_ka>();
        } else {
            mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_kk>();
        }
    }
    next_kk(1);
}

void TrianglePTR::next_aa(int nSamples) {
    const float* freq   = in(0);
    const float* phase  = in(1);
    const float* sync   = in(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double lastsync = mSync;
    for (int i = 0; i < nSamples; ++i) {
        
    }
}

} // namespace TrianglePTR

PluginLoad(TrianglePTRUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<TrianglePTR::TrianglePTR>(ft, "TrianglePTR", false);
}
