// PluginTrianglePTR.cpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#include "SC_PlugIn.hpp"
#include "TrianglePTR.hpp"

static InterfaceTable* ft;

namespace TrianglePTR {

TrianglePTR::TrianglePTR() {
    if (isAudioRateIn(0)) {
        if (isAudioRateIn(1)) {
            if (isAudioRateIn(2)) {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_aaaa>();
                else                    mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_aaak>();
            }
            else {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_aaka>();
                else                    mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_aakk>();
            }
        }
        else {
            if (isAudioRateIn(2)) {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_akaa>();
                else                    mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_akak>();
            }
            else {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_akka>();
                else                    mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_akkk>();
            }
        }
    }
    else {
        if (isAudioRateIn(1)) {
            if (isAudioRateIn(2)) {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_kaaa>();
                else                    mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_kaak>();
            }
            else {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_kaka>();
                else                    mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_kakk>();
            }
        }
        else {
            if (isAudioRateIn(2)) {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_kkaa>();
                else                    mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_kkak>();
            }
            else {
                if (isAudioRateIn(3))   mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_kkka>();
                else                    mCalcFunc = make_calc_function<TrianglePTR,&TrianglePTR::next_kkkk>();
            }
        }
    }
    next_kkkk(1);
}

inline double TrianglePTR::algorithm(double p, double t0, double t2, double t3,
            double w, double a, double b, double c,
            double dc, double p1, double p2, double p3) {
   if (p < w) {
       if (p < t0) {
           // y = B*x - B*DC - 1 - P3*x^4
           return b*p - b*dc - 1.f - p3 * p * p * p * p;
       }
       else if (p < t2) {
           // y = B*x - B*DC - 1 + 2*P3*x^4 - P2*x^3 + 1.5*P1*x^2 - C*x + 0.25*C*T0
           return b*p - b*dc - 1.f + 2.f * p3 * p * p * p * p
               - p2 * p * p * p + 1.5f * p1 * p * p - c * p
               + 0.25f * c * t0;
       }
       else if (p < t3) {
           // y = B*x - B*DC - 1 - P3*x^4 + P2*x^3 - 4.5*P1*x^2 + 7*C*x - 3.75*C*T0
           return b*p - b*dc - 1.f - p3 * p * p * p * p
               + p2 * p * p * p - 4.5f * p1 * p * p
               + 7.f * c * p - 3.75f * c * t0;
       }
       else {
           // y = A*x - A*DC - 1
           return a * p - a * dc - 1.f;
       }
   }
   else {
       double pw = p - w;
       if (pw < t0) {
           // y = A*x - A*DC + 1 + P3*x^4 
           return a*pw - a*dc + 1.f + p3 * pw * pw *pw *pw;
       }
       else if (pw < t2) {
           // y = A*x - A*DC + 1 - 2*P3*x^4 + P2*x^3 - 1.5*P1*x^2 + C*x - 0.25*C*T0
           return a*pw - a*dc + 1.f - 2.f * p3 * pw * pw *pw *pw
               + p2 * pw * pw * pw - 1.5f * p1 * pw * pw + c * pw
               - 0.25f * c * t0;
       }
       else if (pw < t3) {
           // y = A*x - A*DC + 1 + P3*x^4 - P2*x^3 + 4.5*P1*x^2 - 7*c*x + 3.75*C*T0
           return a*pw - a*dc + 1.f + p3 * pw * pw * pw *pw
               - p2 * pw * pw *pw + 4.5f * p1 * pw * pw 
               - 7.f * c * pw + 3.75f * c * t0;
       }
       else {
           // y = B*x - B*DC + 1
           return b * pw - b * dc + 1.f;
       }
   }
}

void TrianglePTR::next_aaaa(int nSamples) {
    const float* freq   = in(0);
    const float* phase  = in(1);
    const float* sync   = in(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double a        = 2.f / w;
        double b        = 2.f / (w - 1);
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync[i];
        phasein   = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_aaak(int nSamples) {
    const float* freq   = in(0);
    const float* phase  = in(1);
    const float* sync   = in(2);
    const float width  = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    double a        = 2.f / w;
    double b        = 2.f / (w - 1);
    double c        = 0.5f * a * b;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync[i];
        phasein   = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_aaka(int nSamples) {
    const float* freq   = in(0);
    const float* phase  = in(1);
    const float sync   = in0(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double a        = 2.f / w;
        double b        = 2.f / (w - 1);
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync;
        phasein   = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_aakk(int nSamples) {
    const float* freq   = in(0);
    const float* phase  = in(1);
    const float sync   = in0(2);
    const float width  = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    double a        = 2.f / w;
    double b        = 2.f / (w - 1);
    double c        = 0.5f * a * b;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync;
        phasein   = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_akaa(int nSamples) {
    const float* freq   = in(0);
    const float phase   = in0(1);
    const float* sync   = in(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double a        = 2.f / w;
        double b        = 2.f / (w - 1);
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync[i];
        phasein   = phase;
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_akak(int nSamples) {
    const float* freq   = in(0);
    const float phase   = in0(1);
    const float* sync   = in(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double a        = 2.f / w;
        double b        = 2.f / (w - 1);
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync[i];
        phasein   = phase;
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_akka(int nSamples) {
    const float* freq   = in(0);
    const float phase   = in0(1);
    const float sync    = in0(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double a        = 2.f / w;
        double b        = 2.f / (w - 1);
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync;
        phasein   = phase;
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_akkk(int nSamples) {
    const float* freq   = in(0);
    const float phase   = in0(1);
    const float sync    = in0(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    double a        = 2.f / w;
    double b        = 2.f / (w - 1);
    double c        = 0.5f * a * b;
    for (int i = 0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f*step;
        double step3    = 3.f*step;
        double samples  = 1.f / step;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync;
        phasein   = phase;
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_kaaa(int nSamples) {
    const float freq    = in0(0);
    const float* phase  = in(1);
    const float* sync   = in(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double dc       = 1.5f * step;
    for (int i = 0; i < nSamples; ++i) {
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double a        = 2.f / w;
        double b        = 2.f / (w - 1);
        double c        = 0.5f * a * b;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync[i];
        phasein   = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_kaak(int nSamples) {
    const float freq    = in0(0);
    const float* phase  = in(1);
    const float* sync   = in(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    double a        = 2.f / w;
    double b        = 2.f / (w - 1);
    double c        = 0.5f * a * b;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double dc       = 1.5f * step;
    double p1       = c * samples;
    double p2       = p1 * samples;
    double p3       = p2 * samples / 12.f;
    for (int i = 0; i < nSamples; ++i) {
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync[i];
        phasein   = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_kaka(int nSamples) {
    const float freq    = in0(0);
    const float* phase  = in(1);
    const float sync    = in0(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double dc       = 1.5f * step;
    for (int i = 0; i < nSamples; ++i) {
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double a        = 2.f / w;
        double b        = 2.f / (w - 1);
        double c        = 0.5f * a * b;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync;
        phasein   = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_kakk(int nSamples) {
    const float freq    = in0(0);
    const float* phase  = in(1);
    const float sync    = in0(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    double a        = 2.f / w;
    double b        = 2.f / (w - 1);
    double c        = 0.5f * a * b;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double dc       = 1.5f * step;
    double p1       = c * samples;
    double p2       = p1 * samples;
    double p3       = p2 * samples / 12.f;
    for (int i = 0; i < nSamples; ++i) {
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync;
        phasein   = phase[i];
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_kkaa(int nSamples) {
    const float freq    = in0(0);
    const float phase   = in0(1);
    const float* sync   = in(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double dc       = 1.5f * step;
    for (int i = 0; i < nSamples; ++i) {
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double a        = 2.f / w;
        double b        = 2.f / (w - 1);
        double c        = 0.5f * a * b;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync[i];
        phasein   = phase;
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_kkak(int nSamples) {
    const float freq    = in0(0);
    const float phase   = in0(1);
    const float* sync   = in(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double a        = 2.f / w;
    double b        = 2.f / (w - 1);
    double c        = 0.5f * a * b;
    double dc       = 1.5f * step;
    double p1       = c * samples;
    double p2       = p1 * samples;
    double p3       = p2 * samples / 12.f;
    for (int i = 0; i < nSamples; ++i) {
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync[i];
        phasein   = phase;
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_kkka(int nSamples) {
    const float freq    = in0(0);
    const float phase   = in0(1);
    const float sync    = in0(2);
    const float* width  = in(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double dc       = 1.5f * step;
    for (int i = 0; i < nSamples; ++i) {
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        double a        = 2.f / w;
        double b        = 2.f / (w - 1);
        double c        = 0.5f * a * b;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync;
        phasein   = phase;
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}

void TrianglePTR::next_kkkk(int nSamples) {
    const float freq    = in0(0);
    const float phase   = in0(1);
    const float sync    = in0(2);
    const float width   = in0(3);
    float* outbuf   = out(0);
    float* syncOut  = out(1);

    double pos      = mPhase;
    double rateinv  = 1 / sampleRate();
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    double a        = 2.f / w;
    double b        = 2.f / (w - 1);
    double c        = 0.5f * a * b;
    double step     = freq * rateinv;
    double step2    = 2.f*step;
    double step3    = 3.f*step;
    double samples  = 1.f / step;
    double dc       = 1.5f * step;
    double p1       = c * samples;
    double p2       = p1 * samples;
    double p3       = p2 * samples / 12.f;
    for (int i = 0; i < nSamples; ++i) {
        double out1      = 0.f;
        double out2      = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
        }
        else {
            pos = pos + (phase - phasein) + step;
            while (pos >= 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
        }
        lastsync  = sync;
        phasein   = phase;
        if (pos < step) out2 = 1.f;
        out1 = algorithm(pos, step, step2, step3, w, a, b, c, dc, p1, p2, p3); 

        syncOut[i] = out2;
        outbuf[i] = out1;
    }
    mPhase      = pos;
    mSync       = lastsync;
    mPhaseIn    = phasein;
}
} // namespace TrianglePTR

PluginLoad(TrianglePTRUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<TrianglePTR::TrianglePTR>(ft, "TrianglePTR", false);
}
