#include <math.h>
#include <string.h>

extern void ufwabo(int isecpendel, double seismkonst, double seismdamp, double voltgain,
                   double devampl, double period, int iampltype, double *ampl,
                   const char *ifilter);
extern void mpdr2(int isecpendel, double seismkonst, double seismdamp,
                  double voltgain, double pdrampl, double period, int iampltyp,
                  double *ampl);

void muk(double epdist, double depth, const char *ifilter, int isecpendel,
         double seismkonst, double seismdamp, double voltgain,
         double stacor, double devampl, double period, double *xmag) {
    static const double rdelt[12] = {4.0f, 6.0f, 13.0f, 16.0f, 28.0f, 87.0f,
                                    114.0f, 120.0f, 134.0f, 141.0f, 146.0f, 180.0f};
    static const double sigar[12] = {6.2f, 7.0f, 7.0f, 5.8f, 6.6f, 7.0f,
                                    7.5f, 7.5f, 6.9f, 7.1f, 6.9f, 6.9f};

    if (isecpendel <= 0 || seismkonst <= 0.0f || seismdamp <= 0.0f ||
        voltgain <= 0.0f || devampl <= 0.0f || period <= 0.0f) {
        *xmag = -13.0f;
        return;
    }

    *xmag = -13.0f;

    double epdistkm;
    if (epdist < 0.0f) {
        epdistkm = -epdist / 360.0f * 40030.0f;
    } else {
        epdistkm = sqrtf(epdist * epdist + depth * depth);
    }

    int iampltype = (epdistkm <= 2200.0f) ? 1 : 2;

    double ampl = 0.0f;
    if (strcmp(ifilter, "DE") == 0 || strcmp(ifilter, "AD") == 0) {
        ufwabo(isecpendel, seismkonst, seismdamp, voltgain,
               devampl, period, iampltype, &ampl, ifilter);
    } else if (strcmp(ifilter, "PD") == 0) {
        mpdr2(isecpendel, seismkonst, seismdamp, voltgain,
              devampl, period, iampltype, &ampl);
    }

    if (ampl <= 1.0e-10f) {
        *xmag = -13.0f;
        return;
    }

    double waampl = 0.0f;
    double boampl = 0.0f;

    if (iampltype == 1) {
        waampl = ampl / 2.0f;
    } else {
        boampl = ampl;
    }

    if (epdistkm > 2200.0f) {
        double delta = epdistkm / 40030.0f * 360.0f;
        double sigma = sigar[11];

        if (delta > 180.0f) {
            delta = 180.0f;
            sigma = sigar[11];
        } else {
            int idx = 0;
            while (idx < 12 && delta > rdelt[idx]) {
                idx++;
            }
            if (idx <= 0) {
                sigma = sigar[0];
            } else if (idx >= 12) {
                sigma = sigar[11];
            } else {
                double a = (sigar[idx] - sigar[idx - 1]) / (rdelt[idx] - rdelt[idx - 1]);
                double b = sigar[idx] - a * rdelt[idx];
                sigma = a * delta + b;
            }
        }

        *xmag = log10f(boampl / period) + sigma + stacor;
        return;
    }

    if (epdistkm > 0.0f && epdistkm <= 60.0f) {
        *xmag = log10f(waampl) + 0.018f * epdistkm + 1.77f + 0.40f;
    }
    if (epdistkm > 60.0f && epdistkm <= 700.0f) {
        *xmag = log10f(waampl) + 0.0038f * epdistkm + 2.62f + 0.40f;
    }
    if (epdistkm > 1100.0f && epdistkm <= 1700.0f) {
        *xmag = log10f(waampl) + 0.0029f * epdistkm + 3.40f + 0.40f - 2.0f;
    }

    if (*xmag == -13.0f) {
        return;
    }

    *xmag += stacor;
}
