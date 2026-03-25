#include <math.h>
#include <complex.h>
#include <string.h>

static double complex b6(double om) {
    const double zpi = 6.28319f;
    const double fg1 = 12.0f;
    const double omg1 = om / (zpi * fg1);
    const double complex j = I;

    double complex gbn1 = 1.0f + 1.2217f * j * omg1 - 0.3887f * omg1 * omg1;
    double complex gbn2 = 1.0f + 0.9686f * j * omg1 - 0.3505f * omg1 * omg1;
    double complex gbn3 = 1.0f + 0.5131f * j * omg1 - 0.2756f * omg1 * omg1;
    return 1.0f / (gbn1 * gbn2 * gbn3);
}

static double complex hp1(double om) {
    const double zpi = 6.28319f;
    const double fg = 0.1f;
    const double tau = 1.0f / (zpi * fg);
    const double complex j = I;
    double complex h = j * om * tau;
    return h / (1.0f + h);
}

static double complex sts373(double om) {
    const double zpi = 6.28319f;
    const double fg1 = 5.0f;
    const double omg1 = zpi * fg1;
    const double complex j = I;
    double complex o = om / omg1;
    double complex gbn1 = 1.0f + 1.3397f * j * o - 0.4889f * o * o;
    double complex gbn2 = 1.0f + 0.7743f * j * o - 0.3890f * o * o;
    return 1.0f / (gbn1 * gbn2);
}

static double complex dev(double om) {
    const double zpi = 6.28319f;
    const double fg = 15.0f;
    const double tau = 1.0f / (zpi * fg);
    const double complex j = I;
    return 1.6f * cpowf(1.0f / (1.0f + j * om * tau), 2.0f);
}

void ufwabo(int isecpendel, double seismkonst, double seismdamp, double voltgain,
            double devampl, double period, int iampltype, double *ampl, const char *ifilter) {
    const double zpi = 6.28319f;
    const double complex j = I;

    double t02 = 0.0f;
    double t01 = 0.0f;
    if (isecpendel == 1) {
        t02 = 11.9f;
        t01 = 3.7f;
    }
    if (isecpendel == 2) {
        t02 = 12.7f;
        t01 = 112.0f;
    }

    double ts = (double)isecpendel;
    double s = seismkonst;
    double hs = seismdamp;

    double om = zpi / period;
    double complex jom = j * om;
    double om2 = om * om;

    double omg2 = zpi / t02;
    double tau2 = 1.0f / omg2;
    double complex g2 = jom * tau2;
    g2 = cpowf(g2 / (1.0f + g2), 2.0f);

    double omg1 = zpi / t01;
    double tau1 = 1.0f / omg1;
    double complex g1 = jom * tau1;
    g1 = g1 / (1.0f + g1);

    double om0 = zpi / ts;
    double complex gseismb = om2 * jom * s / (om0 * om0 + 2.0f * jom * om0 * hs - om2);

    double complex G = g1 * g2 * gseismb / hp1(om);

    if (strcmp(ifilter, "DE") == 0) {
        G = G * sts373(om) * dev(om);
    }

    if (strcmp(ifilter, "AD") == 0) {
        G = G * hp1(om) * b6(om);
        double gain = 16.0f;
        G = G / 10.0f;
        G = G * gain;
    }

    if (iampltype == 2) {
        G = G / 1000.0f;
    } else if (iampltype == 1) {
        double twa = 0.8f;
        double o0w = zpi / twa;
        double hw = 0.78f;
        double vwa = 2800.0f;
        G = G / (om2 * vwa / (o0w * o0w + 2.0f * jom * o0w * hw - om2));
    }

    double amplresp = cabs(G);
    *ampl = devampl / (voltgain * amplresp);
}
