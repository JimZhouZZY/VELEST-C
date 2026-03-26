#include <math.h>
#include <stdio.h>
#include "../include/globals.h"

extern double avrefrres;
extern double abrefrres;
extern int nrrefrres;
extern double avotheres;
extern double abotheres;
extern int nrotheres;
extern double avreflres;
extern double abreflres;
extern int nrreflres;
extern int lmax;
extern int irfrout;
extern int irflout;
extern int iresout;
extern int icoordsystem;
extern double res[IST][IEQ];
extern double stnazires[IST][8];

void resisave(int nrp, int nrpdeep, double rp[3][200], int nobs, int i, int k1, int mll) {
    double xxx, yyy, xhyp, yhyp, xstn, ystn, azi, dist;
    int iazi;

    int obs = nobs - 1;
    int ev = i - 1;
    if (obs < 0 || obs >= IST || ev < 0 || ev >= IEQ) {
        return;
    }
    if (k1 < 0 || k1 >= IST) {
        return;
    }
    if (nrp < 2 || nrp > INRPMAX) {
        return;
    }
    if (nrpdeep < 0 || nrpdeep >= nrp - 1) {
        return;
    }

    if (rp[2][nrpdeep] == rp[2][nrpdeep + 1]) {
        avrefrres += res[obs][ev];
        abrefrres += fabs(res[obs][ev]);
        nrrefrres += 1;

        if (lmax == 10 && irfrout == 1) {
            xxx = rp[0][nrpdeep];
            if (icoordsystem == 2) {
                xxx = -xxx;
            }
            yyy = rp[1][nrpdeep];
            fprintf(stdout, "  RFR%10.3f%10.3f%10.3f        .1\n", xxx, yyy, res[obs][ev]);

            xxx = (rp[0][nrpdeep] + rp[0][nrpdeep + 1]) * 0.5f;
            if (icoordsystem == 2) {
                xxx = -xxx;
            }
            yyy = (rp[1][nrpdeep] + rp[1][nrpdeep + 1]) * 0.5f;
            fprintf(stdout, "  RFR%10.3f%10.3f%10.3f        .1\n", xxx, yyy, res[obs][ev]);

            xxx = rp[0][nrpdeep + 1];
            if (icoordsystem == 2) {
                xxx = -xxx;
            }
            yyy = rp[1][nrpdeep + 1];
            fprintf(stdout, "  RFR%10.3f%10.3f%10.3f        .1\n", xxx, yyy, res[obs][ev]);
        }
    } else {
        if (mll == 0) {
            avotheres += res[obs][ev];
            abotheres += fabs(res[obs][ev]);
            nrotheres += 1;
        } else {
            avreflres += res[obs][ev];
            abreflres += fabs(res[obs][ev]);
            nrreflres += 1;

            if (irflout == 1) {
                if (icoordsystem == 2) {
                    xxx = -rp[0][nrpdeep];
                } else {
                    xxx = rp[0][nrpdeep];
                }
                fprintf(stdout, "  RFL%10.4f%10.4f%10.4f        .1\n", xxx, rp[1][nrpdeep], res[obs][ev]);
            }
        }
    }

    xhyp = rp[0][0];
    yhyp = rp[1][0];
    xstn = rp[0][nrp - 1];
    ystn = rp[1][nrp - 1];

    azi = 57.296f * atan2f(xhyp - xstn, yhyp - ystn);
    if (azi < 0.0) {
        azi += 360.0;
    }

    iazi = 0;
    if (azi >= 0.0 && azi < 90.0) iazi = 1;
    if (azi >= 90.0 && azi < 180.0) iazi = 2;
    if (azi >= 180.0 && azi < 270.0) iazi = 3;
    if (azi >= 270.0 && azi <= 360.0) iazi = 4;

    stnazires[k1][2 * iazi - 2] += res[obs][ev];
    stnazires[k1][2 * iazi - 1] += 1.0f;

    if (iresout == 1) {
        dist = (rp[0][0] - rp[0][nrp - 1]) * (rp[0][0] - rp[0][nrp - 1])
             + (rp[1][0] - rp[1][nrp - 1]) * (rp[1][0] - rp[1][nrp - 1])
             + (rp[2][0] - rp[2][nrp - 1]) * (rp[2][0] - rp[2][nrp - 1]);
        dist = sqrtf(dist);
        fprintf(stdout, " %6.2f  %7.3f\n", dist, res[obs][ev]);
    }
}
