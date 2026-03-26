#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include "../include/globals.h"

extern int legs;
extern int neqs;
extern int nshot;
extern int nltot;
extern int nmod;
extern int nsp;
extern int nsta;
extern int ksta;
extern int icoordsystem;
extern bool single_turbo;
extern int isingle;
extern int ireflector;
extern char reflchar;
extern FILE *fm_ptr;

extern int isconstrain[];
extern int iconstrain[];
extern int map1[];
extern double e[5][IEQ];
extern double b[];
extern double rms[IEQ];
extern double avres[IEQ];
extern double scale[];
extern double vp[ITOTMODELS][INLTOT];
extern double hp[ITOTMODELS][INLTOT];
extern int nplay[INLTOT];
extern int laysum[INLTOT];
extern double ptcor[IST];
extern double stcor[IST];
extern char stn[IST][5];

void nittoutput(double damp) {
    double avdt = 0.0, avdx = 0.0, avdy = 0.0, avdz = 0.0;
    double aavdt = 0.0, aavdx = 0.0, aavdy = 0.0, aavdz = 0.0;
    double cc[IST] = {0.0};

    FILE *out = fm_ptr ? fm_ptr : stdout;

    fprintf(out, "\n");
    fprintf(out, "  eq       ot     x      y      z      rms   avres   dot     dx     dy     dz\n");

    for (int i = 0; i < legs; ++i) {
        if (isingle != 0) {
            if (isconstrain[0] == 1) fprintf(out, " *** nitt<2 --> depth-adjustment := 0.0 for event %d\n", i + 1);
            if (isconstrain[1] == 1) fprintf(out, " *** igap>250 --> depth-adjustment := 0.0 for event %d\n", i + 1);
            if (isconstrain[2] == 1) fprintf(out, " *** depth-adjustment constrained for event %d\n", i + 1);
        }
        if (iconstrain[i] == 1) {
            fprintf(out, " ***** depth constrained for event %5d\n", i + 1);
        }

        int j1 = 4 * i;
        if (i >= neqs) {
            j1 = 3 * neqs + i;
        }

        if (icoordsystem == 2) {
                fprintf(out, " %4d %7.2f%7.2f%7.2f%7.2f%7.2f%7.2f %7.3f%7.3f%7.3f%7.3f\n",
                    i + 1, e[0][i], -e[1][i], e[2][i], e[3][i], rms[i], avres[i],
                    b[j1], -b[j1 + 1], b[j1 + 2], b[j1 + 3]);
        } else {
                fprintf(out, " %4d %7.2f%7.2f%7.2f%7.2f%7.2f%7.2f %7.3f%7.3f%7.3f%7.3f\n",
                    i + 1, e[0][i], e[1][i], e[2][i], e[3][i], rms[i], avres[i],
                    b[j1], b[j1 + 1], b[j1 + 2], b[j1 + 3]);
        }

        avdt += b[j1];
        avdx += b[j1 + 1];
        avdy += b[j1 + 2];
        avdz += b[j1 + 3];
        aavdt += fabs(b[j1]);
        aavdx += fabs(b[j1 + 1]);
        aavdy += fabs(b[j1 + 2]);
        aavdz += fabs(b[j1 + 3]);
    }

    if (icoordsystem == 2) avdx = -avdx;

    avdt /= (double)legs;
    avdx /= (double)legs;
    avdy /= (double)legs;
    avdz /= (double)legs;
    aavdt /= (double)legs;
    aavdx /= (double)legs;
    aavdy /= (double)legs;
    aavdz /= (double)legs;

        fprintf(out, "\n A V E R A G E   of ADJUSTMENTS :                    %7.3f%7.3f%7.3f%7.3f\n",
            avdt, avdx, avdy, avdz);
        fprintf(out, " A V E R A G E   of ABSOLUTE ADJUSTMENTS :           %7.3f%7.3f%7.3f%7.3f\n\n",
            aavdt, aavdx, aavdy, aavdz);

    if (damp != 1.0f) {
        fprintf(out, "\n Step length damping of %7.5f was applied.\n\n", damp);
    } else {
        fprintf(out, "\nNO step length damping applied\n\n");
    }

    if (scale[4] != 0.0) {
        int j1 = 4 * neqs + nshot + nltot;
        int ksta1 = (nsp == 2) ? (ksta / 2) : ksta;

        fprintf(out, "\n Adjusted station corrections:\n");
        fprintf(out, "  stn  ptcor  dpcor\n");

        for (int j = 0; j < nsta; ++j) {
            cc[j] = 0.0;
            if (map1[j] > 0 && map1[j] <= ksta1) {
                int kk1 = j1 + map1[j] - 1;
                cc[j] = b[kk1];
            }
            fprintf(out, "  %4s %7.3f %7.3f\n", stn[j], ptcor[j], cc[j]);
        }

        if (nsp == 2) {
            int ksta2 = ksta - ksta1;
            j1 = 4 * neqs + nshot + nltot + ksta1;
            fprintf(out, " Adjusted station corrections:\n");
            fprintf(out, "  stn  stcor  dscor\n");
            for (int j = 0; j < nsta; ++j) {
                cc[j] = 0.0;
                if (map1[j] > 0 && map1[j] <= ksta2) {
                    int kk1 = j1 + map1[j] - 1;
                    cc[j] = b[kk1];
                }
                fprintf(out, "  %4s %7.3f %7.3f\n", stn[j], stcor[j], cc[j]);
            }
        }
    }

    fprintf(out, "\n");
}
