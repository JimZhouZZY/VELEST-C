#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/globals.h"

extern int isingle;
extern int iturbo;
extern int ifixsolution;
extern int icoordsystem;
extern int inrpmax;
extern int jl;
extern int nl;
extern double delta;
extern double tkj;
extern double w[IST][IEQ];
extern double dtdv[];
extern double dtdr[];
extern double h[];
extern double thk[];
extern double v[];
extern double sphase[IST][IEQ];
extern double res[IST][IEQ];

extern void maxri(int n, const double arr[], double *maxv, int *idx);

void travderiv(const char *raytype, int nl_local, int mll,
               const double v1[], const double vsq1[],
               const double rp[3][200], int nrp,
               double x2[], double y2[], double z2[], double ss[],
               double r1, double r2, int ievent, int inobs) {
    double f[100] = {0.0};
    double zmax = 0.0;
    int jndex = 0;

    if (!(isingle != 0 && iturbo == 1)) {
        for (int ii = 0; ii < nl_local; ++ii) {
            dtdv[ii] = 0.0;
        }
    }

    for (int ii = 0; ii < 3; ++ii) {
        dtdr[ii] = 0.0;
    }

    if (ifixsolution == 9) {
        return;
    }

    int ev = ievent - 1;
    int obs = inobs - 1;
    if (ev < 0 || ev >= IEQ || obs < 0 || obs >= IST) {
        return;
    }

    if (isingle != 0 && w[obs][ev] == 0.0) {
        return;
    }

    if (strcmp(raytype, "direct") == 0) {
        for (int j = 0; j < nrp - 1; ++j) {
            int jx = nrp - j - 1;
            x2[j] = rp[0][j + 1] - rp[0][j];
            y2[j] = rp[1][j + 1] - rp[1][j];
            z2[j] = rp[2][j + 1] - rp[2][j];
            ss[j] = sqrtf(x2[j] * x2[j] + y2[j] * y2[j] + z2[j] * z2[j]);
            if (!(isingle != 0 && iturbo == 1)) {
                dtdv[jx - 1] = -ss[j] / vsq1[jx - 1];
            }
        }
        dtdr[0] = -x2[0] / (v[jl] * ss[0]);
        dtdr[1] = -y2[0] / (v[jl] * ss[0]);
        dtdr[2] = -z2[0] / (v[jl] * ss[0]);
    } else if (strcmp(raytype, "refracted") == 0) {
        for (int j = 0; j < nrp; ++j) {
            z2[j] = rp[2][j];
        }

        maxri(nrp, z2, &zmax, &jndex);

        int jb = nl_local - 1;
        while (jb >= 0 && h[jb] > (zmax + 0.01f)) {
            jb--;
        }

        double dtdd = 1.0f / v[jb];
        dtdr[0] = (r1 / delta) * dtdd;
        dtdr[1] = (r2 / delta) * dtdd;
        dtdr[2] = -sqrtf(vsq1[jb] - vsq1[jl]) / (v1[jb] * v1[jl]);

        if (!(isingle != 0 && iturbo == 1)) {
            for (int j = 0; j < nl_local; ++j) {
                f[j] = 1.0f;
                if (j >= jl) f[j] = 2.0f;
                if (j > jb) f[j] = 0.0;
            }

            for (int j = 0; j < jb; ++j) {
                dtdv[jb] += thk[j] * v1[j] * f[j] /
                            (vsq1[jb] * sqrtf(vsq1[jb] - vsq1[j]));
                dtdv[j] = -thk[j] * v1[jb] * f[j] /
                          (vsq1[j] * sqrtf(vsq1[jb] - vsq1[j]));
            }

            dtdv[jl] += tkj * v1[jb] /
                        (vsq1[jl] * sqrtf(vsq1[jb] - vsq1[jl]));
            dtdv[jb] += -tkj * v1[jl] /
                        (vsq1[jb] * sqrtf(vsq1[jb] - vsq1[jl])) - delta / vsq1[jb];
        }
    } else if (strcmp(raytype, "reflected") == 0) {
        int jx = 0;
        for (int j = 0; j < nrp - 1; ++j) {
            x2[j] = rp[0][j + 1] - rp[0][j];
            y2[j] = rp[1][j + 1] - rp[1][j];
            z2[j] = rp[2][j + 1] - rp[2][j];
            ss[j] = sqrtf(x2[j] * x2[j] + y2[j] * y2[j] + z2[j] * z2[j]);

            if (j == 0) {
                dtdr[0] = -x2[0] / (v[jl] * ss[0]);
                dtdr[1] = -y2[0] / (v[jl] * ss[0]);
                dtdr[2] = -z2[0] / (v[jl] * ss[0]);
            }

            if (isingle != 0 && iturbo == 1) {
                break;
            }

            int idownward = 0;
            if (z2[j] > 0.0) idownward = 1;
            if (j == (nrp - 2) && idownward == 1) idownward = 0;

            if (idownward == 1) {
                jx = jl + j;
                dtdv[jx] = -ss[j] / vsq1[jx];
            } else {
                if (jx == mll) {
                    dtdv[mll] -= ss[j] / vsq1[mll];
                    jx--;
                } else {
                    dtdv[jx] -= ss[j] / vsq1[jx];
                    jx--;
                }
            }
        }

        if (jx != 0) {
            abort();
        }
    } else {
        abort();
    }

    if (ifixsolution == 1) {
        dtdr[2] = 0.0;
    }

    if (nrp > inrpmax) {
        fprintf(stdout, "travderiv>>> nrp greater than inrpmax\n");
        abort();
    }
}
