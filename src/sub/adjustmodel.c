#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include "../include/globals.h"

extern int legs, neqs, nshot, nitt, isingle, itopo, ifixsolution;
extern int nsta, nsp, ksta, nltot, nmod, lowveloclay;
extern int knobs[IEQ], igap[IEQ], isconstrain[3], iconstrain[IEQ];
extern int map1[IST], iyr[IEQ], imo[IEQ], iday[IEQ], ihr[IEQ], imin[IEQ];
extern int nplay[INLTOT], laysum[INLTOT];

extern double zadj, zmin;
extern double e[5][IEQ], b[INVA * 2], pt[IST][IEQ], d[IST][3][IEQ];
extern double delta, scale[7], vp[ITOTMODELS][INLTOT], ptcor[IST], stcor[IST];
extern bool single_turbo;
extern FILE *fm_ptr;
extern char topo1file[81], topo2file[81];

extern int juliam(int iyr, int imo, int idy, int ihr, int imn);
extern void timeclear(int *iyr, int *imo, int *iday, int *ihr, int *imin,
                      double *sec, int *itime);
extern void chtop(double xx, double yy, double *zk, const char *topo1file,
                  const char *topo2file);

void adjustmodel(double damp) {
    (void)damp;

    int jj = -1;

    for (int i = 0; i < legs; ++i) {
        int n = (i >= neqs) ? 1 : 4;

        int year = (iyr[i] < 100) ? (iyr[i] + 1900) : iyr[i];
        int iminold = juliam(year, imo[i], iday[i], ihr[i], imin[i]);
        // Debug: print the shape of e
        fprintf(stderr, "DEBUG: e shape = [%zu][%zu]\n", sizeof(e)/sizeof(e[0]), sizeof(e[0])/sizeof(e[0][0]));
        fprintf(stderr, "DEBUG: i=%d, e[0]=%f, e[1]=%f, e[2]=%f, e[3]=%f, e[4]=%f\n", i, e[0][i], e[1][i], e[2][i], e[3][i], e[4][i]);
        for (int j = 0; j < n; ++j) {
            ++jj;
            e[j][i] = e[j][i] + b[jj];
        }

        if (e[0][i] < 0.0) {
            int itime = iminold;
            timeclear(&iyr[i], &imo[i], &iday[i], &ihr[i], &imin[i], &e[0][i], &itime);
            for (int j = 0; j < knobs[i] && j < IST; ++j) {
                pt[j][i] = pt[j][i] + (double)(iminold - itime) * 60.0;
            }
        }

        for (int j = 0; j < 3; ++j) {
            isconstrain[j] = 0;
        }
        iconstrain[i] = 0;

        if (isingle != 0) {
            if (nitt < 2) {
                isconstrain[0] = 1;
                e[3][i] = e[3][i] - b[3];
                b[3] = 0.0;
            }

            if (igap[0] > 250) {
                double dmin = 999.9f;
                for (int k = 0; k < knobs[0] && k < IST; ++k) {
                    delta = sqrtf((e[1][0] - d[k][0][0]) * (e[1][0] - d[k][0][0]) +
                                  (e[2][0] - d[k][1][0]) * (e[2][0] - d[k][1][0]));
                    if (delta < dmin) {
                        dmin = delta;
                    }
                }
                if (dmin > 15.0f) {
                    isconstrain[1] = 1;
                    if (iconstrain[i] == 1) {
                        iconstrain[i] = 3;
                    }
                    e[3][i] = e[3][i] - b[3];
                    b[3] = 0.0;
                }
            }

            if (fabs(b[3]) > zadj) {
                isconstrain[2] = 1;
                if (b[3] > 0.0) {
                    e[3][i] = e[3][i] - b[3] + zadj;
                    b[3] = zadj;
                } else {
                    e[3][i] = e[3][i] - b[3] - zadj;
                    b[3] = -zadj;
                }
            }
        } else {
            if (i < neqs && fabs(b[jj]) > zadj) {
                iconstrain[i] = 1;
                if (b[jj] > 0.0) {
                    e[3][i] = e[3][i] - b[jj] + zadj;
                    b[jj] = zadj;
                } else {
                    e[3][i] = e[3][i] - b[jj] - zadj;
                    b[jj] = -zadj;
                }
            }
        }

        if (itopo > 0) {
            if (e[3][i] < 0.0) {
                chtop(-e[1][i], e[2][i], &zmin, topo1file, topo2file);
            } else {
                zmin = 0.0;
            }
        }

        if (e[3][i] < zmin || (ifixsolution > 0 && e[3][i] <= 0.0)) {
            b[jj] = b[jj] - (e[3][i] - zmin);
            e[3][i] = zmin;
            iconstrain[i] = 1;
        }
    }

    if (scale[5] != 0.0) {
        if (!single_turbo) {
            FILE *out = fm_ptr ? fm_ptr : stdout;
            fprintf(out, "\n doing velocity adjustments now...\n");
        }

        int j1 = 4 * neqs + nshot; /* 0-based start */

        for (int i = 0; i < nmod; ++i) {
            int k2 = laysum[i];                 /* Fortran 1-based */
            int j11 = j1 + k2 - 1;              /* convert */
            int j22 = j11 + nplay[i] - 1;
            int kj = 1;

            for (int jjj = j11; jjj <= j22; ++jjj) {
                vp[i][kj - 1] = vp[i][kj - 1] + b[jjj];

                if (lowveloclay == 0 && kj > 1) {
                    if (vp[i][kj - 1] < vp[i][kj - 2]) {
                        if (!single_turbo) {
                            FILE *out = fm_ptr ? fm_ptr : stdout;
                            fprintf(out, " WARNING: Tried to introduce a low-velocity-layer! (Layer %2d)\n", kj);
                            fprintf(out, " Setting DVP from %5.2f to 0.0 and VP to vp(layer_above)+0.001\n", b[jjj]);
                        }
                        fprintf(stderr, " WARNING: Tried to introduce a low-velocity-layer! (Layer %2d)\n", kj);
                        fprintf(stderr, " Setting DVP from %5.2f to 0.0 and VP to vp(layer_above)+0.001\n", b[jjj]);
                        b[jjj] = 0.0;
                        vp[i][kj - 1] = vp[i][kj - 2] + 0.001f;
                    }
                }
                kj += 1;
            }
        }
    }

    if (scale[4] != 0.0) {
        if (!single_turbo) {
            FILE *out = fm_ptr ? fm_ptr : stdout;
            fprintf(out, "\n doing station-correction adjustments...\n\n");
        }

        int k1 = 4 * neqs + nshot + nltot + 1; /* Fortran 1-based */
        int ksta1 = ksta;
        if (nsp == 2) {
            ksta1 = ksta / 2;
        }

        for (int j = 0; j < nsta; ++j) {
            if (map1[j] == 0) continue;
            if (map1[j] > ksta1) continue;
            int kk1 = k1 - 1 + map1[j];
            int idx = kk1 - 1;
            ptcor[j] = ptcor[j] + b[idx];
        }

        if (nsp == 2) {
            k1 = 4 * neqs + nshot + nltot + 1 + ksta1;
            int ksta2 = ksta - ksta1;
            for (int j = 0; j < nsta; ++j) {
                if (map1[j] == 0) continue;
                if (map1[j] > ksta2) continue;
                int kk1 = k1 - 1 + map1[j];
                int idx = kk1 - 1;
                stcor[j] = stcor[j] + b[idx];
            }
        }
    }
}
