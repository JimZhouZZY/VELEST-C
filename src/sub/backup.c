#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define IEQ 658

typedef void (*chtop_backup_fn_t)(double x, double y, double *zmin,
                                  const char *topo1file, const char *topo2file);

void backup_model(
    int legs,
    int neqs,
    int nshot,
    int nitt,
    int invertratio,
    int nmod,
    int nltot,
    int nsta,
    int nsp,
    int nsinv,
    int ksta,
    int itopo,
    int ifixsolution,
    int lowveloclay,
    int ireflector,
    char reflchar,
    const int *nplay,
    const int *laysum,
    const int *map1,
    double *b,
    double e[][IEQ],
    double *vp,
    const double *hp,
    double *ptcor,
    double *stcor,
    bool single_turbo,
    FILE *logfp,
    FILE *warnfp,
    const char (*stn)[5],
    const char *topo1file,
    const char *topo2file,
    chtop_backup_fn_t chtop_fn) {
    int i;
    int j;
    int k;
    int jjj = -1;
    double *cc;

    if (legs <= 0 || neqs < 0 || nsta < 0 || nmod < 0 || nltot < 0 ||
        invertratio == 0 || nplay == NULL || laysum == NULL || map1 == NULL ||
        b == NULL || e == NULL || vp == NULL || hp == NULL || ptcor == NULL ||
        stcor == NULL) {
        return;
    }

    if (!single_turbo && logfp == NULL) {
        logfp = stdout;
    }
    if (warnfp == NULL) {
        warnfp = stderr;
    }

    cc = (double *)calloc((size_t)nsta, sizeof(double));
    if (cc == NULL) {
        return;
    }

    for (i = 0; i < legs; ++i) {
        int n = (i < neqs) ? 4 : 1;
        for (k = 0; k < n; ++k) {
            ++jjj;
            b[jjj] /= 2.0f;

            if (k == 3) {
                double zzz = 0.0f;
                e[k][i] -= b[jjj];

                if (itopo > 0 && e[k][i] < 0.0f && chtop_fn != NULL) {
                    chtop_fn(-e[1][i], e[2][i], &zzz, topo1file,
                             topo2file);
                    if (e[k][i] < zzz) {
                        e[k][i] = zzz;
                    }
                    if (ifixsolution > 0) {
                        e[k][i] = zzz;
                    }
                }
            } else {
                e[k][i] -= b[jjj];
            }
        }
    }

    if ((nitt % invertratio) == 0) {
        for (i = 0; i < nmod; ++i) {
            int nl = nplay[i];
            int j1 = 4 * neqs + nshot + laysum[i];

            if (!single_turbo && logfp != NULL) {
                fprintf(logfp, " Velocity readjustments:\n");
                fprintf(logfp, " Velocity model%4d\n", i + 1);
            }

            for (k = 0; k < nl; ++k) {
                int bidx = j1 + k;
                char reflch = (k == (ireflector - 1)) ? reflchar : ' ';

                b[bidx] /= 2.0f;
                vp[i * nltot + k] -= b[bidx];

                if (lowveloclay == 0 && k > 0 &&
                    vp[i * nltot + k] < vp[i * nltot + (k - 1)]) {
                    if (!single_turbo && logfp != NULL) {
                        fprintf(logfp,
                                " WARNING: Tried to introduce a low-velocity-layer! (Layer %2d)\n",
                                k + 1);
                        fprintf(logfp,
                                " Setting DVP from %5.2f to 0.0 and VP to vp(layer_above)+0.001\n",
                                b[bidx]);
                    }
                    fprintf(warnfp,
                            " WARNING: Tried to introduce a low-velocity-layer! (Layer %2d)\n",
                            k + 1);
                    fprintf(warnfp,
                            " Setting DVP from %5.2f to 0.0 and VP to vp(layer_above)+0.001\n",
                            b[bidx]);
                    b[bidx] = 0.0f;
                    vp[i * nltot + k] = vp[i * nltot + (k - 1)] + 0.001f;
                }

                if (!single_turbo && logfp != NULL) {
                    fprintf(logfp, " %7.3f %7.3f %7.3f   %c\n", vp[i * nltot + k],
                            b[bidx], hp[i * nltot + k], reflch);
                }
            }
        }
    }

    if (nsinv != 0 && (nitt % invertratio) == 0) {
        int k1 = 4 * neqs + nshot + nltot;
        int ksta1 = (nsp == 2) ? (ksta / 2) : ksta;
        int k2 = k1 + ksta1;

        for (j = k1; j < k2; ++j) {
            b[j] /= 2.0f;
        }

        for (j = 0; j < nsta; ++j) {
            int m = map1[j];
            cc[j] = 0.0f;
            if (m < 0 || m >= ksta1) {
                continue;
            }
            cc[j] = b[k1 + m];
            ptcor[j] -= cc[j];
        }

        if (!single_turbo && logfp != NULL) {
            for (j = 0; j < nsta; ++j) {
                fprintf(logfp, "  %.4s %7.3f %7.3f", stn != NULL ? stn[j] : "    ",
                        ptcor[j], cc[j]);
                if ((j + 1) % 5 == 0 || j + 1 == nsta) {
                    fprintf(logfp, "\n");
                }
            }
            fprintf(logfp, "\n");
            fprintf(logfp, " Half adjustments made\n");
        }

        if (nsp == 2) {
            int k1s = 4 * neqs + nshot + nltot + ksta1;
            int ksta2 = ksta - ksta1;
            int k2s = k1s + ksta2;

            for (j = k1s; j < k2s; ++j) {
                b[j] /= 2.0f;
            }

            if (!single_turbo && logfp != NULL) {
                fprintf(logfp, " S correction readjustments:\n");
            }

            for (j = 0; j < nsta; ++j) {
                int m = map1[j];
                cc[j] = 0.0f;
                if (m < 0 || m >= ksta2) {
                    continue;
                }
                cc[j] = b[k1s + m];
                stcor[j] -= cc[j];
            }

            if (!single_turbo && logfp != NULL) {
                for (j = 0; j < nsta; ++j) {
                    fprintf(logfp, "  %.4s %7.3f %7.3f", stn != NULL ? stn[j] : "    ",
                            stcor[j], cc[j]);
                    if ((j + 1) % 5 == 0 || j + 1 == nsta) {
                        fprintf(logfp, "\n");
                    }
                }
            }
        }
    }

    if (!single_turbo && logfp != NULL) {
        fprintf(logfp, "\n");
        fprintf(logfp, " Half adjustments made\n");
    }

    free(cc);
}

extern int legs, neqs, nshot, nitt, invertratio, nmod, nltot, nsta, nsp, nsinv, ksta;
extern int itopo, ifixsolution, lowveloclay, ireflector;
extern char reflchar;
extern int nplay[];
extern int laysum[];
extern int map1[];
extern double b[];
extern double e[][658];
extern double vp[][100];
extern double hp[][100];
extern double ptcor[];
extern double stcor[];
extern bool single_turbo;
extern char stn[][5];
extern char topo1file[81], topo2file[81];
extern void chtop(double x, double y, double *zmin, const char *topo1file, const char *topo2file);

void backup(void) {
    backup_model(
        legs,
        neqs,
        nshot,
        nitt,
        invertratio,
        nmod,
        nltot,
        nsta,
        nsp,
        nsinv,
        ksta,
        itopo,
        ifixsolution,
        lowveloclay,
        ireflector,
        reflchar,
        nplay,
        laysum,
        map1,
        b,
        e,
        &vp[0][0],
        &hp[0][0],
        ptcor,
        stcor,
        single_turbo != 0,
        stdout,
        stderr,
        (const char (*)[5])stn,
        topo1file,
        topo2file,
        chtop);
}
