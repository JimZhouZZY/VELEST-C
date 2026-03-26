#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../include/globals.h"

extern int isingle;
extern int iturbo;
extern int legs;
extern int nsp;
extern int nsta;
extern int istaout;
extern int ielev[IST];
extern int model[ITOTMODELS * IST];
extern int knobs[IEQ];
extern int map1[IST];
extern int kpwt[IST][IEQ];
extern int iyr[IEQ], imo[IEQ], iday[IEQ], ihr[IEQ], imin[IEQ];
extern double e[5][IEQ];
extern double d[IST][3][IEQ];
extern double res[IST][IEQ];
extern double tctime[IST][IEQ];
extern double w[IST][IEQ];
extern double sphase[IST][IEQ];
extern double ptcor[IST], stcor[IST];
extern double xla[IST], xlo[IST];
extern double vpvs;
extern double vp[ITOTMODELS][INLTOT];
extern double hp[ITOTMODELS][INLTOT];
extern double vdamp[ITOTMODELS][INLTOT];
extern int nmod;
extern int nplay[INLTOT];
extern char smn[IST][IEQ][5];
extern char stn[IST][5];
extern char fm[81];
extern char stafile[81];

void finalstaresi(void) {
    static const char phz[3] = {'S', 'P', 'm'};
    double del[IST];
    double aa[IST], bb[IST], dd[IST], ee[IST];
    double aas[IST], bbs[IST], dds[IST], ees[IST];
    int icc[IST], iccs[IST];
    char phzz[IST];

    if (isingle == 0) {
        if (iturbo == 0) {
            fprintf(stdout, "\n\n\n\n");
        } else {
            fprintf(stdout, "TURBO-option is set; residuals are NOT printed for each event!\n");
        }
    }

    if (isingle == 0 && iturbo == 0) {
        for (int i = 0; i < legs; ++i) {
            fprintf(stdout, "\n");
            int nin = imin[i];
            int h = ihr[i];
            if (nin < 0) {
                nin += 60;
                h -= 1;
            }
            fprintf(stdout,
                    " Station residuals for event=%4d   %02d%02d%02d %02d%02d %5.2f\n",
                    (isingle == 0) ? (i + 1) : isingle,
                    iyr[i], imo[i], iday[i], h, nin, e[0][i]);
            fprintf(stdout, " sta ph wt  res   ttime delta\n");
            for (int j = 0; j < knobs[i]; ++j) {
                del[j] = sqrtf((e[1][i] - d[j][0][i]) * (e[1][i] - d[j][0][i]) +
                               (e[2][i] - d[j][1][i]) * (e[2][i] - d[j][1][i]));
                phzz[j] = 'p';
                if (sphase[j][i] == 1.0f) phzz[j] = 's';
                if (sphase[j][i] == -1.0f) phzz[j] = 'm';
                fprintf(stdout, " %4s %c %2d %7.3f %6.2f %6.2f\n",
                        smn[j][i], phzz[j], kpwt[j][i], res[j][i], tctime[j][i], del[j]);
            }
        }
    }

    fprintf(stdout, "\n\n station statistics, remember nsp was set to:%2d\n", nsp);
    fprintf(stdout, "\n sta phase nobs avres  avwres    std    wsum    delay\n\n");

    for (int m = 0; m < nsta; ++m) {
        aa[m] = bb[m] = dd[m] = ee[m] = 0.0;
        aas[m] = bbs[m] = dds[m] = ees[m] = 0.0;
        icc[m] = iccs[m] = 0;
    }

    for (int m = 0; m < nsta; ++m) {
        for (int i = 0; i < legs; ++i) {
            int k = knobs[i];
            for (int j = 0; j < k; ++j) {
                if (strcmp(stn[m], smn[j][i]) == 0 && (sphase[j][i] == 0.0 || sphase[j][i] == -1.0f)) {
                    aa[m] += res[j][i] * w[j][i];
                    bb[m] += res[j][i] * res[j][i] * w[j][i] * w[j][i];
                    ee[m] += res[j][i];
                    dd[m] += w[j][i];
                    icc[m] += 1;
                    break;
                }
            }
            if (nsp == 1) continue;
            for (int j = 0; j < k; ++j) {
                if (strcmp(stn[m], smn[j][i]) == 0 && (sphase[j][i] == 1.0f || sphase[j][i] == 2.0f)) {
                    aas[m] += res[j][i] * w[j][i];
                    bbs[m] += res[j][i] * res[j][i] * w[j][i] * w[j][i];
                    ees[m] += res[j][i];
                    dds[m] += w[j][i];
                    iccs[m] += 1;
                    break;
                }
            }
        }
    }

    FILE *stfp = NULL;
    if (istaout > 0) {
        stfp = fopen(stafile, "w");
        if (stfp != NULL) fprintf(stfp, "%s\n", fm);
    }

    for (int m = 0; m < nsta; ++m) {
        if (stfp != NULL) {
            char cns = 'N';
            char cew = 'W';
            double lat = xla[m];
            double lon = xlo[m];
            if (lat < 0.0) { cns = 'S'; lat = -lat; }
            if (lon < 0.0) { cew = 'E'; lon = -lon; }
            fprintf(stfp, "%4s %8.4f%c %9.4f%c %5d %3d %3d %8.4f %8.4f\n",
                    stn[m], lat, cns, lon, cew, ielev[m], model[m], map1[m], ptcor[m], stcor[m]);
        }

        if (dd[m] > 0.0 && icc[m] >= 2) {
            double var = (bb[m] - aa[m] * aa[m] / dd[m]) * (double)icc[m] / (dd[m] * (double)(icc[m] - 1));
            if (var < 0.0) var = 0.0;
            bb[m] = sqrtf(var);
            aa[m] /= dd[m];
            ee[m] /= (double)icc[m];
            fprintf(stdout, " %4s   %c%4d%8.4f%8.4f%8.4f %8.4f%8.4f\n",
                    stn[m], phz[1], icc[m], ee[m], aa[m], bb[m], dd[m], ptcor[m]);
        }

        if (nsp != 1 && dds[m] > 0.0 && iccs[m] >= 2) {
            double var = (bbs[m] - aas[m] * aas[m] / dds[m]) * (double)iccs[m] / (dds[m] * (double)(iccs[m] - 1));
            if (var < 0.0) var = 0.0;
            bbs[m] = sqrtf(var);
            aas[m] /= dds[m];
            ees[m] /= (double)iccs[m];
            double stcor1 = (nsp == 3) ? (ptcor[m] * vpvs) : stcor[m];
            fprintf(stdout, " %4s   %c%4d%8.4f%8.4f%8.4f %8.4f%8.4f\n",
                    stn[m], phz[0], iccs[m], ees[m], aas[m], bbs[m], dds[m], stcor1);
        }
    }

    fprintf(stdout, "\n\n");
    if (stfp != NULL) {
        fprintf(stfp, "\n");
        fclose(stfp);
    }

    if (istaout == 2) {
        FILE *fp = fopen("velout.mod", "w");
        if (fp != NULL) {
            fprintf(fp, "Output model:\n");
            for (int m = 0; m < nmod; ++m) {
                fprintf(fp, "%d\n", nplay[m]);
                for (int i = 0; i < nplay[m]; ++i) {
                    fprintf(fp, "%5.2f     %7.2f  %7.3f\n", vp[m][i], hp[m][i], vdamp[m][i]);
                }
            }
            fclose(fp);
        }
    }
}
