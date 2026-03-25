#include <stdbool.h>
#include <stdio.h>

#define IEQ 658
#define IST 650

extern int nsta, legs, nsp, nstaeff;
extern int knobs[IEQ], istm[IST][IEQ], nactualsta[IST];
extern double sphase[IST][IEQ];
extern char stn[IST][5];
extern bool single_turbo;
extern FILE *fm_ptr;

void actualstations(void) {
    int nobsp[IST] = {0};
    int nobss[IST] = {0};
    int nofreadings = 0;
    int nofreadp = 0;
    int nofreads = 0;

    for (int i = 0; i < nsta; ++i) {
        nactualsta[i] = 0;
        nobsp[i] = 0;
        nobss[i] = 0;
    }

    nstaeff = 0;

    for (int i = 0; i < legs; ++i) {
        for (int k = 0; k < knobs[i]; ++k) {
            int sta = istm[k][i];
            if (sta < 0 || sta >= nsta) {
                continue;
            }
            nactualsta[sta] += 1;
            if (nsp == 2) {
                if (sphase[k][i] == 0.0f) nobsp[sta] += 1;
                if (sphase[k][i] == 1.0f) nobss[sta] += 1;
            }
            nofreadings += 1;
        }
    }

    if (!single_turbo && fm_ptr != NULL) fprintf(fm_ptr, "\n");

    for (int i = 0; i < nsta; ++i) {
        if (nactualsta[i] > 0) {
            if (!single_turbo && fm_ptr != NULL) {
                if (nsp == 2) {
                    fprintf(fm_ptr, " readings for station %.4s : tot=%4d  P:%4d  S:%4d\n",
                            stn[i], nactualsta[i], nobsp[i], nobss[i]);
                } else {
                    fprintf(fm_ptr, " readings for station %.4s :%4d\n", stn[i], nactualsta[i]);
                }
            }
            nstaeff += 1;
        }
        nofreadp += nobsp[i];
        nofreads += nobss[i];
    }

    if (!single_turbo && fm_ptr != NULL) {
        fprintf(fm_ptr, "\n");
        fprintf(fm_ptr, " Total number of stations with readings:%4d\n", nstaeff);
        fprintf(fm_ptr, "\n");
        fprintf(fm_ptr, " Total number of readings: %7d\n", nofreadings);
        fprintf(fm_ptr, " Total number of P readings: %7d\n", nofreadp);
        fprintf(fm_ptr, " Total number of S readings: %7d\n", nofreads);
        fprintf(fm_ptr, "\n");
    }
}