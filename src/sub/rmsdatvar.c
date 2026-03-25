#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#define IEQ 658
#define IST 650

extern int legs;
extern int nitt;
extern int isingle;
extern int nobswithw0;
extern int nvareff;
extern bool single_turbo;
extern int knobs[IEQ];
extern int iyr[IEQ];
extern int imo[IEQ];
extern int iday[IEQ];
extern int ihr[IEQ];
extern int imin[IEQ];
extern double w[IST][IEQ];
extern double res[IST][IEQ];
extern double avres[IEQ];
extern double rms[IEQ];
extern double davar1;
extern double xmsqrs1;
extern FILE *fm_ptr;

void rmsdatvar(void) {
    int knobst = 0;
    double tres = 0.0f;
    FILE *logfp = (fm_ptr != NULL) ? fm_ptr : stdout;

    davar1 = 0.0f;
    xmsqrs1 = 0.0f;

    for (int i = 0; i < legs; ++i) {
        int j2 = knobs[i];
        knobst += j2;

        for (int j = 0; j < j2; ++j) {
            if (w[j][i] <= 0.0f) {
                continue;
            }
            avres[i] += res[j][i] * w[j][i];
            rms[i] += (res[j][i] * res[j][i]) * w[j][i] * w[j][i];
        }

        davar1 += rms[i];
        tres += avres[i];

        if ((j2 - nobswithw0) <= 1) {
            continue;
        }

        rms[i] = sqrtf((rms[i] - (avres[i] * avres[i]) / (j2 - nobswithw0))
                     / ((j2 - nobswithw0) - 1));
        avres[i] = avres[i] / (j2 - nobswithw0);
    }

    if (nitt == 0 && !single_turbo) {
        fprintf(logfp, "\n Total number of observations is: %5d\n\n", knobst);
    }

    if (nitt == 0 && isingle != 0 && !single_turbo) {
        fprintf(logfp, "Number of observations with normalized weight 0.0 : %d\n", nobswithw0);
        fprintf(logfp, "knobs(i)   = %d\n", knobst);
        fprintf(logfp, "nobswithw0 = %d\n", nobswithw0);
    }

    if ((knobst - nobswithw0) > nvareff) {
        davar1 = (davar1 - tres * tres / (knobst - nobswithw0))
               / ((knobst - nobswithw0) - nvareff);
        xmsqrs1 = davar1 * ((knobst - nobswithw0) - nvareff) / (knobst - nobswithw0);
    } else {
        davar1 = 999.99f;
        xmsqrs1 = 999.99f;
    }

    if (!single_turbo && isingle == 0) {
        fprintf(logfp, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
        fprintf(logfp, "Events with  | AVRES | > 1.0 SEC are suspicious !\n\n");

        int found = 0;
        for (int i = 0; i < legs; ++i) {
            if (fabs(avres[i]) > 1.0f) {
                fprintf(logfp,
                        " Event# %3d >>> %02d%02d%02d %02d%02d  AVRES =%6.2f NOBS =%3d\n",
                        i + 1, iyr[i], imo[i], iday[i], ihr[i], imin[i], avres[i], knobs[i]);
                found++;
            }
        }

        if (found > 0) {
            fprintf(logfp, "\n^^^^^^^^^ C H E C K   these events above ^^^^^^^^\n");
            fprintf(logfp, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
        } else {
            fprintf(logfp, "ZERO events of this kind found! (lucky guy!!!)\n");
        }
        fprintf(logfp, "\n");
    }
}
