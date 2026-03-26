#include <math.h>
#include <stdbool.h>
#include <stdio.h>

extern int nitt;
extern int ibackups;
extern int isingle;
extern bool single_turbo;
extern double davar1;
extern double xmsqrs1;
extern FILE *fm_ptr;

extern void gapcalc(int i);
extern void avresistatist(void);

void checksolution(int *istopflag, int *better) {
    static double datvar = 0.0;
    static double xmsqrs2 = 0.0;

    FILE *logfp = (fm_ptr != NULL) ? fm_ptr : stdout;

    int decreasing = 0;
    *better = 1;

    if (nitt == 0) {
        datvar = davar1;
        xmsqrs2 = xmsqrs1;
    }

    if (!single_turbo) {
        fprintf(logfp, "\n");
        if (ibackups == 0) {
            fprintf(logfp, " Iteration nr %2d obtained:\n", nitt);
        } else {
            fprintf(logfp, "(Iteration nr %2d)   BACKUP nr %1d obtained:\n", nitt, ibackups);
        }
        fprintf(logfp, " DATVAR=%12.6f mean sqrd residual= %12.6f  RMS RESIDUAL= %12.6f\n",
                davar1, xmsqrs1, sqrtf(xmsqrs1));
    }

    if (isingle != 0) {
        gapcalc(0);
    }

    if (isingle == 0) {
        avresistatist();
    }

    double varat1 = datvar / davar1;

    if (isingle != 0) {
        if (nitt > 2 && fabs(datvar - davar1) < 1e-6f) {
            *istopflag = 1;
            if (!single_turbo) {
                fprintf(logfp, "Changes in datvar < 1e-6   : STOPPING...\n");
            }
        }
    }

    double varat2 = xmsqrs2 / xmsqrs1;

    if (varat1 >= 0.99f) {
        decreasing = 1;
    } else if (varat1 < 0.99f && varat2 >= 0.99f) {
        if (!single_turbo) {
            fprintf(logfp, " *** WARNING: the data variance has increased.\n");
        }
        decreasing = 1;
    } else if (varat2 >= 0.99f) {
        decreasing = 1;
    }

    if (decreasing == 1) {
        datvar = davar1;
        xmsqrs2 = xmsqrs1;
    } else {
        *better = 0;
    }
}
