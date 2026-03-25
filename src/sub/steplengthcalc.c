#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#define IEQ 658
#define INSHOT 50
#define INLTOT 100
#define IST 650
#define ILIP (4 * IEQ + INSHOT + INLTOT + IST)
#define INVA (ILIP - 1)

extern int nvar;
extern int ibackups;
extern bool single_turbo;
extern double steplen;
extern double b[INVA * 2];
extern FILE *fm_ptr;

void steplengthcalc(void) {
    steplen = 0.0f;
    for (int i = 0; i < nvar; ++i) {
        if (!isfinite(b[i])) {
            FILE *out = fm_ptr ? fm_ptr : stdout;
            fprintf(out, "DEBUG STEPLEN: non-finite b[%d]=%g\n", i, (double)b[i]);
            steplen = NAN;
            return;
        }
        steplen += b[i] * b[i];
    }
    steplen = sqrtf(steplen);

    if (ibackups > 0) {
        steplen = -steplen;
    }

    if (!single_turbo) {
        FILE *out = fm_ptr ? fm_ptr : stdout;
        fprintf(out, "\n");
        fprintf(out, " (Applied) Step length = %7.3f\n", steplen);
    }
}
