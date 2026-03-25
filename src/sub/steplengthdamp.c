#include <math.h>

#define IEQ 658
#define INSHOT 50
#define INLTOT 100
#define IST 650
#define ILIP (4 * IEQ + INSHOT + INLTOT + IST)
#define INVA (ILIP - 1)

extern int neqs;
extern int nshot;
extern int nmod;
extern int nvar;
extern int icount;
extern int nplay[INLTOT];
extern int laysum[INLTOT];
extern double scale[7];
extern double veladj;
extern double b[INVA * 2];

void steplengthdamp(double *damp) {
    double btemp;
    int j11, j22;

    if (!damp) {
        return;
    }

    *damp = 1.0f;
    if (icount == 1) {
        return;
    }
    if (scale[5] == 0.0f) {
        return;
    }

    for (int i = 0; i < nmod; ++i) {
        j11 = 4 * neqs + nshot + laysum[i];
        j22 = j11 + nplay[i] - 1;
        for (int jjj = j11; jjj <= j22; ++jjj) {
            int idx = jjj - 1;
            if (idx < 0 || idx >= nvar) {
                continue;
            }
            if (veladj > fabs(b[idx])) {
                continue;
            }
            btemp = veladj / fabs(b[idx]);
            if (btemp < *damp) {
                *damp = btemp;
            }
        }
    }

    for (int jjj = 0; jjj < nvar; ++jjj) {
        b[jjj] *= *damp;
    }
}
