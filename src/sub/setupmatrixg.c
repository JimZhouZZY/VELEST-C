#include "../include/globals.h"

extern int nvar;
extern int isingle;
extern int neqs;
extern int nshot;
extern int nltot;
extern int ksta;
extern int nsta;
extern int nsp;
extern int nshfix;
extern int nshcor;
extern double zadj;
extern double veladj;
extern double scale[7];
extern double w[IST][IEQ];
extern double sphase[IST][IEQ];
extern double dtdr[3];
extern double dtdv[INLTOT];
extern int map2[INSHOT];
extern int map1[IST];
extern int istm[IST][IEQ];
extern int ifx[IEQ];
extern int iphase[IST][IEQ];
extern int nplay[INLTOT];
extern int laysum[INLTOT];
extern double vdamp[ITOTMODELS][INLTOT];
extern double s[INVA * 2];
extern double res[IST][IEQ];
extern double g[];
extern double rht[];

extern void accunormeqs(double svec[], int nvar_local, double residual, double weight, double gmat[], double rhtvec[]);

void setupmatrixg(int neq, int i) {
    int j, nn, nni, mm, is1, kk1, ng, ksta1, ksta2, k, k2;
    int nl, nf, mf, n, k1, is;

    int ev = neq - 1;
    int obs = i - 1;

    for (j = 0; j < nvar; ++j) {
        s[j] = 0.0;
    }

    if (isingle != 0 && w[obs][ev] == 0.0) {
        goto accumulate;
    }

    nn = 4 * neq - 2;
    nni = nn - 1;
    if (neq > neqs) {
        nni = 3 * neqs + neq;
    }
    mm = nn + 2;
    if (neq > neqs) {
        mm = nni;
    }

    s[nni - 1] = 1.0f * scale[0];
    if (sphase[obs][ev] == 2.0f) {
        s[nni - 1] = 0.0;
    }

    is1 = 0;
    if (neq <= neqs) {
        goto hypo_done;
    }

    s[nni - 1] = 1.0f * scale[6];
    kk1 = map2[neq - neqs - 1];
    if (nshfix == 1 && kk1 != 0) {
        s[nni - 1] = 0.0;
    }
    if (nshcor == 0) {
        goto hypo_done;
    }
    if (kk1 <= 0 || kk1 > ksta) {
        goto hypo_done;
    }
    if (scale[4] == 0.0) {
        goto hypo_done;
    }

    if (nsp == 3 && sphase[obs][ev] == 1.0f) {
        goto hypo_done;
    }

    if (nsp == 2 && sphase[obs][ev] == 1.0f) {
        ng = 4 * neqs + nshot + nltot + (ksta / 2) + 1;
        ksta2 = ksta - (ksta / 2);
        if (kk1 > ksta2) {
            goto hypo_done;
        }
        is1 = ng - 1 + kk1;
    } else {
        ng = 4 * neqs + nshot + nltot + 1;
        ksta1 = ksta;
        if (kk1 > ksta1) {
            goto hypo_done;
        }
        is1 = ng - 1 + kk1;
    }

    s[is1 - 1] = 1.0f * scale[4];

hypo_done:
    k = 0;
    if (zadj == 0.0) {
        dtdr[2] = 0.0;
    }

    if (neq <= neqs) {
        if (ifx[ev] == 1) {
            dtdr[1] = 0.0;
        }
        for (j = nn; j <= mm; ++j) {
            k += 1;
            s[j - 1] = dtdr[k - 1] * scale[k];
        }
    }

    if (scale[5] != 0.0) {
        k2 = iphase[obs][ev];
        if (k2 < 1) {
            k2 = 1;
        }
        if (k2 > ITOTMODELS) {
            k2 = ITOTMODELS;
        }

        nl = nplay[k2 - 1];
        nf = 4 * neqs + nshot + laysum[k2 - 1];
        mf = nf + nl - 1;

        k = 0;
        for (n = nf; n <= mf; ++n) {
            k += 1;
            if (veladj == 0.0) {
                dtdv[k - 1] = 0.0;
            }
            if (vdamp[k2 - 1][k - 1] != 0.0) {
                s[n - 1] = dtdv[k - 1] * scale[5] / vdamp[k2 - 1][k - 1];
            }
        }
    }

    if (scale[4] == 0.0) {
        goto accumulate;
    }

    k1 = istm[obs][ev];
    ksta1 = ksta;
    if (nsp == 2) {
        ksta1 = ksta / 2;
    }

    if (nsp == 2 && sphase[obs][ev] == 1.0f) {
        ng = 4 * neqs + nshot + nltot + (ksta / 2) + 1;
        ksta2 = ksta - (ksta / 2);
        if (k1 < 0 || k1 >= nsta) {
            goto accumulate;
        }
        if (map1[k1] == 0 || map1[k1] > ksta2) {
            goto accumulate;
        }
        is = ng - 1 + map1[k1];
    } else {
        if (nsp == 3 && sphase[obs][ev] == 1.0f) {
            goto accumulate;
        }
        ng = 4 * neqs + nshot + nltot + 1;
        if (k1 < 0 || k1 >= nsta) {
            goto accumulate;
        }
        if (map1[k1] == 0 || map1[k1] > ksta1) {
            goto accumulate;
        }
        is = ng - 1 + map1[k1];
    }

    s[is - 1] = 1.0f * scale[4];

accumulate:
    accunormeqs(s, nvar, res[obs][ev], w[obs][ev], g, rht);
}
