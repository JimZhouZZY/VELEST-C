#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h> 

#define IEQ 658
#define INSHOT 50
#define INLTOT 100
#define IST 650
#define ILIP (4 * IEQ + INSHOT + INLTOT + IST)
#define INVA (ILIP - 1)

extern bool single_turbo;
extern int isingle;
extern int neqs;
extern int nshot;
extern int nltot;
extern int iresolcalc;
extern int ifixsolution;
extern int nvar;
extern int ksta;

extern double g[];
extern double rht[];
extern double gcopy[60][INVA * 2];
extern double scale[];
extern double vdamp[][100];
extern int nplay[];
extern double rc[][4];
extern double covc[][4];
extern double r3[][3];
extern double spread;
extern double sv[];
extern double ale[];
extern double s[];

extern void luelmp(double a[], double b[], int n, double x[]);
extern void fixunt(double b[], int neqs_local, int nshot_local, int nltot_local, int ksta_local,
                   double scale_local[], double vdamp_local[][100], int itotmodels_local,
                   int inltot_local, int nplay1);
extern void spreadd(double r[][4], int n, double *spread_out);
extern void spreadb(double r[][4], int n, double *spread_out);

void resolcovar(double davari) {
    int mef = 4 * neqs + nshot - 4;
    int nef = 4 + nltot;
    if (isingle != 0) {
        nef = 4;
    }

    int n1 = 4 * neqs;
    int n2 = n1 + nshot;
    int n3 = n2 + nltot;

    double rdiag[100] = {0.0f};

    for (int row_index = 0; row_index < nef; ++row_index) {
        for (int col_index = 0; col_index < nvar; ++col_index) {
            rht[col_index] = gcopy[row_index][col_index];
        }

        int line_index = row_index;
        if (row_index >= 4) {
            line_index = row_index + mef;
        }

        luelmp(g, rht, nvar, rht);

        if (isingle == 0) {
            rdiag[row_index] = rht[line_index];
        } else {
            for (int column = 0; column < 4; ++column) {
                rc[line_index][column] = rht[column];
            }
        }

        luelmp(g, rht, nvar, rht);

        int block = (line_index + 1) % 4;
        if (block == 0) {
            block = 4;
        }

        double scale1 = scale[block - 1];
        if (line_index >= n1) {
            scale1 = scale[0];
        }
        if (line_index >= n2) {
            scale1 = scale[5];
        }
        if (line_index >= n3) {
            scale1 = scale[4];
        }

        for (int column = 0; column < nvar; ++column) {
            rht[column] = rht[column] * scale1;
        }

        fixunt(rht, neqs, nshot, nltot, ksta, scale, vdamp, 2, 100, nplay[0]);

        if (isingle != 0) {
            for (int column = 0; column < 4; ++column) {
                covc[line_index][column] = rht[column];
            }
        }

        s[row_index] = sqrtf(fabs(rht[line_index]));
    }

    if (isingle != 0) {
        double spread1 = 0.0f;
        double spread2 = 0.0f;

        spreadd(rc, 4, &spread1);
        spread = spread1;

        if (ifixsolution == 1) {
            for (int row = 0; row < 3; ++row) {
                for (int column = 0; column < 3; ++column) {
                    r3[row][column] = rc[row][column];
                }
            }
            double spread3 = 0.0f;
            spreadd((double (*)[4])r3, 3, &spread3);
            spread = spread3;
        }

        if (ifixsolution == 9) {
            double spread4 = 0.0f;
            spreadd((double (*)[4])rc, 1, &spread4);
            spread = spread4;
        }

        spreadb(rc, 4, &spread2);

        double size = 0.0f;
        for (int index = 0; index < 4; ++index) {
            size += covc[index][index];
        }

        if (!single_turbo) {
            fprintf(stdout, " D-Spread(R) = %6.3f   B-G Spread(R) = %6.3f       Size (C) = %6.3f\n",
                    spread, spread2, size);
        }

        double avresol = 0.0f;
        for (int index = 0; index < nef; ++index) {
            avresol += rc[index][index];
        }
        avresol /= (double)nef;
        (void)avresol;

        fprintf(stdout, "                        OT (sec)   X (km)    Y (km)    Z (km) \n");
        fprintf(stdout, " Sigma (CHD):         %10.4f%10.4f%10.4f%10.4f\n", s[0], s[1], s[2], s[3]);
        fprintf(stdout, " Resolution (CHD):    %10.4f%10.4f%10.4f%10.4f D-spread =%6.3f\n",
                rc[0][0], rc[1][1], rc[2][2], rc[3][3], spread);
        fprintf(stdout, " Data Variance      = %10.4f\n", davari);
        fprintf(stdout, " Singular values:     %10.4f%10.4f%10.4f%10.4f     ALE =%7.3f\n",
                sv[0], sv[1], sv[2], sv[3], ale[0]);
    }

    if (!single_turbo) {
        fprintf(stdout, "\n");
        fprintf(stdout, "Rdiag of selected model parameters:\n");
        for (int index = 0; index < nef; ++index) {
            fprintf(stdout, "%10.4f", rdiag[index]);
        }
        fprintf(stdout, "\n\n");
    }
}
