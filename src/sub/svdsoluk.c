#include <math.h>
#include <stdlib.h>
#include <stdio.h>

extern void lsvdf(double a[][4], int lda, int m, int n, double b[], int ldb,
                  int job, double s[], double wk[][2], int *ier);
extern void alesubr(double s[], int m, double *ale);
extern void matrtran(double a[], int n, int m, double at[]);
extern void matrmult(double a[], int m, int p, double b[], int p2, int n, double c[], int m2, int n2);

void svdsoluk(double ain[][4], double bin[], int m, double eigmin,
              double x[], double s[], double *ale, double cov[], double r[][4]) {
    const int mm = 4;
    const int mm2 = 8;

    if (m > 4) {
        fprintf(stderr, "Error: m (%d) > 4 in svdsoluk\n", m);
        abort();
    }

    double a[4][4] = {{0.0}};
    double b[4] = {0.0};
    double at[4][4] = {{0.0}};
    double wk[8][2] = {{0.0}};

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            a[i][j] = ain[i][j];
        }
        b[i] = bin[i];
        s[i] = 0.0;
        cov[i] = 0.0;
    }

    int ier = 0;
    lsvdf(a, m, m, m, b, m, 1, s, wk, &ier);

    alesubr(s, m, ale);

    if (eigmin < 0.0) {
        return;
    }

    int nfre = 0;
    for (int i = 0; i < m; ++i) {
        x[i] = 0.0;
        if (s[i] > eigmin) {
            nfre++;
        } else {
            s[i] = 0.0;
            for (int j = 0; j < m; ++j) {
                a[j][i] = 0.0;
            }
        }
    }

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < nfre; ++j) {
            x[i] += a[i][j] * b[j] / s[j];
        }
    }

    for (int i = 0; i < m; ++i) {
        double sum = 0.0;
        for (int j = 0; j < nfre; ++j) {
            sum += a[i][j] * a[i][j] / (s[j] * s[j]);
        }
        cov[i] = sum;
    }

    matrtran((double *)a, m, m, (double *)at);
    matrmult((double *)a, m, m, (double *)at, m, m, (double *)r, m, m);

    (void)mm;
    (void)mm2;
}
