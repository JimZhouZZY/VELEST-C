#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

/*
 * Modern lsvdf: compute SVD of A, fill singular values s[], U->a, V^T->b
 * wk is set to zero (interface compatible)
 */
void lsvdf(double a[][4], int ia, int m, int n,
           double b[][4], int ib, int nb,
           double s[], double wk[][2], int *ier) {
    (void)ia;
    (void)ib;

    /*
    // Debug print all input variables
    printf("lsvdf called with:\n");
    printf("  ia = %d, m = %d, n = %d, ib = %d, nb = %d\n", ia, m, n, ib, nb);
    printf("  a (input matrix):\n");
    for (int i = 0; i < m; ++i) {
        printf("    ");
        for (int j = 0; j < n; ++j) {
            printf("%g ", a[i][j]);
        }
        printf("\n");
    }
    printf("  b (input matrix):\n");
    for (int i = 0; i < m; ++i) {
        printf("%g ", b[0][i]);
    }
    printf("\n");
    printf("  s (input, before SVD):\n    ");
    for (int i = 0; i < ((m < n) ? m : n); ++i) {
        printf("%g ", s[i]);
    }
    printf("\n");
    printf("  wk (input, before SVD):\n");
    for (int i = 0; i < n; ++i) {
        printf("    %g %g\n", wk[i][0], wk[i][1]);
    }
    printf("\n");
    */

    *ier = 0;
    if (m <= 0 || n <= 0) {
        *ier = 34;
        return;
    }

    const double zero = 0.0;

    // Flatten A for LAPACK row-major
    double *Aflat = (double*)malloc(sizeof(double) * m * n);
    if (!Aflat) { *ier = 100; return; }
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            Aflat[i*n + j] = a[i][j];

    // Allocate U (m x m) and VT (n x n)
    double *U = (double*)calloc(m * m, sizeof(double));
    double *VT = (double*)calloc(n * n, sizeof(double));
    if (!U || !VT) { *ier = 101; free(Aflat); if(U) free(U); if(VT) free(VT); return; }

    // Workspace array for LAPACK
    double *superb = (double*)malloc(sizeof(double) * (n-1));
    if (!superb) { *ier = 102; free(Aflat); free(U); free(VT); return; }

    // Compute full SVD: A = U * diag(s) * VT
    int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A',
                              m, n, Aflat, n, s, U, m, VT, n, superb);
    free(superb);
    if (info != 0) {
        *ier = info;
        free(Aflat); free(U); free(VT);
        return;
    }

    // wk[][] set to zero to maintain interface
    for (int j = 0; j < n; ++j) {
        wk[j][0] = zero;
        wk[j][1] = zero;
    }

    // Copy U into a[m][4] (row-major)
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            a[i][j] = (j < m) ? U[i*m + j] : zero;

    // Copy VT into b[m][4] (row-major)
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < nb; ++j)
            b[0][i] = (j < n) ? VT[j*n + i] : zero;

    *ier = 33; // success
    /*
    printf("  s (input, after SVD):\n    ");
    for (int i = 0; i < ((m < n) ? m : n); ++i) {
        printf("%g ", s[i]);
    }
    printf("\n");
    printf("  wk (input, after SVD):\n");
    for (int i = 0; i < n; ++i) {
        printf("    %g %g\n", wk[i][0], wk[i][1]);
    }
    printf("\n");
    */

    free(Aflat); free(U); free(VT);
}