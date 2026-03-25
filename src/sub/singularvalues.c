/*
 * Compute singular values (eigenvalues) of symmetric matrix G
 * from Fortran SINGULARVALUES subroutine
 */

#include <stddef.h>
#include <stdio.h>
#include <stdbool.h> 

/* External global variables from vel_com.inc */
extern int nvar;
extern int nitt;
extern int ifixsolution;
extern bool single_turbo;
extern double g[];
extern double rht[];
extern double sv[];
extern double ale[];
extern double covs[];
extern double rs[];

/* External functions */
extern void svdsoluk(double A[4][4], double *RHt_ptr, int n, double scale,
                     double Xsol[4], double SV_ptr[], double *ale_val,
                     double COVs_ptr[], double Rs_ptr[]);
extern void alesubr(double SV_ptr[], int mode, double *ale_val);

void singularvalues(int i) {
    double A[4][4];
    double Xsol[4];
    int k, ii, jj, i1, j;
    
    /* Extract columns/rows from symmetric matrix G where k-th row or column appears */
    for (k = 0; k < nvar; k++) {
        ii = 0;
        jj = 0;
        /* Iterate over compressed symmetric storage (lower triangle + diagonal) */
        for (i1 = 0; i1 < nvar; i1++) {
            for (j = 0; j <= i1; j++) {
                /* Include element if it involves row/column k */
                if (i1 == k || j == k) {
                    A[k][ii] = g[jj];  /* 0-based access to g array */
                    ii++;
                }
                jj++;
            }
        }
    }
    
    /* Perform singular value decomposition */
    svdsoluk(A, rht, nvar, -1.0, Xsol, sv, &ale[i], covs, rs);
    
    /* Apply eigenvalue-based ALE if requested */
    if (ifixsolution == 1) {
        alesubr(sv, 3, &ale[i]);
    }
    if (ifixsolution == 9) {
        alesubr(sv, 1, &ale[i]);
    }
    
    /* Debug output if not in turbo mode */
    if (!single_turbo) {
        fprintf(stderr, "Singular values; iteration #%d\n", nitt);
        fprintf(stderr, "  ");
        for (jj = 0; jj < nvar; jj++) {
            fprintf(stderr, "%10.6f  ", sv[jj]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "ALE = %g\n", ale[i]);
    }
}
