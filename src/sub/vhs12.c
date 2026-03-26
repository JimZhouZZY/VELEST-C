/*
 * Householder matrix transformation (IMSL VHS12)
 * Applies orthogonal Householder transformation to a matrix
 * from Fortran VHS12 subroutine
 */

#include <stddef.h>
#include <stdio.h>
#include <math.h>

void vhs12(int mode, int lp, int l1, int m, double u[], double *up,
           double c[], int incu, int incc, int icv, int ncv) {
    int ij, ilp, il1, im, incr, i2, i3, i4, j;
    double sm, b;
    double one = 1.0f;
    double cl, clinv, sm1;
    
    /* Validate input parameters */
    if (lp <= 0 || lp >= l1 || l1 > m) {
        fprintf(stderr, "vhs12: invalid parameters or insufficient array sizes\n");
        return;
    }
    
    /* Calculate indices into arrays */
    ilp = (lp - 1) * incu;      /* 0-based index for LP element */
    il1 = (l1 - 1) * incu;      /* 0-based index for L1 element */
    im = (m - 1) * incu;        /* 0-based index for M element */
    
    cl = fabs(u[ilp]);
    
    if (mode != 2) {
        /* CONSTRUCT THE TRANSFORMATION */
        for (ij = il1; ij <= im; ij += incu) {
            cl = fmaxf(fabs(u[ij]), cl);
        }
        
        if (cl <= 0.0) {
            return;
        }
        
        clinv = one / cl;
        sm = (u[ilp] * clinv) * (u[ilp] * clinv);
        
        for (ij = il1; ij <= im; ij += incu) {
            sm = sm + (u[ij] * clinv) * (u[ij] * clinv);
        }
        
        sm1 = sm;  /* Convert to single precision */
        cl = cl * sqrtf(sm1);
        
        if (u[ilp] > 0.0) {
            cl = -cl;
        }
        
        *up = u[ilp] - cl;
        u[ilp] = cl;
    }
    
    /* APPLY THE TRANSFORMATION */
    if (ncv <= 0) {
        return;
    }
    
    if (cl <= 0.0) {
        return;
    }
    
    b = (*up) * u[ilp];
    
    /* B must be non-positive; if B == 0, return */
    if (b >= 0.0) {
        return;
    }
    
    b = one / b;
    
    /* Apply transformation I + U*(U^T)/B to C */
    i2 = 1 - icv + incc * (lp - 1);
    incr = incc * (l1 - lp);
    
    for (j = 0; j < ncv; j++) {
        i2 = i2 + icv;
        i3 = i2 + incr;
        i4 = i3;
        
        sm = c[i2] * (*up);
        
        // SEGFAULT here
        fprintf(stderr, "DEBUG: il1=%d, im=%d, incu=%d, len_c=%ld, len_u=%ld\n", il1, im, incu, (long)(sizeof(c)/sizeof(c[0])), (long)(sizeof(u)/sizeof(u[0])));
        for (ij = il1; ij <= im; ij += incu) {
            sm = sm + c[i3] * u[ij];
            i3 = i3 + incc;
        }
        
        if (sm == 0.0) {
            continue;
        }
        
        sm = sm * b;
        c[i2] = c[i2] + sm * (*up);
        
        for (ij = il1; ij <= im; ij += incu) {
            c[i4] = c[i4] + sm * u[ij];
            i4 = i4 + incc;
        }
    }
}
