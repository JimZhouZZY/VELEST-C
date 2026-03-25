/*
 * Ellipsoid to plane coordinate transformation
 * Swiss projection (Bern coordinate system)
 * from Fortran EBELL subroutine
 */

#include <stddef.h>
#include <math.h>

void ebell(double yl, double xb, double *l, double *b, double *my) {
    double a[8], bb[8], c[8], rz[8], qz[8];
    double x, y, p, q, db, dl;
    int i;
    
    /* Convert from meters (Bern 600/200) to internal coordinates */
    x = 1000.0 * (double)xb - 200000.0;
    y = 1000.0 * (double)yl - 600000.0;
    
    /* Coefficients for transformation polynomial */
    a[0] = 1.4623614572021;
    a[1] = 1.225255821052e-7;
    a[2] = 1.3687923002e-14;
    a[3] = 1.971224191e-21;
    a[4] = 2.97898051e-28;
    a[5] = 4.650273e-35;
    a[6] = 0.0;
    a[7] = 0.0;
    
    bb[0] = 3.4564252673326e-2;
    bb[1] = 2.89600437564e-9;
    bb[2] = 4.651046030e-16;
    bb[3] = 6.43850954e-23;
    bb[4] = 9.600412e-30;
    bb[5] = 1.50512e-36;
    bb[6] = 0.0;
    bb[7] = 0.0;
    
    c[0] = 2.2146704979846e-2;
    c[1] = -1.280815253730e-9;
    c[2] = 7.4775676024e-18;
    c[3] = 4.691943327e-24;
    c[4] = -3.6550101e-31;
    c[5] = 3.71615e-39;
    c[6] = 0.0;
    c[7] = 0.0;
    
    /* Compute powers of complex number X+iY */
    rz[0] = x;
    qz[0] = y;
    for (i = 1; i < 6; i++) {
        rz[i] = x * rz[i - 1] - y * qz[i - 1];
        qz[i] = y * rz[i - 1] + x * qz[i - 1];
    }
    rz[6] = x * rz[5] - y * qz[5];
    qz[6] = y * rz[5] + x * qz[5];
    rz[7] = x * rz[6] - y * qz[6];
    qz[7] = y * rz[6] + x * qz[6];
    
    /* Horner's method polynomial evaluation */
    i = 7;
    q = a[i] * rz[i];
    p = a[i] * qz[i];
    *my = 0.0f;
    
    while (i > 0) {
        i--;
        q = q + a[i] * rz[i];
        p = p + a[i] * qz[i];
        *my = *my + (double)(bb[i] * qz[i]);
    }
    
    /* Compute longitude and latitude */
    dl = 3.2343101932327e-2 * p;
    db = q * c[7];
    for (i = 6; i > 0; i--) {
        db = q * (db + c[i]);
    }
    
    *l = (double)((dl + 26782.5) / 3600.0);
    *b = (double)((db + 169028.66) / 3600.0);
}
