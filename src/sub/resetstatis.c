/*
 * Reset statistical variables for ray tracing
 * Called by TRAVELTIME subroutine
 * from Fortran RESETSTATIS subroutine
 */

#include <stddef.h>
#include <stdio.h>
#include "../include/globals.h"

/* External global variables and file descriptors from vel_com.inc */
extern int irflout;
extern int irfrout;
extern int iresout;
extern int nltot;
extern int nsta;
extern FILE *file77;
extern FILE *file78;
extern FILE *file79;
extern int irefllayer[INLTOT];
extern int irefrlayer[INLTOT];
extern double refraylen[INLTOT];
extern double hitlay[INLTOT][3];
extern int noheadwave;
extern double avhraylen;
extern double avvraylen;
extern double sterr;
extern double direrr;
extern double refrerr;
extern double reflerr;
extern double avrefrres;
extern double avotheres;
extern double avreflres;
extern double abrefrres;
extern double abotheres;
extern double abreflres;
extern int nrrefrres;
extern int nrotheres;
extern int nrreflres;
extern double stnazires[IST][8];

void resetstatis(void) {
    int ihitl, jhitl;
    
    /* Rewind output files if active */
    if (irflout == 1 && file77 != NULL) {
        rewind(file77);
    }
    if (irfrout == 1 && file78 != NULL) {
        rewind(file78);
    }
    if (iresout == 1 && file79 != NULL) {
        rewind(file79);
    }
    
    /* Reset layer hit statistics */
    for (ihitl = 0; ihitl < nltot; ihitl++) {
        irefllayer[ihitl] = 0;
        irefrlayer[ihitl] = 0;
        refraylen[ihitl] = 0.0;
        hitlay[ihitl][0] = 0.0;
        hitlay[ihitl][1] = 0.0;
        hitlay[ihitl][2] = 0.0;
    }
    
    /* Reset ray type statistics */
    noheadwave = 0;      /* Number of straight & direct waves */
    avhraylen = 0.0;     /* Average horizontal ray length */
    avvraylen = 0.0;     /* Average vertical ray length */
    sterr = 0.0;         /* Straight/direct wave raytracer error */
    direrr = 0.0;        /* Direct wave raytracer error */
    refrerr = 0.0;       /* Refracted wave raytracer error */
    reflerr = 0.0;       /* Reflected wave raytracer error */
    avrefrres = 0.0;     /* Average residual for refracted rays */
    avotheres = 0.0;     /* Average residual for other rays */
    avreflres = 0.0;     /* Average residual for reflected rays */
    abrefrres = 0.0;     /* Average absolute residual for refracted */
    abotheres = 0.0;     /* Average absolute residual for other */
    abreflres = 0.0;     /* Average absolute residual for reflected */
    nrrefrres = 0;       /* Count of refracted ray residuals */
    nrotheres = 0;       /* Count of other ray residuals */
    nrreflres = 0;       /* Count of reflected ray residuals */
    
    /* Reset station-azimuth residual statistics */
    for (ihitl = 0; ihitl < nsta; ihitl++) {
        for (jhitl = 0; jhitl < 8; jhitl++) {
            stnazires[ihitl][jhitl] = 0.0;
        }
    }
}
