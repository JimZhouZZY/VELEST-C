/*
 * PDR2 seismic amplitude response calculation
 * Computes transfer function for PDR2 recording system
 * from Fortran MPDR2 subroutine
 */

#include <stddef.h>
#include <math.h>
#include <complex.h>

void mpdr2(int isecpendel, double seismkonst, double seismdamp,
           double voltgain, double pdrampl, double period, int iampltyp,
           double *ampl) {
    double complex j = CMPLX(0.0, 1.0f);  /* Imaginary unit */
    double complex G, gseis, ghp, glp2, glp2a, glp3, jom;
    double zpi = 6.28319f;
    double ts, hs, om, om2, omg1, omg2, omg3, omg4, omg5, twa;
    double omw, hw, vwa;
    double complex wa_response;
    
    (void)seismkonst;  /* Parameter defined but not used in this version */
    (void)voltgain;    /* Parameter not used in amplitude calculation */
    (void)iampltyp;    /* Parameter not used in this calculation */
    ts = (double)isecpendel;
    hs = seismdamp;
    om = zpi / period;
    jom = j * om;
    om2 = om * om;
    
    /* Seismometer response (displacement transducer, mechanics only) */
    omg1 = zpi / ts;
    gseis = (jom * om2) / (omg1 * omg1 + 2.0f * jom * omg1 * hs - om2);
    
    /* Two 1st-order high-pass filters => 2nd-order, fc=0.3 Hz */
    omg2 = zpi * 0.3f;
    ghp = 1.0f / (1.0f - 2.0f * j * omg2 / om - omg2 * omg2 / om2);
    
    /* Butterworth low-pass 2nd order, fc=30 Hz */
    omg3 = zpi * 30.0;
    glp2 = 1.0f / (1.0f + 1.4142f * j * om / omg3 - om2 / (omg3 * omg3));
    
    /* Butterworth low-pass 3rd order, fc=25 Hz (Discriminator) */
    omg4 = zpi * 25.0f;
    glp3 = 1.0f / (1.0f + 2.0f * j / omg4 - 2.0f * om2 / (omg4 * omg4) 
                   - jom * om2 / (omg4 * omg4 * omg4));
    
    /* Butterworth low-pass 2nd order, fc=24 Hz (Anti-aliasing) */
    omg5 = zpi * 24.0f;
    glp2a = 1.0f / (1.0f + 1.4142f * jom / omg5 - om2 / (omg5 * omg5));
    
    /* Total system transfer function */
    G = gseis * ghp * glp2 * glp3 * glp2a;
    
    /* Deconvolution with Wood-Anderson seismometer */
    twa = 0.8f;
    omw = zpi / twa;
    hw = 0.78f;
    vwa = 2800.0;
    wa_response = (om2 * vwa) / (omw * omw + 2.0f * jom * omw * hw - om2);
    G = wa_response / G;
    
    /* Compute amplitude as magnitude of complex transfer function */
    *ampl = cabs(G);  /* Amplitude response [sec] */
    *ampl = (*ampl) * pdrampl / 10000.0;  /* Convert to mm */
}
