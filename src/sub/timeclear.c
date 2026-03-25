/*
 * Normalize time components (handle out-of-range seconds/minutes)
 * from Fortran TIMECLEAR subroutine and JULIAM function
 */

#include <stddef.h>
#include <stdio.h>

/* External function */
extern void datum(int itf, int *iyr, int *imo, int *iday, int *ihr, int *imin);

/* JULIAM: Convert year-month-day-hour-minute to minutes from epoch */
static int juliam(int iyr, int imo, int idy, int ihr, int imn) {
    int ky, km, kd;
    int ky4, ky1, ky0, kl, l;
    int kmo[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    int leap = 1;
    int result;
    
    ky = iyr;
    km = imo;
    if (km < 1 || km > 12) {
        fprintf(stderr, "Error: invalid month %d\n", km);
        return -1;
    }
    kd = idy;
    if (km <= 0) km = 1;
    
    /* Calculate day number */
    result = 365 * ky;
    kd = kmo[km - 1] + kd;
    
    /* Calculate leap day adjustments */
    ky4 = ky / 4;
    ky1 = ky / 100;
    ky0 = ky / 1000;
    kl = leap * (ky4 - ky1 + ky0);
    
    l = 0;
    if (ky4 * 4 == ky && (ky1 * 100 != ky || ky0 * 1000 == ky)) {
        l = leap;
    }
    if (l != 0 && km < 3) {
        kl = kl - leap;
    }
    
    result = result + kd + kl;
    
    /* Convert to hours then minutes */
    result = result * 24 + ihr;
    result = result * 60 + imn;
    
    return result;
}

void timeclear(int *iyr, int *imo, int *iday, int *ihr, int *imin,
               double *sec, int *itime) {
    int iyr1;
    double sec1;
    
    /* Handle 2-digit year format */
    if (*iyr < 1900) {
        iyr1 = *iyr + 1900;
    } else {
        iyr1 = *iyr;
    }
    
    /* Convert to Julian minutes */
    *itime = juliam(iyr1, *imo, *iday, *ihr, *imin);
    
    /* Normalize seconds (handle < 0 or > 60) */
    sec1 = *sec;
    while (sec1 < 0.0f || sec1 > 60.0f) {
        if (sec1 < 0.0f) {
            *sec = sec1 + 60.0f;
            *itime = *itime - 1;
        } else if (sec1 > 60.0f) {
            *sec = sec1 - 60.0f;
            *itime = *itime + 1;
        }
        sec1 = *sec;
    }
    
    /* Update date/time components */
    datum(*itime, &iyr1, imo, iday, ihr, imin);
    
    /* Restore original year format if needed */
    if (*iyr < 1900) {
        *iyr = iyr1 - 1900;
    } else {
        *iyr = iyr1;
    }
}
