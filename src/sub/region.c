/*
 * Query geographic region name for given coordinates
 * Supports both Swiss and lat/lon coordinate systems
 * from Fortran REGION subroutine
 */

#include <stddef.h>
#include <string.h>

/* External functions for region lookup */
extern void regread(const char *regnamfile, const char *regkoordfile);
extern void regch(double x, double y, char *place, int *nreg);
extern void geoko(double *x, double *y, double xlat, double xlon, int dir);
extern void regworld(int iii, double xlat, double xlon, char *place, int *nreg);

void region(int ityp, double x, double y, char *cname, int *nreg,
            const char *regnamfile, const char *regkoordfile) {
    static int iregread = 0;
    char place[33];
    double xlat, xlon;
    int iii;
    
    /* Load region data on first call */
    if (iregread == 0) {
        regread(regnamfile, regkoordfile);
        iregread = 1;
    }
    
    *nreg = 0;
    memset(place, 0, sizeof(place));
    
    /* Branch on coordinate type */
    if (ityp == 1) {
        /* Swiss coordinates (x, y) */
        regch(x, y, place, nreg);
        if (*nreg == 0) {
            /* Convert to lat/lon and retry */
            geoko(&xlat, &xlon, x, y, 1);
            goto regworld_lookup;
        } else {
            strncpy(cname, place, 32);
            cname[32] = '\0';
        }
    } else if (ityp == 2 || ityp == 3) {
        /* Lat/lon coordinates */
        xlat = x;
        xlon = y;
        
    regworld_lookup:
        memset(place, 0, sizeof(place));
        iii = (ityp == 3) ? 1 : 0;  /* Extended region lookup if ityp==3 */
        regworld(iii, xlat, xlon, place, nreg);
        
        if (strlen(place) == 0) {
            strncpy(cname, "*****", 32);
        } else {
            strncpy(cname, place, 32);
        }
        cname[32] = '\0';
    } else {
        /* Invalid coordinate type */
        strncpy(cname, "ERROR", 32);
        cname[32] = '\0';
    }
}
