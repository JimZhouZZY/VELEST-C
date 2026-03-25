/*
 * Read seismic station input parameters from file
 * Retrieves instrument and filter information by station and time
 * from Fortran STINP subroutine
 */

#include <stddef.h>
#include <stdio.h>
#include <string.h>

#define MAX_STLINES 600

static char stlin[MAX_STLINES][81];
static int maxstlin = 0;
static int ifirstcall = 0;

/* External function */
extern int juliam(int iyr, int imo, int idy, int ihr, int imn);

void stinp(int itime, char *stn, char *sfreq, int *isper, int *iscon,
           int *isdmp, int *isamp, double *scor, int *ier) {
    int ilin, ierr, isyr1, ismo1, isdy1, ishr1, ismin1;
    int isyr2, ismo2, isdy2, ishr2, ismin2, iscor;
    int istime1, istime2;
    char snam6[7];
    char stn6[7];
    char line_buffer[81];
    FILE *file10 = NULL;
    
    /* First call: read station parameter file */
    if (ifirstcall != 10000001) {
        ifirstcall = 10000001;
        maxstlin = 0;
        
        file10 = fopen("STLIST.DAT", "r");
        if (file10 == NULL) {
            *ier = -6;
            return;
        }
        
        while (maxstlin < MAX_STLINES &&
               fgets(line_buffer, sizeof(line_buffer), file10) != NULL) {
            strncpy(stlin[maxstlin], line_buffer, 80);
            stlin[maxstlin][80] = '\0';
            maxstlin++;
        }
        fclose(file10);
    }
    
    /* Search for station in parameter file */
    strncpy(stn6, stn, 4);
    stn6[4] = '\0';
    
    ilin = 0;
    ierr = 0;
    
    while (1) {
        if (ilin >= maxstlin) {
            *ier = -6;  /* Station not found */
            return;
        }
        
        /* Extract station name (first 4-6 chars) */
        strncpy(snam6, stlin[ilin], 6);
        snam6[6] = '\0';
        
        if (strncmp(snam6, stn6, 4) != 0) {
            ilin++;
            continue;
        }
        
        /* Found station; read parameter lines */
        ilin++;
        if (ilin >= maxstlin) {
            *ier = -6;
            return;
        }
        
        /* Parse parameter line: year mo dy hr mn year mo dy hr mn freq ... */
        if (sscanf(stlin[ilin], "%4d%2d%2d%2d%2d %4d%2d%2d%2d%2d %2s%2d%3d%3d%6d%6lf",
                   &isyr1, &ismo1, &isdy1, &ishr1, &ismin1,
                   &isyr2, &ismo2, &isdy2, &ishr2, &ismin2,
                   sfreq, isper, iscon, isdmp, isamp, scor) < 16) {
            ilin++;
            continue;
        }
        
        /* Check time validity */
        if (isyr1 != 0) {
            istime1 = juliam(isyr1, ismo1, isdy1, ishr1, ismin1);
            istime2 = juliam(isyr2, ismo2, isdy2, ishr2, ismin2);
            
            if (itime < istime1 || itime > istime2) {
                ilin++;
                if (ilin >= maxstlin || stlin[ilin][0] == '\0') {
                    *ier = -6;
                    return;
                }
                continue;
            }
        }
        
        /* Skip to end of parameter block */
        ilin++;
        while (ilin < maxstlin) {
            int test_yr = 0;
            if (sscanf(stlin[ilin], "%4d", &test_yr) != 1 || test_yr == 0) {
                break;
            }
            ilin++;
        }
        
        iscor = (int)(*scor * 1000.0f);
        *ier = 0;
        return;
    }
}
