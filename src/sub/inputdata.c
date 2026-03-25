/*
 * inputdata.c
 * Migrated from inputdata.f (original Fortran, 532 lines)
 * 
 * Reads earthquake and shot phase data, performs validation and filtering.
 * Key Fortran semantics converted to C:
 *  - SAVE variable (iphaseteststopflag) -> static variable
 *  - labeled continue + goto -> while loop with control flow
 *  - 1-based array indices -> 0-based (Fortran e(1,i) -> C e[0][i-1])
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

/* ===== External declarations ===== */
/* These correspond to Fortran COMMON blocks defined in vel_globals.c */

extern int ieq, ist, maxobsperevent, inshot;
extern int nvar, nvareff, kvar, nsta, ksta, legs, lip, neqs, nshot,
           nsinv, nshcor, nshfix, icount;
extern int ittmax, invertratio, iresolcalc, ifixsolution;
extern double zmin, zmininput, delmin, veladj, swtfac, xythet,
            zadj, vthet, othet, zthet, stathet, rmsmin, dmax;
extern double zshift;
extern int icoordsystem, itrial, ised, isingle, itopo, iturbo,
           lowveloclay, iusestacor;
extern int ielev[650];
extern double ztrial;
extern bool single_turbo;

extern int nmod, nsp, nltot;
extern int nplay[100];
extern int ireflector, lmax;
extern double vp[2][100], hp[2][100], vdamp[2][100], thkp[2][100];

extern double ptcor[650], stcor[650];
extern double d[650][3][658];     /* station-event distance [station][xyz][event] */
extern double x[650][3];          /* station coordinates [station][xyz] */
extern int model[2*650];         /* velocity model assignment */
extern int map1[650], map2[650];
extern double xla[650], xlo[650];

extern double pt[650][658];       /* phase time [station][event] */
extern int kpwt[650][658];       /* phase weight [station][event] */
extern int istm[650][658];       /* station number [station][event] */
extern double w[650][658];        /* weight factor [station][event] */
extern int iphase[650][658];     /* model number [station][event] */
extern double sphase[650][658];   /* phase type: 0=P, 1=S, -1=reflected P, 2=S-P */
extern int knobs[658];           /* number of observations per event */
extern double depthsofinput[658];
extern double amx[650], prx[650];
extern int nobswithw0;

extern char smn[650][658][5];    /* station name [station][event] */
extern char stn[650][5];         /* station list [station] */

extern int iyr[658], imo[658], iday[658], ihr[658], imin[658];
extern double e[][658];          /* event: [0]=origin time, [1]=lat, [2]=lon, [3]=depth */
extern double emag[658];
extern int ifx[658];

extern char prmk[650][2];        /* phase remarks for single-event mode */

/* File pointers for I/O (8=earthquakes, 9=shots) */
extern FILE *iunit_fp[10];
extern char phasefile[81], shotfile[81];

extern FILE *fm_ptr;
extern int iabort;

/* Forward declarations of input functions (they expect FILE* pointer) */
void inputcnv(FILE *fp, int *nobs,
              char sta[][5], int *iyr, int *imo, int *iday, int *ihr, int *imin, double sec[],
              char rmk1[], char rmk2[], char cphase[], int iwt[], double amx[], double prx[],
              double *xlat, double *xlon, double *xmagni, double *depth, double *origtime,
              int *i1, int *i2, int *i3, char *eventtype);
void inputarcvel(FILE *fp, int *nobs,
                 char sta[][5], int *iyr, int *imo, int *iday, int *ihr, int *imin, double sec[],
                 char rmk1[], char rmk2[], char cphase[], int iwt[], double amx[], double prx[],
                 double *xlat, double *xlon, double *xmagni, double *depth, double *origtime,
                 int *i1, int *i2, int *i3, char *eventtype);
void inputsed_new(FILE *fp, int *nobs,
                  char sta[][5], int *iyr, int *imo, int *iday, int *ihr, int *imin, double sec[],
                  char rmk1[], char rmk2[], char cphase[], int iwt[], double amx[], double prx[],
                  double *xlat, double *xlon, double *xmagni, double *depth, double *origtime,
                  int *i1, int *i2, int *i3, char *eventtype, int itrial);
void geoko(double *x, double *y, double lat, double lon, int dir);
void sdc(double *x, double *y, double lat, double lon, int dir);

/* Helper: read observation arrays into local buffers */
void inputdata(int event_idx)
{
    static int iphaseteststopflag = 0;  /* SAVE variable: persist across calls */
    
    int nobs, ntobs, jj, nobsread, i1, i2, i3, icc, jjmin1, jk, jjmin, jshot, j, k, iunit, ll, itest, ie;
    int l, nobsp, nobss;
    double wsum, xlat, xlon, alon, depth, z, ss1;
    double xxlat, aalon, xxx;
    
    /* Local arrays: observations per this event */
    char cphase[181];       /* phase type characters */
    char rmk1[181];         /* remark 1 */
    char rmk2[181];         /* remark 2 */
    char sta_local[181][5]; /* station names */
    int ipwt_local[181];    /* phase weights */
    double sec_local[181];   /* seconds */
    
    char sc = 's', ss = 'S';             /* phase type constants */
    char eventtype;
    char cns, cew;                       /* coordinate direction indicators */
    
    int i = event_idx;  /* Convert to 1-based indexing for Fortran logic */
    FILE *fp;           /* File pointer for reading */
    
    /* File opening on first call */
    if (i == 1) {
        if (neqs > 0) {
            if (iunit_fp[8] == NULL) {
                iunit_fp[8] = fopen(phasefile, "r");
                if (!iunit_fp[8]) {
                    fprintf(stderr, "ERROR: Cannot open phase file: %s\n", phasefile);
                    iabort = 1;
                    return;
                }
            }
        }
        if (nshot > 0) {
            if (iunit_fp[9] == NULL) {
                iunit_fp[9] = fopen(shotfile, "r");
                if (!iunit_fp[9]) {
                    fprintf(stderr, "ERROR: Cannot open shot file: %s\n", shotfile);
                    iabort = 1;
                    return;
                }
            }
        }
    }
    
    /* ===== RE-ENTRY LOOP for events with < 3 observations ===== */
    while (1) {
        nobs = 0;
        nobswithw0 = 0;
        ntobs = 0;
        wsum = 0.0;
        
        /* Select input file (earthquakes vs shots) */
        if (i <= neqs) {
            iunit = 8;  /* EQS */
            fp = iunit_fp[8];
        } else {
            iunit = 9;  /* SHOTS */
            fp = iunit_fp[9];
        }
        
        /* Initialize observation arrays */
        xlat = 0.0;
        xlon = 0.0;
        for (jj = 0; jj < 180; jj++) {
            cphase[jj] = ' ';
        }
        
        i1 = 0;  /* dummy i1 switch */
        
        /* ===== CALL INPUT ROUTINE (based on data format) ===== */
        fprintf(stderr, "DEBUG ised=%d\n", ised);
        if (ised == 0) {
            inputcnv(fp, &nobsread,
                      sta_local, &iyr[i-1], &imo[i-1], &iday[i-1], 
                      &ihr[i-1], &imin[i-1], sec_local,
                      rmk1, rmk2, cphase, ipwt_local, amx, prx,
                      &xlat, &alon, &emag[i-1], &depth, &e[0][i-1],
                      &i1, &i2, &i3, &eventtype);
        } else if (ised == 1) {
            inputarcvel(fp, &nobsread,
                        sta_local, &iyr[i-1], &imo[i-1], &iday[i-1],
                        &ihr[i-1], &imin[i-1], sec_local,
                        rmk1, rmk2, cphase, ipwt_local, amx, prx,
                        &xlat, &alon, &emag[i-1], &depth, &e[0][i-1],
                        &i1, &i2, &i3, &eventtype);
        } else if (ised == 2) {
            inputsed_new(fp, &nobsread,
                         sta_local, &iyr[i-1], &imo[i-1], &iday[i-1],
                         &ihr[i-1], &imin[i-1], sec_local,
                         rmk1, rmk2, cphase, ipwt_local, amx, prx,
                         &xlat, &alon, &emag[i-1], &depth, &e[0][i-1],
                         &i1, &i2, &i3, &eventtype, itrial);
        } else {
            fprintf(stderr, "INPUTDATA: ised flag has no supported value!\n");
            exit(1);
        }
        
        /* Check for end-of-file */
        if (nobsread == -1) {
            knobs[i-1] = -1;
            return;
        }
        
        /* ===== Check for events with < 3 observations ===== */
        if (isingle != 0) {
            if (nobsread < 3) {
                printf("Event #%d\n", i);
                printf(" skipped because it has fewer than 3 obs.\n");
                continue;  /* Re-enter the while loop */
            }
        }
        
        if (isingle != 0) {
            fprintf(fm_ptr, "1 E V E N T   N R .   %6d                    0                    0\n", isingle);
            printf("Event #%d\n", isingle);
        }
        
        /* ===== Trial hypocenter from first station (if needed) ===== */
        if (itrial > 0 || (xlat == 0.0 && alon == 0.0)) {
            /* Find first station with weight < 5 on station list */
            if (itrial > 0) {
                jjmin1 = 0;
                for (jj = 0; jj < nobsread; jj++) {
                    if (ipwt_local[jj] < 5) {
                        int found = 0;
                        for (jk = 0; jk < nsta; jk++) {
                            if (strcmp(sta_local[jj], stn[jk]) == 0) {
                                jjmin1 = jj;
                                found = 1;
                                break;
                            }
                        }
                        if (found) break;
                    }
                }
                
                if (jjmin1 >= 0 && jjmin1 < nobsread) {
                    jjmin = jjmin1;
                    /* Find earliest time among valid stations */
                    for (jj = jjmin1; jj < nobsread; jj++) {
                        if (ipwt_local[jj] < 5) {
                            if (sec_local[jj] < sec_local[jjmin1]) {
                                jjmin = jj;
                            }
                        }
                    }
                    
                    /* Find station number for first phase */
                    jjmin1 = 0;
                    for (jj = 0; jj < nsta; jj++) {
                        if (strcmp(sta_local[jjmin], stn[jj]) == 0) {
                            jjmin1 = jj;
                            break;
                        }
                    }
                }
            }
            
            /* Set trial hypocenter */
            if (ifixsolution != 9) {
                xlat = xla[jjmin1] + 0.001;
                alon = xlo[jjmin1] + 0.001;
                if (ifixsolution != 1) {
                    depth = ztrial;
                } else {
                    if (!single_turbo) fprintf(fm_ptr, "DEPTH fixed !!!\n");
                    printf("DEPTH fixed !\n");
                }
            } else {
                if (!single_turbo) fprintf(fm_ptr, "HYPOCENTER fixed !!!\n");
                printf("HYPOCENTER fixed !\n");
                if (icoordsystem == 2 && alon > 0.) alon = -alon;
            }
        }
        
        /* Adjust depth if fixed at surface */
        if (ifixsolution != 0 && itopo > 0 && depth <= 0.0) {
            depth = -3.0;
        }
        
        /* Set depth or shot parameter */
        if (i > neqs) {
            jshot = i - neqs;  /* SHOTS */
            map2[jshot-1] = icc;
            z = depth;
        } else {
            z = depth + zshift;  /* EQS */
        }
        
        /* ===== DEBUG: Output icoordsystem ===== */
        fprintf(stderr, "[DEBUG] icoordsystem = %d\n", icoordsystem);

        /* ===== Transform LAT/LON to Cartesian coordinates ===== */
        if (icoordsystem == 2) {
            geoko(&e[1][i-1], &e[2][i-1], xlat, -alon, -1);
            e[1][i-1] = -e[1][i-1];
        } else {
            sdc(&e[1][i-1], &e[2][i-1], xlat, alon, -1);
        }
        e[3][i-1] = z;
        depthsofinput[i-1] = z;
        
        /* ===== OBSERVATION LOOP ===== */
        for (j = 0; j < nobsread; j++) {
            /* Weight cutoff checks */
            if (isingle == 0) {
                if (ipwt_local[j] >= 4) {  /* Skip weights >= 4 */
                    continue;
                }
            } else {
                if (ipwt_local[j] >= 6) {  /* Skip weights >= 6 */
                    continue;
                }
            }
            
            /* Skip S-waves if single P-wave model (nsp==1) */
            if (nsp == 1 && cphase[j] == 's') continue;
            if (nsp == 1 && cphase[j] == 'S') continue;
            if (nsp == 1 && cphase[j] == '-') continue;  /* Skip S-P phases */
            
            /* ===== PHASE TYPE VALIDATION ===== */
            /* Complex nested checks from Fortran lines 195-228 */
            if (cphase[j] != 'p' && cphase[j] != 'P' && ipwt_local[j] < 5) {
                if (cphase[j] != 's' && cphase[j] != 'S') {
                    if (cphase[j] != '-') {
                        if (cphase[j] != 'm' && cphase[j] != 'M') {
                            /* Unknown phase type */
                            if (!single_turbo) {
                                fprintf(fm_ptr, "WARNING:\n");
                                fprintf(fm_ptr, "what phase is this ?  %c ???\n", cphase[j]);
                            }
                            printf("what phase is this ?  %c ???\n", cphase[j]);
                            printf("\n");
                            if (isingle > 0) {
                                fprintf(fm_ptr, " DELETED: %s unknown phase is: %c\n", 
                                        sta_local[j], cphase[j]);
                            }
                            continue;  /* Skip this observation */
                        } else {
                            /* Reflected P phase (m or M) */
                            if (ireflector == 0) {
                                if (!single_turbo) {
                                    fprintf(fm_ptr, "WARNING:\n");
                                    fprintf(fm_ptr, "subr. INPUTDATA >>> Phase is : %c\n", cphase[j]);
                                    fprintf(fm_ptr, "but ireflector is: %d\n", ireflector);
                                    fprintf(fm_ptr, "Phase therefore neglected !!\n");
                                }
                                printf("subr. INPUTDATA >>> Phase is : %c\n", cphase[j]);
                                printf("but ireflector is: %d\n", ireflector);
                                printf("Phase therefore neglected !!\n");
                                printf("\n");
                                continue;  /* Skip this observation */
                            }
                        }
                    }
                }
            }
            
            /* ===== STATION NAME VALIDATION ===== */
            k = -1;  /* station index not found */
            for (k = 0; k < nsta; k++) {
                if (strcmp(sta_local[j], stn[k]) == 0) {
                    break;
                }
            }
            
            if (k >= nsta) {  /* Station not found */
                if (!single_turbo) {
                    if (isingle == 0) fprintf(fm_ptr, "WARNING:     Event # %d\n", i);
                    if (isingle > 0) fprintf(fm_ptr, "WARNING:     Event # %d\n", isingle);
                    fprintf(fm_ptr, " WARNING:  Station: >>>%s<<< not found in stationlist!\n", sta_local[j]);
                    fprintf(fm_ptr, "Phase therefore skipped\n");
                }
                printf(" WARNING:  Station: >>>%s<<< not found in stationlist!\n", sta_local[j]);
                printf("Phase therefore skipped\n");
                printf("\n");
                if (isingle > 0) {
                    fprintf(fm_ptr, " DELETED: %s not on station-list\n", sta_local[j]);
                }
                continue;  /* Skip this observation */
            }
            
            /* ===== DISTANCE CHECK ===== */
            ss1 = sqrt((x[k][0] - e[1][i-1]) * (x[k][0] - e[1][i-1]) +
                       (x[k][1] - e[2][i-1]) * (x[k][1] - e[2][i-1]));
            if (ss1 > dmax) {
                if (!single_turbo) {
                    fprintf(fm_ptr, "WARNING:\n");
                    fprintf(fm_ptr, " epicentral distance: %6.1f > dmax (%6.1f) ==> skipping phase !\n", ss1, dmax);
                }
                if (isingle != 0) {
                    printf(" epicentral distance: %6.1f > dmax (%6.1f) ==> skipping phase !\n", ss1, dmax);
                    fprintf(fm_ptr, " DELETED: %s epicentral-distance too large\n", sta_local[j]);
                }
                continue;  /* Skip this observation */
            }
            
            /* ===== DUPLICATE PHASE CHECK ===== */
            if (ipwt_local[j] != 5) {  /* Skip weight 5 phases from duplicate test */
                for (ll = 0; ll < nobs; ll++) {
                    itest = 9;  /* Default: no match */
                    
                    /* Skip if this is S-phase under nsp==1 */
                    if (nsp == 1 && (cphase[j] == 's' || cphase[j] == 'S')) {
                        continue;
                    }
                    
                    /* Check for duplicate */
                    if (k == istm[ll][i-1] && kpwt[ll][i-1] < 4) {
                        if (cphase[j] == 'p') itest = 0;
                        if (cphase[j] == 'P') itest = 0;
                        if (cphase[j] == 's') itest = 1;
                        if (cphase[j] == 'S') itest = 1;
                        if (cphase[j] == 'm') itest = -1;
                        if (cphase[j] == 'M') itest = -1;
                        if (cphase[j] == '-') itest = 2;
                        
                        if ((int)sphase[ll][i-1] == itest) {  /* Same phase type found */
                            if (isingle == 0) ie = i;
                            if (isingle > 0) ie = isingle;
                            
                            if (!single_turbo) {
                                fprintf(fm_ptr, "WARNING:\n");
                                fprintf(fm_ptr, "PHASETEST: POSSIBLE ERROR in phaselist !!\n");
                                fprintf(fm_ptr, "---> %02d %02d %02d %02d %02d\n",
                                        iyr[i-1], imo[i-1], iday[i-1], ihr[i-1], imin[i-1]);
                                fprintf(fm_ptr, "Event=%d Obs-nr. = %d >>> Station %s & Phase = %c already occured!\n",
                                        ie, j, sta_local[j], cphase[j]);
                            }
                            printf("PHASETEST: POSSIBLE ERROR in phaselist !!\n");
                            printf("---> %02d %02d %02d %02d %02d\n",
                                   iyr[i-1], imo[i-1], iday[i-1], ihr[i-1], imin[i-1]);
                            printf("Event=%d Obs-nr. = %d >>> Station %s & Phase = %c already occured!\n",
                                   ie, j, sta_local[j], cphase[j]);
                            printf("\n");
                            
                            if (isingle == 0) {
                                if (!single_turbo) {
                                    fprintf(fm_ptr, "subr. INPUTDATA >>> program will stop\n");
                                    fprintf(fm_ptr, "\n");
                                }
                                iphaseteststopflag = 1;
                            } else {
                                if (!single_turbo) {
                                    fprintf(fm_ptr, "INPUTDATA>>> nevertheless, program continues\n");
                                    fprintf(fm_ptr, "\n");
                                }
                            }
                        }
                    }
                }
            }
            
            /* ===== ACCEPT OBSERVATION (label 888) ===== */
            nobs++;
            
            /* Copy amplitude and period readings */
            amx[nobs-1] = amx[j];
            prx[nobs-1] = prx[j];
            
            if (isingle > 0) {
                if (ised == 2) {
                    prmk[nobs-1][0] = rmk1[j];
                    prmk[nobs-1][1] = rmk2[j];
                } else {
                    prmk[nobs-1][0] = ' ';
                    prmk[nobs-1][1] = ' ';
                }
            }
            
            /* Store station coordinates */
            for (l = 0; l < 3; l++) {
                d[nobs-1][l][i-1] = x[k][l];
            }
            
            /* Store observation data */
            pt[nobs-1][i-1] = sec_local[j];
            kpwt[nobs-1][i-1] = ipwt_local[j];
            istm[nobs-1][i-1] = k;
            strncpy(smn[nobs-1][i-1], sta_local[j], 4);
            smn[nobs-1][i-1][4] = '\0';
            
            /* ===== Handle P and S models (nsp==1 vs nsp==2) ===== */
            if (nsp != 2) {  /* nsp == 1: single velocity model */
                iphase[nobs-1][i-1] = model[0];
                
                /* Set phase type code */
                sphase[nobs-1][i-1] = 0.0;  /* Default: P-phase */
                if (cphase[j] == 's' || cphase[j] == 'S') {
                    sphase[nobs-1][i-1] = 1.0;  /* S-phase */
                }
                if (cphase[j] == 'm' || cphase[j] == 'M') {
                    sphase[nobs-1][i-1] = -1.0;  /* Reflected P */
                }
                if (cphase[j] == '-') {
                    sphase[nobs-1][i-1] = 2.0;  /* S-P */
                }
                
                /* Calculate weight */
                w[nobs-1][i-1] = 1.0 / (double)pow(2.0, 2.0 * ipwt_local[j]);
                
                if (cphase[j] == 's' || cphase[j] == 'S' || cphase[j] == '-') {
                    w[nobs-1][i-1] = swtfac * w[nobs-1][i-1];
                }
                
                if (ipwt_local[j] > 4) {
                    w[nobs-1][i-1] = 0.0;  /* Necessary for single-event mode */
                    nobswithw0++;
                }
                
                wsum += w[nobs-1][i-1];
            } else {
                /* nsp == 2: separate P and S velocity models */
                if (cphase[j] == 's' || cphase[j] == 'S' || cphase[j] == '-') {
                    /* S-phase or S-P phase */
                    sphase[nobs-1][i-1] = 1.0;
                    if (cphase[j] == '-') {
                        sphase[nobs-1][i-1] = 2.0;  /* S-P phase */
                    }
                    iphase[nobs-1][i-1] = model[1];  /* S-wave model */
                    w[nobs-1][i-1] = swtfac * 1.0 / (double)pow(2.0, 2.0 * ipwt_local[j]);
                } else {
                    /* P-phase */
                    sphase[nobs-1][i-1] = 0.0;
                    if (cphase[j] == 'm' || cphase[j] == 'M') {
                        sphase[nobs-1][i-1] = -1.0;  /* Reflected P */
                    }
                    iphase[nobs-1][i-1] = model[0];  /* P-wave model */
                    w[nobs-1][i-1] = 1.0 / (double)pow(2.0, 2.0 * ipwt_local[j]);
                }
                
                if (ipwt_local[j] >= 4) {
                    w[nobs-1][i-1] = 0.0;  /* Necessary for single-event mode */
                    nobswithw0++;
                }
                
                wsum += w[nobs-1][i-1];
            }
        }  /* End of observation loop (j) */
        
        /* ===== EVENT FINALIZATION ===== */
        knobs[i-1] = nobs;
        
        /* Count P and S observations separately (for nsp==2) */
        if (nsp == 2) {
            nobsp = 0;
            nobss = 0;
            for (j = 0; j < nobs; j++) {
                if (sphase[j][i-1] == 0.0) nobsp++;
                if (sphase[j][i-1] == 1.0) nobss++;
            }
        }
        
        /* ===== PRINT EVENT SUMMARY ===== */
        if (!single_turbo) {
            if (xlat < 0.0) {
                cns = 'S';
                xxlat = -xlat;
            } else {
                cns = 'N';
                xxlat = xlat;
            }
            if (alon < 0.0) {
                cew = 'E';
                aalon = -alon;
            } else {
                cew = 'W';
                aalon = alon;
            }
            if (icoordsystem == 2) {
                xxx = -e[1][i-1];
            } else {
                xxx = e[1][i-1];
            }
            
            if (i <= neqs) {  /* EQS */
                if (nsp == 2) {
                    fprintf(fm_ptr, "%3d %02d %02d %02d %02d %02d %5.2f %7.4f%c %8.4f%c %6.2f %7.2f %7.2f %7.2f %5.2f %2d %4d    %5d %5d\n",
                            i, iyr[i-1], imo[i-1], iday[i-1], ihr[i-1], imin[i-1], e[0][i-1],
                            xxlat, cns, aalon, cew, (double)depthsofinput[i-1], xxx, e[2][i-1], e[3][i-1],
                            emag[i-1], ifx[i-1], knobs[i-1] - nobswithw0, nobsp, nobss);
                } else {
                    fprintf(fm_ptr, "%3d %02d %02d %02d %02d %02d %5.2f %7.4f%c %8.4f%c %6.2f %7.2f %7.2f %7.2f %5.2f %2d %4d\n",
                            i, iyr[i-1], imo[i-1], iday[i-1], ihr[i-1], imin[i-1], e[0][i-1],
                            xxlat, cns, aalon, cew, (double)depthsofinput[i-1], xxx, e[2][i-1], e[3][i-1],
                            emag[i-1], ifx[i-1], knobs[i-1] - nobswithw0);
                }
            } else {  /* SHOTS */
                fprintf(fm_ptr, "%3d %02d %02d %02d %02d %02d %5.2f %7.4f%c %8.4f%c %6.2f %7.2f %7.2f %7.2f %5.2f %2d %4d\n",
                        i, iyr[i-1], imo[i-1], iday[i-1], ihr[i-1], imin[i-1], e[0][i-1],
                        xxlat, cns, aalon, cew, (double)depthsofinput[i-1], xxx, e[2][i-1], e[3][i-1],
                        emag[i-1], map2[i-neqs-1], knobs[i-1] - nobswithw0);
            }
        }
        
        /* ===== WEIGHT NORMALIZATION ===== */
        if (isingle > 0) {
            if ((knobs[i-1] - nobswithw0) < nvar) {
                return;
            }
            if (wsum <= 0.0) {
                return;
            }
        }
        
        for (j = 0; j < nobs; j++) {
            w[j][i-1] = w[j][i-1] * ((double)(nobs - nobswithw0) / wsum);
        }
        
        /* ===== CHECK FOR PHASE TEST ERROR AT end ===== */
        if (iphaseteststopflag == 1 && i == (neqs + nshot)) {
            fprintf(stderr, "INPUTDATA >>> PHASE-TEST: FATAL ERROR in phaselist !\n");
            exit(1);
        }
        
        /* ===== EXIT FROM RE-ENTRY LOOP ===== */
        break;  /* Only executed when event is accepted (>= 3 obs) */
    }  /* End of while(1) re-entry loop */
    
    return;
}
