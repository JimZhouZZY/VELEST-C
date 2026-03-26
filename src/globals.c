#include <stdbool.h>
#include <stdio.h>
#include "include/globals.h"

int nvar = 0, nvareff = 0, kvar = 0, nsta = 0, ksta = 0, legs = 0, lip = 0, neqs = 0;
int nsinv = 0, nshcor = 0, nshfix = 0, icount = 0, nshot = 0;

int ittmax = 0, invertratio = 0, iresolcalc = 0, ifixsolution = 0;
double zmin = 0.0, zmininput = 0.0, delmin = 0.0, veladj = 0.0, swtfac = 0.0, xythet = 0.0;
double zadj = 0.0, vthet = 0.0, othet = 0.0, zthet = 0.0, stathet = 0.0, rmsmin = 0.0, dmax = 0.0;

int icnvout = 0, istaout = 0, ismpout = 0, irayout = 0, idrvout = 0, ialeout = 0, idspout = 0;
int irflout = 0, irfrout = 0, iresout = 0;

double zshift = 0.0;

bool single_turbo = false;
int icoordsystem = 0, itrial = 0, ised = 0, isingle = 0, itopo = 0, iturbo = 0, lowveloclay = 0;
int ielev[IST] = {0}, iuseelev = 0, iusestacorr = 0;
double ztrial = 0.0;

int nmod = 0, nsp = 0, nltot = 0, nplay[INLTOT] = {0}, laysum[INLTOT] = {0};
double v[INLTOT] = {0.0}, vsq[INLTOT] = {0.0}, h[INLTOT] = {0.0}, thk[INLTOT] = {0.0}, tkj = 0.0, delta = 0.0;
int jl = 0;

double vp[ITOTMODELS][INLTOT] = {{0.0}};
double hp[ITOTMODELS][INLTOT] = {{0.0}};
double vdamp[ITOTMODELS][INLTOT] = {{0.0}};
double thkp[ITOTMODELS][INLTOT] = {{0.0}};
double vpvs = 0.0;
int ireflector = 0, lmax = 0;
char reflchar = ' ';

double ptcor[IST] = {0.0}, stcor[IST] = {0.0};
double d[IST][3][IEQ] = {{{0.0}}};
int model[ITOTMODELS * IST] = {0};
double x[IST][3] = {{0.0}};
int map1[IST] = {0}, map2[INSHOT] = {0};
double xla[IST] = {0.0}, xlo[IST] = {0.0};

char prmk[IST][2] = {{0}};
double pt[IST][IEQ] = {{0.0}}, w[IST][IEQ] = {{0.0}}, sphase[IST][IEQ] = {{0.0}};
double depthsofinput[IEQ] = {0.0}, tcalc[IST] = {0.0}, amx[IST] = {0.0}, prx[IST] = {0.0};
int kpwt[IST][IEQ] = {{0}}, istm[IST][IEQ] = {{0}}, iphase[IST][IEQ] = {{0}}, knobs[IEQ] = {0};
int iain[IST] = {0}, nobswithw0 = 0, idelta[IST] = {0}, nactualsta[IST] = {0}, nstaeff = 0;
double tctime[IST][IEQ] = {{0.0}};

char smn[IST][IEQ][5] = {{{0}}}, stn[IST][5] = {{0}}, blank[5] = "    ", blank0[5] = "0   ";

char regionname[33] = {0};
int iyr[IEQ] = {0}, imo[IEQ] = {0}, iday[IEQ] = {0}, ihr[IEQ] = {0}, imin[IEQ] = {0}, ifx[IEQ] = {0};
int isconstrain[3] = {0}, iconstrain[IEQ] = {0}, nmag = 0, nreg = 0;
double e[5][IEQ] = {{0.0}}, emag[IEQ] = {0.0}, xmagni[IST] = {0.0}, xmagnitude = 0.0, sdxmagnitude = 0.0;
int igap[IEQ] = {0};

double rht[INVA * 2] = {0.0}, res[IST][IEQ] = {{0.0}}, rms[IEQ] = {0.0}, b[INVA * 2] = {0.0}, s[INVA * 2] = {0.0};
double avres[IEQ] = {0.0}, steplen = 0.0;
char smpline[81] = {0};

double scale[7] = {0.0};
double g[IKVA * 2] = {0.0};
double gcopy[60][INVA * 2] = {{0.0}};
double dtdv[INLTOT] = {0.0}, dtdr[3] = {0.0};
double gg2[MAXOBSPEREVENT][4] = {{0.0}}, ggt[4][MAXOBSPEREVENT] = {{0.0}}, ggg[4][MAXOBSPEREVENT] = {{0.0}};
double gtg[4][4] = {{0.0}}, ggti[4][4] = {{0.0}}, drm[MAXOBSPEREVENT][MAXOBSPEREVENT] = {{0.0}};

double rc[4][4] = {{0.0}}, rs[4][4] = {{0.0}}, sv[4] = {0.0};
double covc[4][4] = {{0.0}}, covs[4] = {0.0}, ale[IEQ] = {0.0}, r3[3][3] = {{0.0}};
double gg[100][100] = {{0.0}};

int irefrlayer[INLTOT] = {0}, ihypoclayer[INLTOT] = {0}, noheadwave = 0, irefllayer[INLTOT] = {0};
double refraylen[INLTOT] = {0.0}, avhraylen = 0.0, avvraylen = 0.0, hitlay[INLTOT][3] = {{0.0}};
double sterr = 0.0, direrr = 0.0, refrerr = 0.0, reflerr = 0.0;

double avrefrres = 0.0, avotheres = 0.0, avreflres = 0.0, abrefrres = 0.0, abotheres = 0.0, abreflres = 0.0;
int nrrefrres = 0, nrotheres = 0, nrreflres = 0;
double stnazires[IST][8] = {{0.0}};

int ibackups = 0;
double davar1 = 0.0, xmsqrs1 = 0.0;
double spread = 0.0;

char modelfilename[81] = {0}, stationfilename[81] = {0}, seismofilename[81] = {0}, scratchfilename[81] = {0};
char phasefile[81] = {0}, shotfile[81] = {0}, topo1file[81] = {0}, topo2file[81] = {0}, regnamfile[81] = {0}, regkoordfile[81] = {0};
char velfile[81] = {0}, cnvfile[81] = {0}, rayfile[81] = {0}, outfile[81] = {0}, smpfile[81] = {0}, stafile[81] = {0};
char drvfile[81] = {0}, alefile[81] = {0}, dsprfile[81] = {0}, rflfile[81] = {0}, rfrfile[81] = {0}, resfile[81] = {0};

char fm[81] = {0}, headerline[3][81] = {{0}};
int iabort = 0, nitt = 0, nittc = 0, istopflag = 0;

double olat = 0.0, olon = 0.0;
double rearth = 0.0, ellip = 0.0, rlatc = 0.0, rad = 0.0;
double aa = 0.0, bb = 0.0, bc = 0.0, sint = 0.0, cost = 0.0, rotate = 0.0;

int indfe[10000] = {0}, lt50[10000] = {0}, lt25[10000] = {0};
char irname[400000] = {0};
double xnrlon[200000] = {0.0};

int inrpmax = (2 * INLTOT + 2);
int itotmodels = ITOTMODELS;
int inltot = INLTOT;

FILE *file77 = NULL, *file78 = NULL, *file79 = NULL;
/* File pointers for I/O (iunit 8=earthquakes, 9=shots) */
FILE *iunit_fp[10] = {NULL};

/* File pointer for report output (unit 16) */
FILE *fm_ptr = NULL;
FILE *fm_cnvfile = NULL;