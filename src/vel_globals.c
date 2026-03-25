#include <stdbool.h>
#include <stdio.h>

#define IEQ 658
#define INSHOT 50
#define ITOTMODELS 2
#define INLTOT 100
#define IST 650
#define MAXOBSPEREVENT 180
#define ILIP (4 * IEQ + INSHOT + INLTOT + IST)
#define INVA (ILIP - 1)
#define IKVA ((INVA * ILIP) / 2)

int nvar = 0, nvareff = 0, kvar = 0, nsta = 0, ksta = 0, legs = 0, lip = 0, neqs = 0;
int nsinv = 0, nshcor = 0, nshfix = 0, icount = 0, nshot = 0;

int ittmax = 0, invertratio = 0, iresolcalc = 0, ifixsolution = 0;
double zmin = 0.0f, zmininput = 0.0f, delmin = 0.0f, veladj = 0.0f, swtfac = 0.0f, xythet = 0.0f;
double zadj = 0.0f, vthet = 0.0f, othet = 0.0f, zthet = 0.0f, stathet = 0.0f, rmsmin = 0.0f, dmax = 0.0f;

int icnvout = 0, istaout = 0, ismpout = 0, irayout = 0, idrvout = 0, ialeout = 0, idspout = 0;
int irflout = 0, irfrout = 0, iresout = 0;

double zshift = 0.0f;

bool single_turbo = false;
int icoordsystem = 0, itrial = 0, ised = 0, isingle = 0, itopo = 0, iturbo = 0, lowveloclay = 0;
int ielev[IST] = {0}, iuseelev = 0, iusestacorr = 0;
double ztrial = 0.0f;

int nmod = 0, nsp = 0, nltot = 0, nplay[INLTOT] = {0}, laysum[INLTOT] = {0};
double v[INLTOT] = {0.0f}, vsq[INLTOT] = {0.0f}, h[INLTOT] = {0.0f}, thk[INLTOT] = {0.0f}, tkj = 0.0f, delta = 0.0f;
int jl = 0;

double vp[ITOTMODELS][INLTOT] = {{0.0f}};
double hp[ITOTMODELS][INLTOT] = {{0.0f}};
double vdamp[ITOTMODELS][INLTOT] = {{0.0f}};
double thkp[ITOTMODELS][INLTOT] = {{0.0f}};
double vpvs = 0.0f;
int ireflector = 0, lmax = 0;
char reflchar = ' ';

double ptcor[IST] = {0.0f}, stcor[IST] = {0.0f};
double d[IST][3][IEQ] = {{{0.0f}}};
int model[ITOTMODELS * IST] = {0};
double x[IST][3] = {{0.0f}};
int map1[IST] = {0}, map2[INSHOT] = {0};
double xla[IST] = {0.0f}, xlo[IST] = {0.0f};

char prmk[IST][2] = {{0}};
double pt[IST][IEQ] = {{0.0f}}, w[IST][IEQ] = {{0.0f}}, sphase[IST][IEQ] = {{0.0f}};
double depthsofinput[IEQ] = {0.0f}, tcalc[IST] = {0.0f}, amx[IST] = {0.0f}, prx[IST] = {0.0f};
int kpwt[IST][IEQ] = {{0}}, istm[IST][IEQ] = {{0}}, iphase[IST][IEQ] = {{0}}, knobs[IEQ] = {0};
int iain[IST] = {0}, nobswithw0 = 0, idelta[IST] = {0}, nactualsta[IST] = {0}, nstaeff = 0;
double tctime[IST][IEQ] = {{0.0f}};

char smn[IST][IEQ][5] = {{{0}}}, stn[IST][5] = {{0}}, blank[5] = "    ", blank0[5] = "0   ";

char regionname[33] = {0};
int iyr[IEQ] = {0}, imo[IEQ] = {0}, iday[IEQ] = {0}, ihr[IEQ] = {0}, imin[IEQ] = {0}, ifx[IEQ] = {0};
int isconstrain[3] = {0}, iconstrain[IEQ] = {0}, nmag = 0, nreg = 0;
double e[5][IEQ] = {{0.0f}}, emag[IEQ] = {0.0f}, xmagni[IST] = {0.0f}, xmagnitude = 0.0f, sdxmagnitude = 0.0f;
int igap[IEQ] = {0};

double rht[INVA * 2] = {0.0f}, res[IST][IEQ] = {{0.0f}}, rms[IEQ] = {0.0f}, b[INVA * 2] = {0.0f}, s[INVA * 2] = {0.0f};
double avres[IEQ] = {0.0f}, steplen = 0.0f;
char smpline[81] = {0};

double scale[7] = {0.0f};
double g[IKVA * 2] = {0.0f};
double gcopy[60][INVA * 2] = {{0.0f}};
double dtdv[INLTOT] = {0.0f}, dtdr[3] = {0.0f};
double gg2[MAXOBSPEREVENT][4] = {{0.0f}}, ggt[4][MAXOBSPEREVENT] = {{0.0f}}, ggg[4][MAXOBSPEREVENT] = {{0.0f}};
double gtg[4][4] = {{0.0f}}, ggti[4][4] = {{0.0f}}, drm[MAXOBSPEREVENT][MAXOBSPEREVENT] = {{0.0f}};

double rc[4][4] = {{0.0f}}, rs[4][4] = {{0.0f}}, sv[4] = {0.0f};
double covc[4][4] = {{0.0f}}, covs[4] = {0.0f}, ale[IEQ] = {0.0f}, r3[3][3] = {{0.0f}};
double gg[100][100] = {{0.0f}};

int irefrlayer[INLTOT] = {0}, ihypoclayer[INLTOT] = {0}, noheadwave = 0, irefllayer[INLTOT] = {0};
double refraylen[INLTOT] = {0.0f}, avhraylen = 0.0f, avvraylen = 0.0f, hitlay[INLTOT][3] = {{0.0f}};
double sterr = 0.0f, direrr = 0.0f, refrerr = 0.0f, reflerr = 0.0f;

double avrefrres = 0.0f, avotheres = 0.0f, avreflres = 0.0f, abrefrres = 0.0f, abotheres = 0.0f, abreflres = 0.0f;
int nrrefrres = 0, nrotheres = 0, nrreflres = 0;
double stnazires[IST][8] = {{0.0f}};

int ibackups = 0;
double davar1 = 0.0f, xmsqrs1 = 0.0f;
double spread = 0.0f;

char modelfilename[81] = {0}, stationfilename[81] = {0}, seismofilename[81] = {0}, scratchfilename[81] = {0};
char phasefile[81] = {0}, shotfile[81] = {0}, topo1file[81] = {0}, topo2file[81] = {0}, regnamfile[81] = {0}, regkoordfile[81] = {0};
char velfile[81] = {0}, cnvfile[81] = {0}, rayfile[81] = {0}, outfile[81] = {0}, smpfile[81] = {0}, stafile[81] = {0};
char drvfile[81] = {0}, alefile[81] = {0}, dsprfile[81] = {0}, rflfile[81] = {0}, rfrfile[81] = {0}, resfile[81] = {0};

char fm[81] = {0}, headerline[3][81] = {{0}};
int iabort = 0, nitt = 0, nittc = 0, istopflag = 0;

double olat = 0.0f, olon = 0.0f;
double rearth = 0.0, ellip = 0.0, rlatc = 0.0, rad = 0.0;
double aa = 0.0f, bb = 0.0f, bc = 0.0f, sint = 0.0f, cost = 0.0f, rotate = 0.0f;

int indfe[10000] = {0}, lt50[10000] = {0}, lt25[10000] = {0};
char irname[400000] = {0};
double xnrlon[200000] = {0.0f};

int inrpmax = (2 * INLTOT + 2);
int itotmodels = ITOTMODELS;
int inltot = INLTOT;

FILE *file77 = NULL, *file78 = NULL, *file79 = NULL;
/* File pointers for I/O (iunit 8=earthquakes, 9=shots) */
FILE *iunit_fp[10] = {NULL};

/* File pointer for report output (unit 16) */
FILE *fm_ptr = NULL;
FILE *fm_cnvfile = NULL;