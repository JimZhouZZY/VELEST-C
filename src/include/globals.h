#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdio.h>
#include <stdbool.h>

/* ------------------     Constants    ------------------ */

#define IEQ 658
#define INSHOT 50
#define ITOTMODELS 2
#define INLTOT 100
#define IST 650
#define MAX_STLINES 600
#define MAXOBSPEREVENT 180
#define INRPMAX (2 * INLTOT + 2)
#define ILIP (4 * IEQ + INSHOT + INLTOT + IST)
#define INVA (ILIP - 1)
#define IKVA ((INVA * ILIP) / 2)

/* ------------------ Global Variables ------------------ */

/* integers */
extern int nvar, nvareff, kvar, nsta, ksta, legs, lip, neqs;
extern int nsinv, nshcor, nshfix, icount, nshot;
extern int ittmax, invertratio, iresolcalc, ifixsolution;
extern int icnvout, istaout, ismpout, irayout, idrvout, ialeout, idspout;
extern int irflout, irfrout, iresout;
extern int icoordsystem, itrial, ised, isingle, itopo, iturbo, lowveloclay;
extern int ielev[IST], iuseelev, iusestacorr;
extern int nmod, nsp, nltot, nplay[INLTOT], laysum[INLTOT];
extern int jl;
extern int ireflector, lmax;
extern int model[ITOTMODELS * IST];
extern int map1[IST], map2[INSHOT];
extern int kpwt[IST][IEQ], istm[IST][IEQ], iphase[IST][IEQ], knobs[IEQ];
extern int iain[IST], nobswithw0, idelta[IST], nactualsta[IST], nstaeff;
extern int iyr[IEQ], imo[IEQ], iday[IEQ], ihr[IEQ], imin[IEQ], ifx[IEQ];
extern int isconstrain[3], iconstrain[IEQ], nmag, nreg;
extern int igap[IEQ];
extern int irefrlayer[INLTOT], ihypoclayer[INLTOT], noheadwave, irefllayer[INLTOT];
extern int ibackups;
extern int iabort, nitt, nittc, istopflag;
extern int indfe[10000], lt50[10000], lt25[10000];
extern int inrpmax, itotmodels, inltot;

/* doubles */
extern double zmin, zmininput, delmin, veladj, swtfac, xythet;
extern double zadj, vthet, othet, zthet, stathet, rmsmin, dmax;
extern double zshift, ztrial;
extern double v[INLTOT], vsq[INLTOT], h[INLTOT], thk[INLTOT], tkj, delta;
extern double vp[ITOTMODELS][INLTOT], hp[ITOTMODELS][INLTOT], vdamp[ITOTMODELS][INLTOT], thkp[ITOTMODELS][INLTOT];
extern double vpvs;
extern double ptcor[IST], stcor[IST];
extern double d[IST][3][IEQ];
extern double x[IST][3];
extern double xla[IST], xlo[IST];
extern double pt[IST][IEQ], w[IST][IEQ], sphase[IST][IEQ];
extern double depthsofinput[IEQ], tcalc[IST], amx[IST], prx[IST];
extern double tctime[IST][IEQ];
extern double xmagni[IST], xmagnitude, sdxmagnitude;
extern double rht[INVA*2], res[IST][IEQ], rms[IEQ], b[INVA*2], s[INVA*2];
extern double avres[IEQ], steplen;
extern double scale[7], g[IKVA*2], gcopy[60][INVA*2];
extern double dtdv[INLTOT], dtdr[3];
extern double gg2[MAXOBSPEREVENT][4], ggt[4][MAXOBSPEREVENT], ggg[4][MAXOBSPEREVENT];
extern double gtg[4][4], ggti[4][4], drm[MAXOBSPEREVENT][MAXOBSPEREVENT];
extern double rc[4][4], rs[4][4], sv[4];
extern double covc[4][4], covs[4], ale[IEQ], r3[3][3], gg[100][100];
extern double refraylen[INLTOT], avhraylen, avvraylen, hitlay[INLTOT][3];
extern double sterr, direrr, refrerr, reflerr;
extern double avrefrres, avotheres, avreflres, abrefrres, abotheres, abreflres;
extern double stnazires[IST][8];
extern double davar1, xmsqrs1, spread;
extern double olat, olon;
extern double rearth, ellip, rlatc, rad;
extern double aa, bb, bc, sint, cost, rotate;
extern double xnrlon[200000];

/* chars */
extern char reflchar;
extern char prmk[IST][2];
extern char smn[IST][IEQ][5], stn[IST][5], blank[5], blank0[5];
extern char regionname[33];
extern char smpline[81];
extern char modelfilename[81], stationfilename[81], seismofilename[81], scratchfilename[81];
extern char phasefile[81], shotfile[81], topo1file[81], topo2file[81], regnamfile[81], regkoordfile[81];
extern char velfile[81], cnvfile[81], rayfile[81], outfile[81], smpfile[81], stafile[81];
extern char drvfile[81], alefile[81], dsprfile[81], rflfile[81], rfrfile[81], resfile[81];
extern char fm[81], headerline[3][81];
extern char irname[400000];

/* bools */
extern bool single_turbo;

/* FILE pointers */
extern FILE *file77, *file78, *file79;
extern FILE *iunit_fp[10];
extern FILE *fm_ptr, *fm_cnvfile;

#endif /* GLOBALS_H */