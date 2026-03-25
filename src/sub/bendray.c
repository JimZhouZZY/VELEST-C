#include <math.h>
#include <stdio.h>

extern int icoordsystem;
extern const char *topo1file;
extern const char *topo2file;
extern void chtop(double x, double y, double *z, const char *topo1, const char *topo2);

void bendray(double rp[3][200], int nrp, const char *staname, double vtop, double *ttt) {
    if (icoordsystem != 2) {
        return;
    }

    double xydist = sqrtf((rp[0][nrp - 1] - rp[0][0]) * (rp[0][nrp - 1] - rp[0][0]) +
                         (rp[1][nrp - 1] - rp[1][0]) * (rp[1][nrp - 1] - rp[1][0]));
    if (xydist > 10.0f) {
        return;
    }

    double takeoff_num = sqrtf((rp[0][1] - rp[0][0]) * (rp[0][1] - rp[0][0]) +
                              (rp[1][1] - rp[1][0]) * (rp[1][1] - rp[1][0]));
    double takeoff_den = sqrtf((rp[0][1] - rp[0][0]) * (rp[0][1] - rp[0][0]) +
                              (rp[1][1] - rp[1][0]) * (rp[1][1] - rp[1][0]) +
                              (rp[2][1] - rp[2][0]) * (rp[2][1] - rp[2][0]));
    double takeoff_angle = 57.296f * asinf(takeoff_num / takeoff_den);
    if ((rp[2][1] - rp[2][0]) < 0.0f) {
        takeoff_angle = 180.0f - takeoff_angle;
    }

    int nrpm1 = nrp - 2;
    double arrive_num = sqrtf((rp[0][nrp - 1] - rp[0][nrpm1]) * (rp[0][nrp - 1] - rp[0][nrpm1]) +
                             (rp[1][nrp - 1] - rp[1][nrpm1]) * (rp[1][nrp - 1] - rp[1][nrpm1]));
    double arrive_den = sqrtf((rp[0][nrp - 1] - rp[0][nrpm1]) * (rp[0][nrp - 1] - rp[0][nrpm1]) +
                             (rp[1][nrp - 1] - rp[1][nrpm1]) * (rp[1][nrp - 1] - rp[1][nrpm1]) +
                             (rp[2][nrp - 1] - rp[2][nrpm1]) * (rp[2][nrp - 1] - rp[2][nrpm1]));
    double arrive_angle = 57.296f * asinf(arrive_num / arrive_den);
    if ((rp[2][1] - rp[2][0]) < 0.0f) {
        arrive_angle = 180.0f - arrive_angle;
    }

    int hypocintop = 0;
    double deltat1 = 0.0f;
    double deltat2 = 0.0f;

    double dx1 = 0.0f, dy1 = 0.0f, dz1 = 0.0f, xyz1 = 0.0f;
    if (rp[2][0] < 0.0f) {
        hypocintop = 1;
        dx1 = rp[0][1] - rp[0][0];
        dy1 = rp[1][1] - rp[1][0];
        dz1 = rp[2][1] - rp[2][0];
        xyz1 = sqrtf(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
    }

    if (hypocintop && takeoff_angle <= 45.0f) {
        hypocintop = 0;
    }
    if (!hypocintop && arrive_angle >= 135.0f) {
        return;
    }

    double dx2 = rp[0][nrp - 1] - rp[0][nrpm1];
    double dy2 = rp[1][nrp - 1] - rp[1][nrpm1];
    double dz2 = rp[2][nrp - 1] - rp[2][nrpm1];
    double xyz2 = sqrtf(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);

    if (hypocintop) {
        dx1 /= 4.0f;
        dy1 /= 4.0f;
        dz1 /= 4.0f;
    }
    dx2 /= 4.0f;
    dy2 /= 4.0f;
    dz2 /= 4.0f;

    double rpn[3][10] = {{0.0f}};

    if (hypocintop) {
        for (int j = 0; j < 4; ++j) {
            rpn[0][j] = rp[0][0] + j * dx1;
            rpn[1][j] = rp[1][0] + j * dy1;
            rpn[2][j] = rp[2][0] + j * dz1;
        }
        rpn[0][4] = rp[0][1];
        rpn[1][4] = rp[1][1];
        rpn[2][4] = rp[2][1];
    }

    for (int j = 5; j <= 8; ++j) {
        rpn[0][j] = rp[0][nrpm1] + (j - 5) * dx2;
        rpn[1][j] = rp[1][nrpm1] + (j - 5) * dy2;
        rpn[2][j] = rp[2][nrpm1] + (j - 5) * dz2;
    }
    rpn[0][9] = rp[0][nrp - 1];
    rpn[1][9] = rp[1][nrp - 1];
    rpn[2][9] = rp[2][nrp - 1];

    if (hypocintop) {
        for (int j = 1; j <= 4; ++j) {
            double zzz = 0.0f;
            chtop(-rpn[0][j], rpn[1][j], &zzz, topo1file, topo2file);
            zzz += 0.1f;
            if (rpn[2][j] < zzz) {
                rpn[2][j] = zzz;
            }
        }
    }

    for (int j = 5; j <= 8; ++j) {
        double zzz = 0.0f;
        chtop(-rpn[0][j], rpn[1][j], &zzz, topo1file, topo2file);
        zzz += 0.1f;
        if (rpn[2][j] < zzz) {
            rpn[2][j] = zzz;
        }
    }

    double ttt1new = 0.0f;
    if (hypocintop) {
        for (int j = 1; j <= 4; ++j) {
            int jm1 = j - 1;
            double xyz1n = sqrtf((rpn[0][j] - rpn[0][jm1]) * (rpn[0][j] - rpn[0][jm1]) +
                               (rpn[1][j] - rpn[1][jm1]) * (rpn[1][j] - rpn[1][jm1]) +
                               (rpn[2][j] - rpn[2][jm1]) * (rpn[2][j] - rpn[2][jm1]));
            ttt1new += xyz1n / vtop;
        }
    }

    double ttt2new = 0.0f;
    for (int j = 6; j <= 9; ++j) {
        int jm1 = j - 1;
        double xyz2n = sqrtf((rpn[0][j] - rpn[0][jm1]) * (rpn[0][j] - rpn[0][jm1]) +
                           (rpn[1][j] - rpn[1][jm1]) * (rpn[1][j] - rpn[1][jm1]) +
                           (rpn[2][j] - rpn[2][jm1]) * (rpn[2][j] - rpn[2][jm1]));
        ttt2new += xyz2n / vtop;
    }

    if (hypocintop) {
        double ttt1old = xyz1 / vtop;
        deltat1 = ttt1new - ttt1old;
    }

    double ttt2old = xyz2 / vtop;
    deltat2 = ttt2new - ttt2old;

    if (fabs(deltat1) > 1e-5f || fabs(deltat2) > 1e-5f) {
        fprintf(stdout, " BENDRAY>>> ray bended below surface! Station: %.4s\n", staname);
        fprintf(stdout, " dt1=%6.3f   dt2=%6.3f\n", deltat1, deltat2);
    }

    *ttt += deltat1 + deltat2;
}
