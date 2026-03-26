void ebell(double yl, double xb, double *l, double *b, double *my);
void elleb(double l, double b, double *x, double *y);

void geoko(double *x, double *y, double xlat, double xlon, int i) {
    double seichmy = 0.0;
    double lat = 0.0;
    double lon = 0.0;

    if (!x || !y) {
        return;
    }
    if (i != 1 && i != -1) {
        return;
    }

    if (i == 1) {
        ebell(xlat, xlon, &lon, &lat, &seichmy);
        *x = lat;
        *y = lon;
    } else {
        elleb(xlon, xlat, x, y);
    }
}
