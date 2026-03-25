int maxri(int n, const double *x, double *xmax, int *jndex) {
    int i;

    if (n <= 0 || x == 0 || xmax == 0 || jndex == 0) {
        return -1;
    }

    *jndex = 0;
    for (i = 0; i < n; ++i) {
        if (x[*jndex] <= x[i]) {
            *jndex = i;
        }
    }

    *xmax = x[*jndex];
    return 0;
}
