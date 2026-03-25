int maxii(int n, const int *nx, int *imax, int *jndex) {
    int i;

    if (n <= 0 || nx == 0 || imax == 0 || jndex == 0) {
        return -1;
    }

    *jndex = 0;
    for (i = 0; i < n; ++i) {
        if (nx[*jndex] <= nx[i]) {
            *jndex = i;
        }
    }

    *imax = nx[*jndex];
    return 0;
}
