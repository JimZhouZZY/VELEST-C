void rside(double *rhs, const double *s, int n, const double *res,
           const double *w) {
    int i;
    double rw;

    if (rhs == 0 || s == 0 || res == 0 || w == 0 || n <= 0) {
        return;
    }

    rw = (*res) * (*w);

    for (i = 0; i < n; ++i) {
        if (s[i] == 0.0) {
            continue;
        }
        rhs[i] += rw * s[i];
    }
}
