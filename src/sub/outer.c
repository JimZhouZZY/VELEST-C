void outer(double *g, const double *s, int n, const double *w) {
    int i;
    int j;
    double weight;

    if (g == 0 || s == 0 || w == 0 || n <= 0) {
        return;
    }

    weight = *w;

    for (i = 0; i < n; ++i) {
        double a;
        if (s[i] == 0.0) {
            continue;
        }

        a = s[i] * weight;

        for (j = 0; j <= i; ++j) {
            int nsym;
            if (s[j] == 0.0) {
                continue;
            }

            nsym = (i * (i + 1)) / 2 + j;
            g[nsym] += s[j] * a;
        }
    }
}
