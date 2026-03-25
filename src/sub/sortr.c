void sortr(double *x, int n) {
    int gap;

    if (x == 0 || n <= 1) {
        return;
    }

    for (gap = n / 2; gap > 0; gap /= 2) {
        int i;
        for (i = gap; i < n; ++i) {
            double value = x[i];
            int j = i;
            while (j >= gap && x[j - gap] > value) {
                x[j] = x[j - gap];
                j -= gap;
            }
            x[j] = value;
        }
    }
}
