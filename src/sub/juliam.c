int juliam(int iyr, int imo, int idy, int ihr, int imn) {
    static const int kmo[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

    int ky = iyr;
    int km = (imo <= 0) ? 1 : imo;
    int kd = kmo[km - 1] + idy;

    int ky4 = ky / 4;
    int ky1 = ky / 100;
    int ky0 = ky / 1000;
    int kl = (ky4 - ky1 + ky0);

    int leap = (ky4 * 4 == ky && (ky1 * 100 != ky || ky0 * 1000 == ky)) ? 1 : 0;
    if (leap && km < 3) {
        kl -= 1;
    }

    int out = 365 * ky + kd + kl;
    out = out * 24 + ihr;
    out = out * 60 + imn;
    return out;
}
