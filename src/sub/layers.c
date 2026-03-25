#include <math.h>
#include <stdlib.h>

void layers(int nx, int ny, int nz,
            const double x[], const double y[], const double z[],
            const double vel[][100][100],
            double xe, double ye, double xr, double yr,
            int nl, double v[], double vsq[]) {
    (void)z;
    (void)nz;

    double weight[4] = {0.0f};
    double slowness_sum[100] = {0.0f};

    int event_i = 0;
    while (event_i < nx && x[event_i] < xe) event_i++;
    event_i -= 1;
    if (event_i < 0) event_i = 0;

    int event_j = 0;
    while (event_j < ny && y[event_j] < ye) event_j++;
    event_j -= 1;
    if (event_j < 0) event_j = 0;

    int recv_i = 0;
    while (recv_i < nx && x[recv_i] < xr) recv_i++;
    recv_i -= 1;
    if (recv_i < 0) recv_i = 0;

    int recv_j = 0;
    while (recv_j < ny && y[recv_j] < yr) recv_j++;
    recv_j -= 1;
    if (recv_j < 0) recv_j = 0;

    int min_i = (event_i < recv_i) ? event_i : recv_i;
    int min_j = (event_j < recv_j) ? event_j : recv_j;

    int num_segments = abs(recv_i - event_i) + abs(recv_j - event_j) + 1;
    int num_points = num_segments + 1;

    double x_segment = (xr - xe) / (double)num_segments;
    double y_segment = (yr - ye) / (double)num_segments;

    for (int layer_index = 0; layer_index < nl; ++layer_index) {
        slowness_sum[layer_index] = 0.0f;
    }

    for (int point_index = 0; point_index < num_points; ++point_index) {
        double x_point = xe + point_index * x_segment;
        double y_point = ye + point_index * y_segment;

        int i_scan = min_i;
        while (i_scan < nx && x_point > x[i_scan]) i_scan++;
        int x_cell = i_scan - 1;
        if (x_cell < 0) x_cell = 0;
        if (x_cell >= nx - 1) x_cell = nx - 2;

        int j_scan = min_j;
        while (j_scan < ny && y_point > y[j_scan]) j_scan++;
        int y_cell = j_scan - 1;
        if (y_cell < 0) y_cell = 0;
        if (y_cell >= ny - 1) y_cell = ny - 2;

        double a = x_point - x[x_cell];
        double b = x[x_cell + 1] - x_point;
        double c = y_point - y[y_cell];
        double d = y[y_cell + 1] - y_point;
        double denominator = (a + b) * (c + d);

        weight[0] = b * d / denominator;
        weight[1] = a * d / denominator;
        weight[2] = b * c / denominator;
        weight[3] = a * c / denominator;

        for (int layer_index = 0; layer_index < nl; ++layer_index) {
            double velocity_point = weight[0] * vel[x_cell][y_cell][layer_index]
                                 + weight[1] * vel[x_cell + 1][y_cell][layer_index]
                                 + weight[2] * vel[x_cell][y_cell + 1][layer_index]
                                 + weight[3] * vel[x_cell + 1][y_cell + 1][layer_index];
            slowness_sum[layer_index] += 1.0f / velocity_point;
        }
    }

    for (int layer_index = 0; layer_index < nl; ++layer_index) {
        v[layer_index] = (double)num_points / slowness_sum[layer_index];
        vsq[layer_index] = v[layer_index] * v[layer_index];
    }
}
