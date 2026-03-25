#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

void lsvdb(double d[], double e[], int n,
           double v[][100], int iv, int nrv,
           double c[][100], int ic, int ncc, int *ier) {
    (void)iv;
    (void)ic;

    *ier = 0;
    if (n <= 0) return;

    // 构造原 bidiagonal 矩阵 B
    double *B = (double*)calloc(n * n, sizeof(double));
    if (!B) {
        *ier = 100; // 内存分配失败
        return;
    }

    // B[i*n + j] 为行主序
    for (int i = 0; i < n; ++i) {
        B[i*n + i] = d[i]; // 对角线
        if (i < n-1) B[i*n + i+1] = e[i+1]; // 次对角线
    }

    // U 和 VT
    double *U = (double*)calloc(n * n, sizeof(double));
    double *VT = (double*)calloc(n * n, sizeof(double));
    if (!U || !VT) {
        *ier = 101;
        free(B);
        if (U) free(U);
        if (VT) free(VT);
        return;
    }

    // 工作数组
    double *superb = (double*)malloc(sizeof(double) * (n-1));
    if (!superb) {
        *ier = 102;
        free(B); free(U); free(VT);
        return;
    }

    // 调用 LAPACKE_dgesvd
    int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', n, n,
                              B, n, d, U, n, VT, n, superb);
    free(superb);

    if (info != 0) {
        *ier = info;
        free(B); free(U); free(VT);
        return;
    }

    // 复制左奇异向量到 v[nrv][100]
    if (nrv > 0) {
        for (int i = 0; i < nrv; ++i)
            for (int j = 0; j < n; ++j)
                v[i][j] = U[i*n + j];
    }

    // 复制右奇异向量到 c[ncc][100]
    if (ncc > 0) {
        for (int i = 0; i < ncc; ++i)
            for (int j = 0; j < n; ++j)
                c[i][j] = VT[i*n + j];
    }

    free(B); free(U); free(VT);

    // 输出 d[] 已经是奇异值（LAPACK 默认升序，可能需要降序）
    // 将 d[] 降序排列，同时交换 v 和 c 对应列
    for (int i = 0; i < n-1; ++i) {
        int k = i;
        double t = d[i];
        for (int j = i+1; j < n; ++j) {
            if (d[j] > t) {
                t = d[j];
                k = j;
            }
        }
        if (k != i) {
            // 交换奇异值
            d[k] = d[i];
            d[i] = t;

            // 交换左奇异向量
            if (nrv > 0) {
                for (int r = 0; r < nrv; ++r) {
                    double tmp = v[r][i];
                    v[r][i] = v[r][k];
                    v[r][k] = tmp;
                }
            }

            // 交换右奇异向量
            if (ncc > 0) {
                for (int r = 0; r < ncc; ++r) {
                    double tmp = c[r][i];
                    c[r][i] = c[r][k];
                    c[r][k] = tmp;
                }
            }
        }
    }

    *ier = 33; // 成功
}