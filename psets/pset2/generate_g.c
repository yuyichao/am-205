// compile using GCC 5.2.0 with the following command line
// gcc -fopenmp -fcilkplus -mavx2 -fPIC -Ofast generate_g.c -o libgen_g.so
//     -shared -Wall -Wextra

#include <stddef.h>
#include <cilk/cilk.h>

void
fill_cilk(double *__restrict__ ary, size_t size)
{
    cilk_for (size_t i = 0;i < size;i++) {
        ary[i] = 0;
    }
}

void
fill_omp(double *__restrict__ ary, size_t size)
{
#pragma omp parallel for
    for (size_t i = 0;i < size;i++) {
        ary[i] = 0;
    }
}

void
generate_g(double *__restrict__ ptr, size_t n)
{
    for (size_t i = 0;i < n - 1;i++) {
        double *ptr2 = ptr + i * n;
#pragma simd
        for (size_t j = 0;j < n;j++) {
            ptr2[j] = j < i ? 0 : (j == i ? 1 : -1);
        }
    }
#pragma simd
    for (size_t j = 0;j < n;j++) {
        ptr[(n - 1) * n + j] = 1;
    }
}

void
generate_g_cilk(double *__restrict__ ptr, size_t n)
{
    cilk_for (size_t i = 0;i < n - 1;i++) {
        double *ptr2 = ptr + i * n;
#pragma simd
        for (size_t j = 0;j < n;j++) {
            ptr2[j] = j < i ? 0 : (j == i ? 1 : -1);
        }
    }
    cilk_for (size_t j = 0;j < n;j++) {
        ptr[(n - 1) * n + j] = 1;
    }
}

void
generate_g_omp(double *__restrict__ ptr, size_t n)
{
#pragma omp parallel for
    for (size_t i = 0;i < n - 1;i++) {
        double *ptr2 = ptr + i * n;
#pragma simd
        for (size_t j = 0;j < n;j++) {
            ptr2[j] = j < i ? 0 : (j == i ? 1 : -1);
        }
    }
#pragma omp parallel for
    for (size_t j = 0;j < n;j++) {
        ptr[(n - 1) * n + j] = 1;
    }
}
