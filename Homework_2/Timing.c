#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Matrix.h"

// Function prototype for time_matvec
double time_matvec(void (*matvec_func)(double *, const Matrix *, const double *),
                   const Matrix *A, const double *x, int num_samples);

int main(int argc, char *argv[]) {
    int m, n, num_samples;

    // Check if command line arguments are provided
    if (argc > 3) {
        m = atoi(argv[1]);
        n = atoi(argv[2]);
        num_samples = atoi(argv[3]);
    } else {
        m = n = 1000;
        num_samples = 10;
    }

    // Initialize matrices and vectors
    Matrix A_cont = init_matrix_contiguous(m, n);
    Matrix A_non_cont = init_matrix_non_contiguous(m, n);
    double *x = (double *)malloc(sizeof(double) * n);

    // For simplicity, filling the matrices and vector with random values
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            A_cont.data[i][j] = rand() % 10;
            A_non_cont.data[i][j] = rand() % 10;
        }
    }
    for (int j = 0; j < n; j++) {
        x[j] = rand() % 10;
    }

    // Timing evaluations
    printf("Average time for matvec (contiguous): %f seconds\n", 
           time_matvec(matvec, &A_cont, x, num_samples));
    printf("Average time for matvec (non-contiguous): %f seconds\n", 
           time_matvec(matvec, &A_non_cont, x, num_samples));
    printf("Average time for matvec_transpose (contiguous): %f seconds\n", 
           time_matvec(matvec_transpose, &A_cont, x, num_samples));
    printf("Average time for matvec_transpose (non-contiguous): %f seconds\n", 
           time_matvec(matvec_transpose, &A_non_cont, x, num_samples));

    // Clean up
    free_matrix_contiguous(&A_cont);
    free_matrix_non_contiguous(&A_non_cont);
    free(x);

    return 0;
}

double time_matvec(void (*matvec_func)(double *, const Matrix *, const double *),
                   const Matrix *A, const double *x, int num_samples) {
    double *b = (double *)malloc(sizeof(double) * A->m);
    clock_t start, end;
    double total_time = 0;

    for (int i = 0; i < num_samples; i++) {
        start = clock();
        matvec_func(b, A, x);
        end = clock();
        total_time += (double)(end - start) / CLOCKS_PER_SEC;
    }

    free(b);
    return total_time / num_samples;
}
