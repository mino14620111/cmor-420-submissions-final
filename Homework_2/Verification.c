#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Assuming the struct and function declarations are in a file named "matrix_functions.h"
#include "Matrix.h"

double compute_error(double* b_exact, double* b_computed, int length) {
    double error = 0.0;
    for (int i = 0; i < length; i++) {
        error += fabs(b_exact[i] - b_computed[i]);
    }
    return error;
}


int main() {
    int m = 2;
    int n = 2;

    // Pre-defined matrix and vector for easy verification
    double A_values[2][2] = {
        {2.0, 3.0},
        {4.0, 5.0}
    };
    double x_values[2] = {1.0, 2.0};

    // Expected results
    double b_exact[2] = {8.0, 14.0};
    double b_exact_transpose[2] = {10.0, 13.0};

    // Initialize matrices and vectors
    Matrix A_cont = init_matrix_contiguous(m, n);
    Matrix A_non_cont = init_matrix_non_contiguous(m, n);
    
    // Assign values to A
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            A_cont.data[i][j] = A_values[i][j];
            A_non_cont.data[i][j] = A_values[i][j];
        }
    }

    // Compute b = A*x and b = transpose(A)*x
    double b_cont[2] = {0};
    double b_non_cont[2] = {0};
    double b_cont_transpose[2] = {0};
    double b_non_cont_transpose[2] = {0};

    matvec(b_cont, &A_cont, x_values);
    matvec(b_non_cont, &A_non_cont, x_values);
    matvec_transpose(b_cont_transpose, &A_cont, x_values);
    matvec_transpose(b_non_cont_transpose, &A_non_cont, x_values);

    // Compute and print errors
    printf("The error in computing b=Ax for matvec with contiguous matrix storage is %.5e\n",
           compute_error(b_exact, b_cont, m));
    printf("The error in computing b=Ax for matvec with non-contiguous matrix storage is %.5e\n",
           compute_error(b_exact, b_non_cont, m));
    printf("The error in computing b=A^Tx for matvec_transpose with contiguous matrix storage is %.5e\n",
           compute_error(b_exact_transpose, b_cont_transpose, m));
    printf("The error in computing b=A^Tx for matvec_transpose with non-contiguous matrix storage is %.5e\n",
           compute_error(b_exact_transpose, b_non_cont_transpose, m));

    // Free allocated memory
    free_matrix_contiguous(&A_cont);
    free_matrix_non_contiguous(&A_non_cont);

    return 0;
}
