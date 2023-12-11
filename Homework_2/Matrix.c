#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int m;
    int n;
    double **data;
} Matrix;

Matrix init_matrix_contiguous(int m, int n);
Matrix init_matrix_non_contiguous(int m, int n);
void free_matrix_contiguous(Matrix *matrix);
void free_matrix_non_contiguous(Matrix *matrix);

int main(void) {
    Matrix matrix1 = init_matrix_contiguous(10, 8);
    Matrix matrix2 = init_matrix_non_contiguous(10, 8);
    
    // Test: Assign values to matrix1.data and matrix2.data

    free_matrix_contiguous(&matrix1);
    free_matrix_non_contiguous(&matrix2);

    return 0;
}

// Initialize a matrix with contiguous memory allocation
Matrix init_matrix_contiguous(int m, int n) {
    Matrix matrix;
    matrix.m = m;
    matrix.n = n;

    matrix.data = (double **) malloc(sizeof(double *) * m);
    matrix.data[0] = (double *) malloc(sizeof(double) * m * n);
    for (int i = 1; i < m; ++i) {
        matrix.data[i] = matrix.data[i-1] + n;
    }

    return matrix;
}

// Initialize a matrix with non-contiguous memory allocation
Matrix init_matrix_non_contiguous(int m, int n) {
    Matrix matrix;
    matrix.m = m;
    matrix.n = n;

    matrix.data = (double **) malloc(sizeof(double *) * m);
    for (int i = 0; i < m; ++i) {
        matrix.data[i] = (double *) malloc(sizeof(double) * n);
    }

    return matrix;
}

// Free allocated memory for a matrix which is stored contiguously in memory
void free_matrix_contiguous(Matrix *matrix) {
    if (matrix && matrix->data) {
        free(matrix->data[0]);
        free(matrix->data);
    }
}

// Free allocated memory for a matrix which is stored non-contiguously in memory
void free_matrix_non_contiguous(Matrix *matrix) {
    if (matrix && matrix->data) {
        for (int i = 0; i < matrix->m; ++i) {
            free(matrix->data[i]);
        }
        free(matrix->data);
    }
}

// Matrix-vector multiplication: b = A*x
void matvec(double * b, const Matrix * A, const double * x) {
    for(int i = 0; i < A->m; i++) {
        b[i] = 0; // initialize to zero
        for(int j = 0; j < A->n; j++) {
            b[i] += A->data[i][j] * x[j];
        }
    }
}

// Transpose matrix-vector multiplication: b = transpose(A)*x
void matvec_transpose(double * b, const Matrix * A, const double * x) {
    for(int j = 0; j < A->n; j++) {
        b[j] = 0; // initialize to zero
        for(int i = 0; i < A->m; i++) {
            b[j] += A->data[i][j] * x[i];
        }
    }
}
