#include <stdbool.h>
// axpy calculations with scalar y
void axpy_scalar_y(double solution[], double a, double x[], double y, int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        solution[i] = a * x[i] + ((y == 0) ? 0 : y);
    }
}

// axpy calculations with vector y
void axpy_vector_y(double solution[], double a, double x[], double y[], int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        solution[i] = a * x[i] + y[i];
    }
}

// overloaded axpy function
void axpy(double solution[], double a, double x[], double y[], int N) {
    if (y == NULL) { // If y is NULL, treat it as a scalar
        axpy_scalar_y(solution, a, x, 0, N);
    } else { // Otherwise, treat y as a vector
        axpy_vector_y(solution, a, x, y, N);
    }
}
// axpy calculations with vector y
void axpby(double solution[], double a, double x[],double b, double y[], int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        solution[i] = a * x[i] + b*y[i];
    }
}