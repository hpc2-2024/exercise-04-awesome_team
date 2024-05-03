#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "axpby.h"

/*! The function f of the exercise sheet*/
double fun(double x, double y){
    return sin(y*M_PI)*sin(x*M_PI)*2.0*M_PI*M_PI;
}

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

/*! Calculating the dot (scalar) product of 2 vectors */
double dot(double v[], double w[], int size) {
    double sum = 0;
    int i;

    #pragma omp parallel for private(i) reduction(+:sum) 
    for (i=0;i<size;i++){
        sum += v[i]*w[i];
    } 

    return sum;
}

/*! Implementation of a matrix free multiplication with 5-star stencil*/
void mfMult(int N, double r[], double y[], double h){
    for (int i=1;i<N+1;i++){
        for (int j=1;j<N+1;j++){
            y[(N+2)*i+j]=4*r[(N+2)*i+j]-r[(N+2)*(i+1)+j] -r[(N+2)*i+j-1]-r[(N+2)*(i-1)+j]-r[(N+2)*i+j+1];
            y[(N+2)*i+j]=1.0/(h*h)*y[(N+2)*i+j];
        }
    }
}

/*! Displaying a vector */
void vec_print(int N, double vec[], char name[]){
    printf("\n %s \n",name);
    for (int i=0;i<N+2;i++){
        for (int j=0;j<N+2;j++){
            printf("%f ",vec[(N+2)*i+j]);
        }
        printf("\n");
    }
    printf("\n");
}



int main(int argc, char** argv){
    int preconditioner=0;
    int number_of_iterations = 0;
    // if you want to use a precondtioner first argument 1
    if (argc>1){
        preconditioner=atoi(argv[1]);
    }
    printf("%d \n",preconditioner);
    // initilize variables
    int N = 6; // N^2 is the number of inner points in our lattice
    int N2 = (N+2)*(N+2);
    double epsilon = 0.000001;
    double h = 1.0/(N+1);
    // init x0
    double alpha = 0;
    double beta=0;
    double err0;
    double *x,*p,*r,*b,*m, *z;
    x=(double *)malloc((N+2)*(N+2)*sizeof(double));
    z=(double *)malloc((N+2)*(N+2)*sizeof(double));
    p=(double *)malloc((N+2)*(N+2)*sizeof(double));
    r=(double *)malloc((N+2)*(N+2)*sizeof(double));
    b=(double *)malloc((N+2)*(N+2)*sizeof(double));
    m=(double *)malloc((N+2)*(N+2)*sizeof(double));
   
    // fill ghost layer with zeros (and everything else also 0)
    for (int i = 0;i<N+2;i++) {
        for (int j = 0;j<N+2;j++){
            x[(N+2)*i+j]=0;
            p[(N+2)*i+j]=0;
            r[(N+2)*i+j]=0;
            b[(N+2)*i+j]=0;
            m[(N+2)*i+j]=0;
        }

    }
    // inner points of x_0,b
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            x[(N+2)*i+j]=0;
            b[(N+2)*i+j]=fun(i*h,j*h);
        }

    }

    // Pre loop calculations ( calculating residuum )
    mfMult(N,x,r,h);
    axpy(r,-1,b,r,(N+2)*(N+2));
    axpy(p,-1,r,0,(N+2)*(N+2));
    #pragma omp parallel for
    for (int i=0;i<(N+2)*(N+2);i++){
        p[i]=r[i]*(-1);
    }
    err0 = sqrt(dot(r,r,N2));
    double old_r_dot, new_r_dot;
    old_r_dot = err0;

    do
    {
        number_of_iterations+=1;

        mfMult(N,p,m,h);
        alpha=old_r_dot/dot(p,m,N2);

        // update x,r
        #pragma omp parallel for
        for (int i=0;i<(N+2)*(N+2);i++){
            x[i]+=alpha*p[i];
            r[i]+=alpha*m[i];
        }
        //update p
        new_r_dot = dot(r,r,N2);
        beta = new_r_dot/old_r_dot;
        #pragma omp parallel for
        for (int i=0;i<(N+2)*(N+2);i++){
            p[i]=-r[i]+beta*p[i];
        }
        old_r_dot = new_r_dot;

    } while (sqrt(new_r_dot)/err0 >= epsilon);

    vec_print(N,x,"vector x");
    printf("Number of iterations: %d\n",number_of_iterations);


    // free allocated memory
    free(x);
    free(p);
    free(r);
    free(b);
    free(m);
    free(z);
    return 0;
}