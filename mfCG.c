#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include "utils.h"
#include "axpby.h"
#include "ilu.h"

/*! The function f of the exercise sheet*/
double fun(double x, double y){
    return sin(y*M_PI)*sin(x*M_PI)*2.0*M_PI*M_PI;
}

double fun_solution(double x, double y){
    return sin(x*M_PI)*sin(y*M_PI);
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

void apply_precon(double a[][5],double r[],double temp[],double z[],int N){
    forward_solve(a,temp,r,N);
    backward_solve(a,z,temp,N);
}

int main(int argc, char** argv){
    int preconditioner=0;
    int number_of_iterations = 0;
    // if you want to use a precondtioner first argument 1
    if (argc>1){
        preconditioner=atoi(argv[1]);
    }
    printf("Preconditioner: %d \n",preconditioner);
    // initilize variables
    int N = 6; // N^2 is the number of inner points in our lattice
    int N2 = (N+2)*(N+2);
    double epsilon = 0.000001;
    double h = 1.0/(N+1);
    // init x0
    double alpha = 0;
    double beta=0;
    double err0, errk;
    double *x,*p,*r,*b,*m, *z,*temp;
    x=(double *)malloc((N+2)*(N+2)*sizeof(double));
    z=(double *)malloc((N+2)*(N+2)*sizeof(double));
    p=(double *)malloc((N+2)*(N+2)*sizeof(double));
    r=(double *)malloc((N+2)*(N+2)*sizeof(double));
    b=(double *)malloc((N+2)*(N+2)*sizeof(double));
    m=(double *)malloc((N+2)*(N+2)*sizeof(double));
    temp=(double *)malloc((N+2)*(N+2)*sizeof(double));

    // Creation of matrix a
    double a[5][N*N];
    lapl_matrix(a,N);
    ilu(a,N,0.001,10000);

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

    int seed = 123456;
    // inner points of x_0,b
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            // randomly initialize x with values in (0,1)
            srand(seed + i);
            double r = (double)rand() / (double)RAND_MAX;
            x[(N+2)*i+j]=r;

            b[(N+2)*i+j]=fun(i*h,j*h);
        }

    }

    // Pre loop calculations ( calculating residuum )
    mfMult(N,x,r,h); // r = Ax
    axpy(r,-1,b,r,(N+2)*(N+2)); // r = Ax -b (together with last line)
    if (preconditioner == 1) {
        //precondition z = M^-1r
        apply_precon(a,r,temp,z,N);
    }
    else {
        z = r; // no preconditioning of the residuum
    }
    axpy(p,-1,z,0,(N+2)*(N+2)); // p = -r (first conjugated gradient direction)

    err0 = sqrt(dot(z,z,N2)); // for break condition

    double old_r_dot, new_r_dot;
    old_r_dot = dot(r,z,N2); // dot product for loop

    do
    {
        number_of_iterations+=1;

        mfMult(N,p,m,h); // Ap
        alpha=old_r_dot/dot(p,m,N2); // rz/pAp

        // update x,r
        axpy(x,alpha,p,x,N2); //x = x + alpha*p
        axpy(r,alpha,m,r,N2); // r=r + alpha*Ap

        if (preconditioner==1) {
            // precondition r: z = Mr
            apply_precon(a,r,temp,z,N);
        }

        //update p
        new_r_dot = dot(r,r,N2);
        beta = new_r_dot/old_r_dot;
        axpby(p,-1,r,beta,p,(N+2)*(N+2));

        old_r_dot = new_r_dot;

        //break criteria
        if (preconditioner==0){
            errk=sqrt(new_r_dot);
        }
        else {
            errk=sqrt(dot(z,z,N2));
        }

    } while (errk/err0 >= epsilon);

    vec_print(N,x,"vector x");
    printf("Number of iterations: %d\n",number_of_iterations);

    // Compute the relative absolute difference 
    double abs_diff = 0.0;
    double abs_sol = 0.0;
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            double x_sol = fun_solution(i*h,j*h);
            abs_diff += fabs(x[(N+2)*i+j] - x_sol);
            abs_sol += x_sol;
        }
    }
    printf("Relative absolute difference between exact and approx. solution: %f%%\n", 100.0 * abs_diff / abs_sol);

    // free allocated memory
    free(x);
    free(p);
    free(r);
    free(b);
    free(m);
    if (preconditioner==1){
        free(z);
    }
    return 0;
}