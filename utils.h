#include <stdio.h>
#include <math.h>

double norm_L1(double *x, double *y, int N){
    double norm = 0;
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            norm += fabs(x[i * N + j] - y[i * N + j]);
        }
    }

    return norm / (N*N);
}

double norm_L2(double *x, double *y, int N){
    double norm = 0;
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            norm += (x[i * N + j] - y[i * N + j]) * (x[i * N + j] - y[i * N + j]);
        }
    }

    return sqrt(norm) / (N*N);
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

/*! Displaying a vector with ghostlayer*/
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

void print_1dim(int N, double vec[],char name[]){
    printf("\n %s \n",name);
    for (int i=0;i<N;i++){
        printf("%f\n",vec[i]);
    }
}

void print_2dim(int N, double vec[][5],char name[]){
    printf("\n %s \n",name);
        for (int i=0;i<N*N;i++){
            for (int j=0;j<5;j++){
                printf("%f ",vec[i][j]);
        }
        printf("\n");
    }
}