
/*!iterative ilu for laplace matrix*/
void ilu(double **l,double u[], int N,double epsilon){
    int iteration_count = 0;
    while (iteration_count < 10){
        iteration_count += 1;
        
        #pragma omp parallel for
        for (int i = 0;i<N*N;i++){
            if (i-N >= 0){
                l[i][0]=-1/u[i-N];
            }
            // we have to skip some li's on the second diagonal
            if (i%N>0) {
                l[i][1]=-1/u[i-1];
            }
            u[i]=4+l[i][0]+l[i][1];
        }

    }
}