void lapl_matrix(double a[][5], int N){

    for (int i=0;i<N*N;i++) {
        //diagonal
        a[i][2]=4;

        if (i%N==0){
            a[i][1]=0;
        }
        else {
            a[i][1]=-1;
        }
        if (i%N==N-1){
            a[i][3]=0;
        }
        else {
            a[i][3]=-1;
        }
    }
    //outer -1 diagonals
    for (int i=N;i<N*N;i++){
        a[i][0]=-1;
    }
    for (int i=0;i<N*N-N;i++){
        a[i][4]=-1;
    }
    for (int i=0;i<N;i++){
        a[i][0]=0;
    }
    for (int i=N*N-N;i<N*N;i++){
        a[i][4]=0;
    }
}

/*!iterative ilu for laplace matrix*/
void ilu(double l[][2],double u[], int N,double epsilon,int max_it){
    int iteration_count = 0;
    while (iteration_count < max_it){
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