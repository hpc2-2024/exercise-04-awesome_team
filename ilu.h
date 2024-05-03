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

void forward_solve(double a[][5],double y[], double b[],int N){
    for (int i=0;i<N*N;i++){
        if (i==0){
            y[i]=b[i];
        }
        else {
            if (i-N>=0){
                y[i]=b[i]-y[i-N]*a[i][0]-a[i][1]*y[i-1];
            }
            else {
                y[i]=b[i]-y[i-1]*a[i][1];
            }
        }
    }
}

void backward_solve(double a[][5],double x[],double y[],int N){
    for (int i=N*N-1;i>=0;i--){
        if (i==N*N-1){
            x[i]=y[i]/a[i][2];
        }
        else {
            if (i+N<=N*N-1){
                x[i]=(y[i]-a[i][3]*x[i+1]-a[i][4]*x[i+N])/a[i][2];
            }
            else {
                x[i]=(y[i]-a[i][3]*x[i+1])/a[i][2];
            }
        }
    }

}

void precondition_ilu(double a[][5],double x[],double y[],double b[],int N){
    // forward solve
    forward_solve(a,y,b,N);
    
    //backward solve
    backward_solve(a,x,y,N);

}