#include <stdio.h>
#include <stdbool.h>
#include "axpby.h"
#include "ilu.h"
#include "utils.h"

void test_axpby(){
    double x[] = {0,0,0};
    double e1[] = {1,0,0};
    double e2[] = {0,1,0};
    double e3[] = {0,0,1};
    double a = 2;

    //TC1
    bool fail = false;
    axpby(x,a,e1,3,e3,3);
    double expected_result[] = {2,0,3};
    for (int i=0;i<3;i++){
        if (expected_result[i]!=x[i]){
            fail = true;
        }
    }
    if (fail==true){
        printf("Testcase 1 has failed\n");
    }
    else {
        printf("Test case 1 has succeded\n");
    }

    //TC2
    fail = false;
    axpy(x,a,e1,0,3);
    double exp_res2[] = {2,0,0};
    for (int i=0;i<3;i++){
        if (exp_res2[i]!=x[i]){
            fail = true;
        }
    }
    if (fail==true){
        printf("Testcase 2 has failed\n");
    }
    else {
        printf("Test case 2 has succeded\n");
    }



    //TC3

}

void test_ilu(){
    //
    int N = 4;
    double l[16][2] = {{0,0},{0,-1},{0,-1},{0,-1},{-1,0},{-1,-1},{-1,-1},{-1,-1},{-1,0},{-1,-1},{-1,-1},{-1,-1},{-1,0},{-1,-1},{-1,-1},{-1,-1}};
    double u[16] = {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};

    ilu(l,u,N,0,2);

    print_1dim(N*N,u,"u diagonal:");
}

void test_lapl(){
    int N=3;
    double a[9][5];

    lapl_matrix(a,N);

    double exp_matrix[9][5]={{0,0,4,-1,-1},{0,-1,4,-1,-1},{0,-1,4,0,-1},{-1,0,4,-1,-1},
        {-1,-1,4,-1,-1},{-1,-1,4,0,-1},{-1,0,4,-1,-0},{-1,-1,4,-1,0},{-1,-1,4,0,0}};
    
    bool fail = false;
    for (int i=0;i<5;i++){
        for (int j=0;j<N*N;j++){
            if (a[j][i]!=exp_matrix[j][i]){
                fail=true;
            }
        }
    }
    if (fail==true){
        printf("Test: Creation of lapl matrix has failed.\n");
    }
    else {
        printf("Test: Creation of lapl matrix has succeded.\n");
    }
    
}

int main(int argc, char **argv){
    //test_axpby();    
    //test_ilu();
    test_lapl();
}