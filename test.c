#include <stdio.h>
#include <stdbool.h>
#include "axpby.h"

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

int main(int argc, char **argv){
    test_axpby();    
}