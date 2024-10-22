#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dif_form_generators.h"
#include "sum_form_generators.h"

void test_difference_generators(){
    float* arr = malloc(sizeof(float)*9);
    gen_amc_dif(arr, 8); 
    for (size_t a = 0; a < 9; a++) printf("%f ", arr[a]);
    printf("\n");
    gen_abp_dif(arr, 8, NULL); 
    for (size_t a = 0; a < 9; a++) printf("%f ", arr[a]);
    printf("\n");
    gen_cc_dif(arr, 8, NULL); 
    for (size_t a = 0; a < 9; a++) printf("%f ", arr[a]);
    printf("\n");
    gen_sp_dif(arr, 8, NULL, NULL); 
    for (size_t a = 0; a < 9; a++) printf("%f ", arr[a]);
    free(arr);
}
void test_summed_generators(){
    float* arr = malloc(sizeof(float)*9);
    gen_amc_sum(arr, 8); 
    for (size_t a = 0; a < 8; a++) printf("%f ", arr[a]);
    printf("\n");
    gen_abp_sum(arr, 8, NULL); 
    for (size_t a = 0; a < 8; a++) printf("%f ", arr[a]);
    printf("\n");
    gen_cc_sum(arr, 8, NULL); 
    for (size_t a = 0; a < 8; a++) printf("%f ", arr[a]);
    printf("\n");
    gen_sp_sum(arr, 8, NULL, NULL); 
    for (size_t a = 0; a < 8; a++) printf("%f ", arr[a]);
    free(arr);
}

//Out size: sizeof(float)*(degree+2)*(degree+1)
//amc <- gen_amc_sum
//abp <- gen_abp_sum
void gen_gjc_dif(float* out, size_t degree, float* amc, float* abp){
    size_t width = degree+1;
    size_t height = degree+2;
    if (!out || !degree) return;
    if (!amc){
        amc = malloc(sizeof(float)*(degree+1));
        gen_amc_sum(amc, degree+1);
    }
    memcpy(out+width*(height-2), amc, sizeof(float)*(degree+1));
    if (!abp) gen_abp_sum(out+width*(height-1), degree+1, amc);
    else memcpy(out+width*(height-1), abp, sizeof(float)*(degree+1));
    for (size_t a = height-3; a <= degree; a--){
        for (size_t b = 0; b <= degree; b++){
            if (b == 0) out[a*width] = out[(a+1)*width];
            else out[a*width+b] = out[(a+1)*width+b] - out[(a+1)*width+b-1];
        }
    }
}

//Out size: sizeof(float)*(degree+2)*(degree+1)
//cc <- gen_cc_sum
//sp <- gen_sp_sum
//amc <- gen_amc_sum
void gen_gjp_dif(float* out, size_t degree, float* cc, float* sp, float* amc){
    size_t width = degree+1;
    size_t height = degree+2;
    if (!out || !degree) return;
    if (!cc){
        cc = malloc(sizeof(float)*(degree+1));
        gen_cc_sum(cc, degree+1, amc);
    }
    memset(out, 0, sizeof(float)*width*height);
    memcpy(out+width*(height-2), cc, sizeof(float)*(degree+1));
    //if (!abp) gen_abp_sum(out+width*(height-1), degree+1, amc);
    //else memcpy(out+width*(height-1), abp, sizeof(float)*(degree+1));
    //for (size_t a = height-3; a <= degree; a--){
    //    for (size_t b = 0; b <= degree; b++){
    //        if (b == 0) out[a*width] = out[(a+1)*width];
    //        else out[a*width+b] = out[(a+1)*width+b] - out[(a+1)*width+b-1];
    //    }
    //}
}

int main(int argc, char** argv){
    float* arr = malloc(sizeof(float)*10*9);
    gen_gjp_dif(arr, 8, NULL, NULL, NULL);
    for (size_t a = 0; a < 10; a++){
        for (size_t b = 0; b < 9; b++){
            printf("%f ", arr[a*9+b]);
        }
        printf("\n");
    }
    /*
    gen_gjc_dif(arr, 8, NULL, NULL);
    for (size_t a = 0; a < 10; a++){
        for (size_t b = 0; b < 9; b++){
            printf("%f ", arr[a*9+b]);
        }
        printf("\n");
    }
    */
}
