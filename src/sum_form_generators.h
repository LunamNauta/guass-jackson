#ifndef SUM_FORM_GENERATORS
#define SUM_FORM_GENERATORS

#include <stdlib.h>

//Out size: sizeof(float)*degree
void gen_amc_sum(float* out, size_t degree){
    if (!out || !degree) return;
    for (size_t a = 1; a < degree+1; a++){
        out[a-1] = 0.0f;
        for (size_t b = 0; b < a; b++){
            float v = b == 0 ? 1.0f : out[b-1];
            out[a-1] -= v/(a-b+1);
        }
    }
}

//Out size: sizeof(float)*degree
//amc <- gen_amc_sum
void gen_abp_sum(float* out, size_t degree, float* amc){
    if (!out || !degree) return;
    int gen_amc = 0;
    if (!amc){
        amc = malloc(sizeof(float)*degree);
        gen_amc_sum(amc, degree);
        gen_amc = 1;
    }
    for (size_t a = 1; a < degree+1; a++){
        out[a-1] = 0.0f;
        for (size_t b = 0; b <= a; b++){
            float v = b == 0 ? 1.0f : amc[b-1];
            out[a-1] += v;
        }
    }
    if (gen_amc) free(amc);
}

//Out size: sizeof(float)*degree
//amc <- gen_amc_sum
void gen_cc_sum(float* out, size_t degree, float* amc){
    if (!out || !degree) return;
    int gen_amc = 0;
    if (!amc){
        amc = malloc(sizeof(float)*degree);
        gen_amc_sum(amc, degree);
        gen_amc = 1;
    }
    for (size_t a = 2; a < degree+1; a++){
        out[a-2] = 0.0f;
        for (size_t b = 0; b <= a; b++){
            float v1 = b == 0 ? 1.0f : (b == 1 ? -1.0f : amc[b-2]);
            float v2 = a == b ? 1.0f : ((a - b == 1) ? -1.0f : amc[a-b-2]);
            out[a-2] += v1*v2;
        }
    }
    if (gen_amc) free(amc);
}

//Out size: sizeof(float)*degree
//cc <- gen_cc_sum
//amc <- gen_amc_sum
void gen_sp_sum(float* out, size_t degree, float* cc, float* amc){
    if (!out || !degree) return;
    int gen_cc = 0;
    if (!cc){
        cc = malloc(sizeof(float)*degree);
        gen_cc_sum(cc, degree, amc);
        gen_cc = 1;
    }
    for (size_t a = 1; a < degree+1; a++){
        out[a-1] = 0.0f;
        for (size_t b = 0; b <= a; b++){
            float v = b == 0 ? 1.0f : cc[b-1];
            out[a-1] += v;
        }
    }
    if (gen_cc) free(cc);
}

#endif
