#ifndef DIF_FORM_GENERATORS
#define DIF_FORM_GENERATORS

#include <stdlib.h>

//Out size: sizeof(float)*(degree+1)
void gen_amc_dif(float* out, size_t degree){
    if (!out || !degree) return;
    out[0] = 1.0f;
    for (size_t a = 1; a < degree+1; a++){
        out[a] = 0.0f;
        for (size_t b = 0; b < a; b++){
            out[a] -= out[b]/(a-b+1);
        }
    }
}

//Out size: sizeof(float)*(degree+1)
//amc <- gen_amc_dif
void gen_abp_dif(float* out, size_t degree, float* amc){
    if (!out || !degree) return;
    int gen_amc = 0;
    if (!amc){
        amc = malloc(sizeof(float)*(degree+1));
        gen_amc_dif(amc, degree);
        gen_amc = 1;
    }
    for (size_t a = 0; a < degree+1; a++){
        out[a] = 0.0f;
        for (size_t b = 0; b <= a; b++){
            out[a] += amc[b];
        }
    }
    if (gen_amc) free(amc);
}

//Out size: sizeof(float)*(degree+1)
//amc <- gen_amc_dif
void gen_cc_dif(float* out, size_t degree, float* amc){
    if (!out || !degree) return;
    int gen_amc = 0;
    if (!amc){
        amc = malloc(sizeof(float)*(degree+1));
        gen_amc_dif(amc, degree);
        gen_amc = 1;
    }
    for (size_t a = 0; a < degree+1; a++){
        out[a] = 0.0f;
        for (size_t b = 0; b <= a; b++){
            out[a] += amc[b]*amc[a-b];
        }
    }
    if (gen_amc) free(amc);
}

//Out size: sizeof(float)*(degree+1)
//cc <- gen_cc_dif
//amc <- gen_amc_dif
void gen_sp_dif(float* out, size_t degree, float* cc, float* amc){
    if (!out || !degree) return;
    int gen_cc = 0;
    if (!cc){
        cc = malloc(sizeof(float)*(degree+1));
        gen_cc_dif(cc, degree, amc);
        gen_cc = 1;
    }
    for (size_t a = 0; a < degree+1; a++){
        out[a] = 0.0f;
        for (size_t b = 0; b <= a; b++){
            out[a] += cc[b];
        }
    }
    if (gen_cc) free(cc);
}

#endif
