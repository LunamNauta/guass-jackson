#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void gen_amc_dif(float* out, size_t degree){
    if (!degree || !out) return;
    out[0] = 1.0f;
    for (size_t a = 1; a < degree; a++){
        out[a] = 0.0f;
        for (size_t b = 0; b < a; b++){
            out[a] -= out[b]/(a-b+1);
        }
    }
}
void gen_abp_dif(float* out, size_t degree, float* amc){
    if (!degree || !out) return;
    int amc_gen = 0;
    if (!amc){
        amc = malloc(sizeof(float)*degree);
        if (!amc) return;
        gen_amc_dif(amc, degree);
        amc_gen = 1;
    }
    for (size_t a = 0; a < degree; a++){
        out[a] = 0.0f;
        for (size_t b = 0; b <= a; b++){
            out[a] += amc[b];
        }
    }
    if (amc_gen) free(amc);
}

void gen_cc_dif(float* out, size_t degree, float* amc){
    if (!degree || !out) return;
    int amc_gen = 0;
    if (!amc){
        amc = malloc(sizeof(float)*degree);
        if (!amc) return;
        gen_amc_dif(amc, degree);
        amc_gen = 1;
    }
    for (size_t a = 0; a < degree; a++){
        out[a] = 0.0f;
        for (size_t b = 0; b <= a; b++){
            out[a] += amc[b]*amc[a-b];
        }
    }
    if (amc_gen) free(amc);
}
void gen_sp_dif(float* out, size_t degree, float* cc, float* amc){
    if (!degree || !out) return;
    int cc_gen = 0;
    if (!cc){
        cc = malloc(sizeof(float)*degree);
        if (!cc) return;
        gen_cc_dif(cc, degree, amc);
        cc_gen = 1;
    }
    for (size_t a = 0; a < degree; a++){
        out[a] = 0.0f;
        for (size_t b = 0; b <= a; b++){
            out[a] += cc[b];
        }
    }
    if (cc_gen) free(cc);
}

void gen_amc_sum(float* out, size_t degree){
    if (!degree || !out) return;
    for (size_t a = 1; a < degree; a++){
        out[a-1] = 0.0f;
        for (size_t b = 0; b < a; b++){
            float v = b == 0 ? 1.0f : out[b-1];
            out[a-1] -= v/(a-b+1);
        }
    }
}
void gen_abp_sum(float* out, size_t degree, float* amc){
    if (!degree || !out) return;
    int amc_gen = 0;
    if (!amc){
        amc = malloc(sizeof(float)*(degree-1));
        if (!amc) return;
        gen_amc_sum(amc, degree);
        amc_gen = 1;
    }
    for (size_t a = 1; a < degree; a++){
        out[a-1] = 0.0f;
        for (size_t b = 0; b <= a; b++){
            float v = b == 0 ? 1.0f : amc[b-1];
            out[a-1] += v;
        }
    }
    if (amc_gen) free(amc);
}
void gen_cc_sum(float* out, size_t degree, float* amc){
    if (!degree || !out) return;
    int amc_gen = 0;
    if (!amc){
        amc = malloc(sizeof(float)*(degree-1));
        if (!amc) return;
        gen_amc_sum(amc, degree);
        amc_gen = 1;
    }
    for (size_t a = 1; a < degree; a++){
        out[a-1] = 0.0f;
        for (size_t b = 0; b <= a; b++){
            float v1 = b == 0 ? 1.0f : amc[b-1];
            float v2 = a == b ? 1.0f : amc[a-b-1];
            out[a-1] += v1*v2;
        }
    }
    if (amc_gen) free(amc);
}
void gen_sp_sum(float* out, size_t degree, float* cc, float* amc){
    if (!degree || !out) return;
    int cc_gen = 0;
    if (!cc){
        cc = malloc(sizeof(float)*(degree-1));
        if (!cc) return;
        gen_cc_sum(cc, degree, amc);
        cc_gen = 1;
    }
    for (size_t a = 1; a < degree; a++){
        out[a-1] = 0.0f;
        for (size_t b = 0; b <= a; b++){
            float v = b == 0 ? 1.0f : cc[b-1];
            out[a-1] += v;
        }
    }
    if (cc_gen) free(cc);
}

void gen_gjc(float* out, size_t degree, float* abp, float* amc){
    if (!degree || !out) return;
    memset(out, 0, sizeof(float)*(degree+2)*(degree+1)); //TODO: Remove
    int amc_gen = 0;
    if (!amc){
        amc = malloc(sizeof(float)*(degree+1));
        if (!amc) return;
        gen_amc_sum(amc, degree+2);
        amc_gen = 1;
    }
    memcpy(out+degree*(degree+1), amc, sizeof(float)*(degree+1));
    if (!abp) gen_abp_sum(out+(degree+1)*(degree+1), degree+2, amc);
    else memcpy(out+(degree+1)*(degree+1), abp, sizeof(float)*(degree+1));
    for (size_t a = 0; a < degree+2; a++){
        for (size_t b = 0; b < degree; b++){
            //out[a*(degree+2)+b]
        }
    }
    if (amc_gen) free(amc);
}

void Test_Difference_Coefficients(){
    printf("Adams-Moulton Corrector (Dif) Coefficients:\n");
    float* amc_dif = malloc(sizeof(float)*9);
    gen_amc_dif(amc_dif, 9);
    for (size_t a = 0; a < 9; a++) printf("%f ", amc_dif[a]);
    printf("\n\n");

    printf("Adams-Bashforth Predictor (Dif) Coefficients:\n");
    float* abp_dif = malloc(sizeof(float)*9);
    gen_abp_dif(abp_dif, 9, amc_dif);
    for (size_t a = 0; a < 9; a++) printf("%f ", abp_dif[a]);
    printf("\n\n");

    printf("Cowell Corrector (Dif) Coefficients:\n");
    float* cc_dif = malloc(sizeof(float)*9);
    gen_cc_dif(cc_dif, 9, amc_dif);
    for (size_t a = 0; a < 9; a++) printf("%f ", cc_dif[a]);
    printf("\n\n");

    printf("Stormer Predictor (Dif) Coefficients:\n");
    float* sp_dif = malloc(sizeof(float)*9);
    gen_sp_dif(sp_dif, 9, NULL, NULL);
    for (size_t a = 0; a < 9; a++) printf("%f ", sp_dif[a]);
    printf("\n\n");

    free(amc_dif);
    free(abp_dif);
    free(cc_dif);
    free(sp_dif);
}
void Test_Summed_Coefficients(){
    printf("Adams-Moulton Corrector (Sum) Coefficients:\n");
    float* amc_sum = malloc(sizeof(float)*8);
    gen_amc_sum(amc_sum, 9);
    for (size_t a = 0; a < 8; a++) printf("%f ", amc_sum[a]);
    printf("\n\n");

    printf("Adams-Bashforth Predictor (Sum) Coefficients:\n");
    float* abp_sum = malloc(sizeof(float)*8);
    gen_abp_sum(abp_sum, 9, amc_sum);
    for (size_t a = 0; a < 8; a++) printf("%f ", abp_sum[a]);
    printf("\n\n");

    printf("Cowell Corrector (Sum) Coefficients:\n");
    float* cc_sum = malloc(sizeof(float)*8);
    gen_cc_sum(cc_sum, 9, amc_sum);
    for (size_t a = 0; a < 8; a++) printf("%f ", cc_sum[a]);
    printf("\n\n");

    printf("Stormer Predictor (Sum) Coefficients:\n");
    float* sp_sum = malloc(sizeof(float)*8);
    gen_sp_sum(sp_sum, 9, cc_sum, amc_sum);
    for (size_t a = 0; a < 8; a++) printf("%f ", sp_sum[a]);
    printf("\n\n");

    free(amc_sum);
    free(abp_sum);
    free(cc_sum);
    free(sp_sum);
}

int main(int argc, char** argv){
    printf("Gauss Jackson Corrector Coefficients:\n");
    float* gjc = malloc(sizeof(float)*10*9);
    gen_gjc(gjc, 8, NULL, NULL);
    for (size_t a = 0; a < 10; a++){
        for (size_t b = 0; b < 9; b++){
            printf("%f ", gjc[a*9+b]);
        }
        printf("\n");
    }
    free(gjc);
}
