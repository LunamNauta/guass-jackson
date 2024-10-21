#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float* gen_amc_dif_coe(float* c, size_t count){
    if (!count) return NULL;
    if (!c) c = malloc(sizeof(float)*count);
    if (!c) return NULL;
    c[0] = 1.0f;
    for (size_t a = 1; a < count; a++){
        c[a] = 0.0f;
        for (size_t b = 0; b < a; b++){
            c[a] -= c[b]/(a-b+1);
        }
    }
    return c;
}
float* gen_amc_sum_coe(float* c, size_t count){
    if (!count) return NULL;
    if (!c) c = malloc(sizeof(float)*(count-1));
    if (!c) return NULL;
    for (size_t a = 0; a < count; a++){
        c[a] = 0.0f;
        for (size_t b = 0; b < a; b++){
            float v = b == 0 ? 1.0f : c[b];
            c[a] -= v/(a-b+1);
        }
    }
    return c;
}
float* gen_abp_dif_coe(float* y, size_t count, float* c){
    if (!count) return NULL;
    if (!c) c = gen_amc_dif_coe(c, count);
    if (!c) return NULL;
    if (!y) y = malloc(sizeof(float)*count);
    if (!y) return NULL;
    for (size_t a = 0; a < count; a++){
        y[a] = 0.0f;
        for (size_t b = 0; b < a+1; b++){
            y[a] += c[b];
        }
    }
    return y;
}
float* gen_cc_dif_coe(float* q, size_t count, float* c){
    if (!count) return NULL;
    if (!c) c = gen_amc_dif_coe(c, count);
    if (!c) return NULL;
    if (!q) q = malloc(sizeof(float)*count);
    if (!q) return NULL;
    for (size_t a = 0; a < count; a++){
        q[a] = 0.0f;
        for (size_t b = 0; b < a+1; b++){
            q[a] += c[b]*c[a-b];
        }
    }
    return q;
}
float* gen_sp_dif_coe(float* l, size_t count, float* q, float* c){
    if (!count) return NULL;
    if (!q) q = gen_cc_dif_coe(q, count, c);
    if (!q) return NULL;
    if (!l) l = malloc(sizeof(float)*count);
    if (!l) return NULL;
    for (size_t a = 0; a < count; a++){
        l[a] = 0.0f;
        for (size_t b = 0; b < a+1; b++){
            l[a] += q[b];
        }
    }
    return l;
}
float* gen_gjp_coe(float* b, size_t count, float* y, float* c){
    if (!count) return NULL;
    if (!b) b = malloc(sizeof(float)*(count+2)*(count+1));
    if (!b) return NULL;
    if (!c) c = gen_amc_dif_coe(b+(count-1)*count, count+1);
    else memcpy(b+(count-1)*count, c, sizeof(float)*count+1);
    if (!c) return NULL;
    if (!y) y = gen_abp_dif_coe(b+count*count, count+1, c);
    else memcpy(b+count*count, y, sizeof(float)*count+1);
    if (!y) return NULL;
    return b;
}

int main(int argc, char** argv){
    float* gjp = gen_amc_sum_coe(NULL, 9);
    for (size_t a = 0; a < 8; a++){
        printf("%f\n", gjp[a]);
    }
    /*
    for (size_t a = 0; a < 10; a++){
        for (size_t b = 0; b < 10; b++) printf("%f ", gjp[a*9+b]);
        printf("\n");
    }
    */

    return 0;
}
