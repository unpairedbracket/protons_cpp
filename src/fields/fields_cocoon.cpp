#include "fields_cocoon.h"

void CocoonField::initFields() {
    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        this->E[j] = { 0.0, 0.0, 0.0 };
        this->B[j] = { 0.0, 0.0, 0.0 };
    }
}

void CocoonField::getFields(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        if(state->running[j])
        {
            double Rscaled = sqrt(state->pos[j].x*state->pos[j].x + state->pos[j].y*state->pos[j].y) / this->r_scale;
            double Zscaled = (state->pos[j].z - 10 * this->z_scale/2)/this->z_scale;
            double BphiOnR = this->B_strength * exp(-Rscaled*Rscaled - Zscaled*Zscaled) / this->r_scale;

            this->B[j].x = - BphiOnR * state->pos[j].y;
            this->B[j].y = + BphiOnR * state->pos[j].x;
        }
    }
}

