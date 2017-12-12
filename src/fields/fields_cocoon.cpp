#include "fields_cocoon.h"

using std::abs;

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
            double RscaledSquared = state->pos[j].x*state->pos[j].x + state->pos[j].y*state->pos[j].y / (this->r_scale*this->r_scale);
            double Zscaled = state->pos[j].z/this->z_scale;
            double BphiOnR = this->B_strength * exp(-RscaledSquared - Zscaled*Zscaled) / this->r_scale;

            this->B[j].x = - BphiOnR * state->pos[j].y;
            this->B[j].y = + BphiOnR * state->pos[j].x;
        }
    }
}

void CocoonField::invalidatePositions(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            if((abs(state->pos[j].z) > this->leeway * this->z_scale && state->vel[j].z * state->pos[j].z > 0)
            || (state->pos[j].x*state->pos[j].x + state->pos[j].y*state->pos[j].y > this->leeway*this->r_scale*this->leeway*this->r_scale && state->vel[j].x * state->pos[j].x + state->vel[j].y * state->pos[j].y > 0)
            ) {
                state->running[j] = false;
                #pragma omp atomic
                state->N_running--;
            }
        }
    }
}

