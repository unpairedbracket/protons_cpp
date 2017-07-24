#include "fields_cocoon.h"

void CocoonField::initFields() {
    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        this->E[j].x = 0.0;
        this->E[j].y = 0.0;
        this->E[j].z = 0.0;

        this->B[j].x = 0.0;
        this->B[j].y = 0.0;
        this->B[j].z = 0.0;
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

CocoonField* initCocoonField(double strength, double radius, double length, long N) {
    CocoonField* field = new CocoonField();
    field->N = N;
    field->E = new Vector3[N];

    field->B = new Vector3[N];

    field->B_strength = strength;
    field->r_scale = radius;
    field->z_scale = length;
    
    return field;
}
