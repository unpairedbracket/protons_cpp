#include "fields_cocoon.h"

void CocoonField::initFields() {
    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        this->Ex[j] = 0.0;
        this->Ey[j] = 0.0;
        this->Ez[j] = 0.0;

        this->Bx[j] = 0.0;
        this->By[j] = 0.0;
        this->Bz[j] = 0.0;
    }
}

void CocoonField::getFields(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        if(state->running[j])
        {
            double Rscaled = sqrt(state->posX[j]*state->posX[j] + state->posY[j]*state->posY[j]) / this->r_scale;
            double Zscaled = (state->posZ[j] - 10 * this->z_scale/2)/this->z_scale;
            double BphiOnR = this->B_strength * exp(-Rscaled*Rscaled - Zscaled*Zscaled) / this->r_scale;

            this->Bx[j] = - BphiOnR * state->posY[j];
            this->By[j] = + BphiOnR * state->posX[j];
        }
    }
}

CocoonField* initCocoonField(double strength, double radius, double length, long N) {
    CocoonField* field = new CocoonField();
    field->N = N;
    field->Ex = new double[N];
    field->Ey = new double[N];
    field->Ez = new double[N];

    field->Bx = new double[N];
    field->By = new double[N];
    field->Bz = new double[N];

    field->B_strength = strength;
    field->r_scale = radius;
    field->z_scale = length;
    
    return field;
}
