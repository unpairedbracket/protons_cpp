#include "accel.h"

#include <cmath>

#include "physical_constants.h"

void accel(ParticleState* state, FieldStructure* fields) {
    fields->getFields(state);
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state->acc[j].x = state->particleInfo->qmratio * ( fields->E[j].x + state->vel[j].y * fields->B[j].z - state->vel[j].z * fields->B[j].y );
            state->acc[j].y = state->particleInfo->qmratio * ( fields->E[j].y + state->vel[j].z * fields->B[j].x - state->vel[j].x * fields->B[j].z );
            state->acc[j].z = state->particleInfo->qmratio * ( fields->E[j].z + state->vel[j].x * fields->B[j].y - state->vel[j].y * fields->B[j].x );
        }
    }
}

void accel_relativistic(ParticleState* state, FieldStructure* fields) {
    accel(state, fields);

    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            double beta_x = state->vel[j].x/c;
            double beta_y = state->vel[j].y/c;
            double beta_z = state->vel[j].z/c;

            double invgamma = sqrt(1 - (beta_x*beta_x + beta_y*beta_y + beta_z*beta_z));

            double a_dot_beta = state->acc[j].x*beta_x + state->acc[j].y*beta_y + state->acc[j].z*beta_z;
            state->acc[j].x = (state->acc[j].x - beta_x * a_dot_beta) * invgamma;
            state->acc[j].y = (state->acc[j].y - beta_y * a_dot_beta) * invgamma;
            state->acc[j].z = (state->acc[j].z - beta_z * a_dot_beta) * invgamma;
        }
    }
}

