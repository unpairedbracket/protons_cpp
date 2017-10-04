#include "detector.h"

void ParticleDetector::init(ParticleInfo* particle, double distance) {
    this->particleInfo = particle;
    this->distance = distance;
}

void ParticleDetector::finalPush(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        state->pos[j].x += (distance - state->pos[j].z) * state->vel[j].x / state->vel[j].z;
        state->pos[j].y += (distance - state->pos[j].z) * state->vel[j].y / state->vel[j].z;
        state->pos[j].z += (distance - state->pos[j].z);
    }
}

