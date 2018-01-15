#include "structs.h"

#include <cstring>
#include <cstdio>

void initParticleState(ParticleState* state, ParticleInfo particleType, long N) {
    state->particleInfo = particleType;
    
    state->pos = new Vector3[N];
    state->vel = new Vector3[N];
    state->acc = new Vector3[N];

    state->running = new bool[N];

    memset(state->running, true, N*sizeof(bool));

    state->N = N;
    state->N_running = N;
}

void shadowParticleState(ParticleState* state, ParticleState* other) {
    state->particleInfo = other->particleInfo;
    
    state->pos = new Vector3[other->N];
    state->vel = new Vector3[other->N];
    state->acc = new Vector3[other->N];

    state->running = other->running;

    state->N = other->N;
    state->N_running = 0;
}

