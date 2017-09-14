#include "structs.h"

void initParticle(ParticleInfo* particle) {
    if(particle->mass == 0)
    {
        printf("Attempt to use massless particle. This will not work well.\n");
        printf("Make sure your particle is intialised\n");
        particle->qmratio = 0;
        return;
    }

    if(particle->charge == 0) {
        printf("Attempt to use neutral particle. Nothing will happen!");
    }

    particle->qmratio = particle->charge/particle->mass;
}

void initParticleState(ParticleState* state, ParticleInfo* particleType, long N) {
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

void initDetector(ParticleDetector* detector, ParticleInfo* particle, double distance) {
    detector->particleInfo = particle;
    detector->distance = distance;
}

