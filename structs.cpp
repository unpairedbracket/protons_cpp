#include "structs.h"

void initParticle(ParticleInfo* particle, double mass, double charge) {
    particle->mass = mass;
    particle->charge = charge;
    particle->qmratio = charge/mass;
}

void initParticleState(ParticleState* state, ParticleSource* source) {
    long N = source->x_extent * source->y_extent;
    state->particleInfo = source->particleInfo;
    
    state->pos = new Vector3[N];

    state->vel = new Vector3[N];

    state->running = new bool[N];

    memset(state->running, true, N*sizeof(bool));

    state->N = N;
    state->N_running = N;
}

void initSource(ParticleSource* source, ParticleInfo* particle, double distance, double divergence, double energy, long x_extent, long y_extent) {
    source->particleInfo = particle;
    source->distance = distance;
    source->divergence = divergence;
    source->energy = energy;

    source->x_extent = x_extent;
    source->y_extent = y_extent;
}

void initDetector(ParticleDetector* detector, ParticleInfo* particle, double distance) {
    detector->particleInfo = particle;
    detector->distance = distance;
}

