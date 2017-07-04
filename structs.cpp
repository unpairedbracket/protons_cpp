#include "structs.h"

ParticleInfo* makeParticle(double mass, double charge) {
    ParticleInfo* particle = new ParticleInfo();
    particle->mass = mass;
    particle->charge = charge;
    particle->qmratio = charge/mass;
    return particle;
}

ParticleState* makeParticleState(ParticleSource* source) {
    ParticleState* state = new ParticleState();

    long N = source->x_extent * source->y_extent;
    state->particleInfo = source->particleInfo;
    
    state->posX = new double[N];
    state->posY = new double[N];
    state->posZ = new double[N];

    state->velX = new double[N];
    state->velY = new double[N];
    state->velZ = new double[N];

    state->N = N;
    return state;
}

ParticleSource* makeSource(ParticleInfo* particle, double distance, double divergence, double energy, long x_extent, long y_extent) {
    ParticleSource* source = new ParticleSource();
    source->particleInfo = particle;
    source->distance = distance;
    source->divergence = divergence;
    source->energy = energy;

    source->x_extent = x_extent;
    source->y_extent = y_extent;
    return source;
}

ParticleDetector* makeDetector(ParticleInfo* particle, double distance) {
    ParticleDetector* detector = new ParticleDetector();
    detector->particleInfo = particle;
    detector->distance = distance;
    return detector;
}
