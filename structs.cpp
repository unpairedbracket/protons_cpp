#include "structs.h"

ParticleInfo* makeParticle(double mass, double charge) {
    ParticleInfo* particle = new ParticleInfo();
    particle->mass = mass;
    particle->charge = charge;
    particle->qmratio = charge/mass;
    return particle;
}

ParticleSource* makeSource(ParticleInfo* particle, double distance, double divergence, double energy) {
    ParticleSource* source = new ParticleSource();
    source->particleInfo = particle;
    source->distance = distance;
    source->divergence = divergence;
    source->energy = energy;
    return source;
}

ParticleDetector* makeDetector(ParticleInfo* particle, double distance) {
    ParticleDetector* detector = new ParticleDetector();
    detector->particleInfo = particle;
    detector->distance = distance;
    return detector;
}
