#pragma once

struct ParticleInfo {
    double mass;
    double charge;

    double qmratio;
};

struct ParticleSource {
    ParticleInfo *particleInfo;
    double distance;
    double divergence;
    double energy;
};

struct ParticleDetector {
    ParticleInfo *particleInfo;
    double distance;
};

ParticleInfo* makeParticle(double mass, double charge);
ParticleSource* makeSource(ParticleInfo* particle, double distance, double divergence, double energy);
ParticleDetector* makeDetector(ParticleInfo* particle, double distance);
