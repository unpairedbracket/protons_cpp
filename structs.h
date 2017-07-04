#pragma once

struct ParticleInfo {
    double mass;
    double charge;

    double qmratio;
};

struct ParticleState {
    ParticleInfo *particleInfo;
    double *posX, *posY, *posZ;
    double *velX, *velY, *velZ;
    long N;
};

struct ParticleSource {
    ParticleInfo *particleInfo;
    double distance;
    double divergence;
    double energy;

    long x_extent, y_extent;
};

struct ParticleDetector {
    ParticleInfo *particleInfo;
    double distance;
};

ParticleInfo* makeParticle(double mass, double charge);
ParticleState* makeParticleState(ParticleSource* source);
ParticleSource* makeSource(ParticleInfo* particle, double distance, double divergence, double energy, long x_extent, long y_extent);
ParticleDetector* makeDetector(ParticleInfo* particle, double distance);
