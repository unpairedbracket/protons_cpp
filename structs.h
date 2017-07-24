#pragma once

#include <cstring>

struct ParticleInfo {
    double mass;
    double charge;

    double qmratio;
};

struct ParticleState {
    ParticleInfo *particleInfo;
    double *posX, *posY, *posZ;
    double *velX, *velY, *velZ;
    bool *running;
    long N;
    long N_running;
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

void initParticle(ParticleInfo* particle, double mass, double charge);
void initParticleState(ParticleState* state, ParticleSource* source);
void initSource(ParticleSource* source, ParticleInfo* particle, double distance, double divergence, double energy, long x_extent, long y_extent);
void initDetector(ParticleDetector* detector, ParticleInfo* particle, double distance);
