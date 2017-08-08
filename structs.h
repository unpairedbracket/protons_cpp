#pragma once

#include <cstring>
#include <cstdio>

struct Vector3 {
    double x, y, z;
};

struct ParticleInfo {
    double mass;
    double charge;

    double qmratio;
};

struct ParticleState {
    ParticleInfo *particleInfo;
    Vector3 *pos;
    Vector3 *vel;
    Vector3 *acc;
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

void initParticle(ParticleInfo* particle);
void initParticleState(ParticleState* state, ParticleSource* source);
void shadowParticleState(ParticleState* state, ParticleState* other);
void initSource(ParticleSource* source, ParticleInfo* particle, double distance, double divergence, double energy, long x_extent, long y_extent);
void initDetector(ParticleDetector* detector, ParticleInfo* particle, double distance);
