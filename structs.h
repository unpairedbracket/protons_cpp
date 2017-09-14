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

struct ParticleDetector {
    ParticleInfo *particleInfo;
    double distance;
};

void initParticle(ParticleInfo* particle);
void initParticleState(ParticleState* state, ParticleInfo* particleType, long N);
void shadowParticleState(ParticleState* state, ParticleState* other);
void initDetector(ParticleDetector* detector, ParticleInfo* particle, double distance);
