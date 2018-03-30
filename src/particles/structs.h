#pragma once

struct Vector3 {
    double x, y, z;
};

struct ParticleInfo {
    double mass;
    double charge;

    double qmratio;
};

struct ParticleState {
    ParticleInfo particleInfo;
    Vector3 *pos;
    Vector3 *vel;
    Vector3 *acc;

    Vector3 *intE;
    Vector3 *intB;

    bool *running;
    long N;
    long N_running;
};

void initParticleState(ParticleState* state, ParticleInfo particleType, long N);
void setAllRunning(ParticleState* state);
void shadowParticleState(ParticleState* state, ParticleState* other);
