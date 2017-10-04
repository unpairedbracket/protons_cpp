#pragma once

#include "../particles/structs.h"

struct ParticleDetector {
    ParticleInfo *particleInfo;
    double distance;

    void init(ParticleInfo* particle, double distance);
    void finalPush(ParticleState* state);
    virtual void output(ParticleState* state) = 0;
};

struct DetectorNoop : ParticleDetector {
    void output(ParticleState* state) override {};
};
