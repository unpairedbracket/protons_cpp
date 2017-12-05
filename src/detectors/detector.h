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

struct DetectorTextFile : ParticleDetector {
    void output(ParticleState* state) override;
};

struct DetectorHDF5 : ParticleDetector {
    void output(ParticleState* state) override;
};
