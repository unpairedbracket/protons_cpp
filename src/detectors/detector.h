#pragma once

#include "../particles/structs.h"

struct ParticleDetector {
    ParticleInfo *particleInfo;
    double distance;

    void init(ParticleInfo* particle, double distance);
    void finalPush(ParticleState* state);
    virtual void detect(ParticleState* state) = 0;
    virtual void output() = 0;
};

struct DetectorTextFile : ParticleDetector {
    ParticleState* state;
    void detect(ParticleState* state) override;
    void output() override;
};

struct DetectorHDF5 : ParticleDetector {
    ParticleState* state;
    void detect(ParticleState* state) override;
    void output() override;
};

struct DetectorFluence : ParticleDetector {
    int detectorPixels[2];
    double detectorSize[2];
    double* detectorArray;

    void detect(ParticleState* state) override;
    void output() override;
};
