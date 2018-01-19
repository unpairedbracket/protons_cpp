#pragma once

#include "../particles/structs.h"

struct ParticleDetector {
    ParticleInfo *particleInfo;
    double distance;

    void init(ParticleInfo* particle, double distance);
    void finalPush(ParticleState* state);
    virtual void detectUndeviated(ParticleState* state) = 0;
    virtual void detect(ParticleState* state) = 0;
    virtual void output() = 0;
};

struct DetectorTextFile : ParticleDetector {
    ParticleState* state;
    void detectUndeviated(ParticleState* state) override;
    void detect(ParticleState* state) override;
    void output() override;
};

struct DetectorHDF5 : ParticleDetector {
    ParticleState* state;
    void detectUndeviated(ParticleState* state) override;
    void detect(ParticleState* state) override;
    void output() override;
};

struct DetectorFluence : ParticleDetector {
    int detectorPixels[2];
    double detectorSize[2];
    double* detectorArray;
    double* nullDetectorArray;

    void detectUndeviated(ParticleState* state) override;
    void detect(ParticleState* state) override;
    void output() override;
};
