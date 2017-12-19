#pragma once

#include "../particles/structs.h"
#include "../sources/source.h"

struct Interpolator {
    int iterations;

    ParticleSource* interpSource;
    ParticleState* interpParticles;

    virtual void setSamplePoints(ParticleState* samplePoints) = 0;
    virtual void setSampleValues(ParticleState* sampleValues) = 0;
    virtual void initState(ParticleInfo* type) = 0;
    virtual void interpolate() = 0;
};

