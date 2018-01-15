#pragma once

#include "../particles/structs.h"
#include "../sources/source.h"

struct Interpolator {
    int iterations;

    ParticleSource* interpSource;
    ParticleState* interpParticles;

    void initState(ParticleInfo type);

    virtual void setSamplePoints(ParticleState* samplePoints) = 0;
    virtual void setSampleValues(ParticleState* sampleValues) = 0;

    virtual void interpolate() = 0;
};

