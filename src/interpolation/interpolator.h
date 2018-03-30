#pragma once

#include "../particles/structs.h"
#include "../sources/source.h"

struct Interpolator {
    virtual void setSamplePoints(ParticleState* samplePoints) = 0;
    virtual void setSampleValues(ParticleState* sampleValues) = 0;

    virtual void interpolate(ParticleState* probeState) = 0;
};

