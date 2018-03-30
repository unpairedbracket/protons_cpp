#pragma once

#include "interpolator.h"

#include "../particles/structs.h"
#include "../sources/source.h"

struct BilinearInterpolator : Interpolator {
    int n_cells[2];
    // The size here is tangent of the (full) x/y opening angle
    double size[2];

    ParticleState sampleVelocities;
    ParticleState sampleDeflections;

    void setSamplePoints(ParticleState* samplePoints) override;
    void setSampleValues(ParticleState* sampleValues) override;

    void interpolate(ParticleState* probeState) override;
};

