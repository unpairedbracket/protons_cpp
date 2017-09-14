#pragma once

#include <cassert>

#include "math.h"
#include "structs.h"

struct ParticleSource {
    virtual ParticleState* genParticleState(ParticleInfo* particleType) = 0;
};

struct SquareSource : ParticleSource {
    double distance;
    double divergence;
    double energy;

    long x_extent, y_extent;

    ParticleState* genParticleState(ParticleInfo* particleType) override;
};

struct HelixSource : ParticleSource {
    double distance;
    double divergence;
    double energy;
    double dphi;

    long N;

    ParticleState* genParticleState(ParticleInfo* particleType) override;
};

void initSource(SquareSource* source, double distance, double divergence, double energy, long x_extent, long y_extent);

void initSource(HelixSource* source, double distance, double divergence, double energy, double pitch, long N);

void initPosPointSource(ParticleState* particles, double distance);
