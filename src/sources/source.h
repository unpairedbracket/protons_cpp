#pragma once

#include "../particles/structs.h"

struct ParticleSource {
    double distance;
    double divergence;

    virtual ParticleState* genParticleState(ParticleInfo* particleType) = 0;
};

struct SquareSource : ParticleSource {
    double energy;

    long x_extent, y_extent;

    ParticleState* genParticleState(ParticleInfo* particleType) override;
};

struct HelixSource : ParticleSource {
    double energy;
    double dphi;

    long N;

    ParticleState* genParticleState(ParticleInfo* particleType) override;
};

void initSource(SquareSource* source, double distance, double divergence, double energy, long x_extent, long y_extent);

void initSource(HelixSource* source, double distance, double divergence, double energy, double pitch, long N);

void initPosPointSource(ParticleState* particles, double distance);
