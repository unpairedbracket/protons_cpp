#pragma once

#include <string>

#include "fields.h"

struct FlashField : FieldStructure {
    std::string filename;

    float* magx;
    float* magy;
    float* magz;

    int* blockIn;
    int* cellIn;

    int nb[3];
    int nblocks;
    Vector3* bounds_min;
    Vector3* bounds_max;
    Vector3 origin;
    double xlim[2], ylim[2], zlim[2];

    void initFields() override;
    void getFields(ParticleState* state) override;
    void invalidatePositions(ParticleState* state) override;

    bool findBlock(Vector3 position, int &block);
};
