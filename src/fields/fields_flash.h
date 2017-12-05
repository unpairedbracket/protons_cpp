#pragma once

#include <cmath>
#include <string>

#include "fields.h"

struct FlashField : FieldStructure {
    std::string filename;

    float* bounds;
    float* magx;
    float* magy;
    float* magz;

    int* blockIn;
    int* cellIn;

    int nb[3];
    int nblocks;
    Vector3* bounds_min;
    Vector3* bounds_max;

    void initFields() override;
    void getFields(ParticleState* state) override;
};
