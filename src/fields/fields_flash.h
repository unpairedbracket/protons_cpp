#pragma once

#include <cmath>
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

    void initFields() override;
    bool findBlock(Vector3 position, int &block);
    void getFields(ParticleState* state) override;
};
