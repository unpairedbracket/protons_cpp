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

    void initFields() override;
    void getFields(ParticleState* state) override;
};
