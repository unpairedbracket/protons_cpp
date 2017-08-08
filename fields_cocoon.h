#pragma once

#include <cmath>

#include "fields.h"

struct CocoonField : FieldStructure {
    double B_strength;
    double r_scale, z_scale;

    void initFields() override;
    void getFields(ParticleState* state) override;
};
