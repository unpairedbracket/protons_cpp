#pragma once

#include "fields.h"

struct CocoonField : FieldStructure {
    double leeway = 5.0;

    double B_strength;
    double r_scale, z_scale;
    double min_z_local, max_z_local, max_r_local;

    void initFields() override;
    void getFields(ParticleState* state) override;
    void invalidatePositions(ParticleState* state) override;
};
