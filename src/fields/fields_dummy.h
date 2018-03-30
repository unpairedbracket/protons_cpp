#pragma once

#include "fields.h"

struct DummyField : FieldStructure {
    void initFields() override { this->min_z = 0; };
    void getFields(ParticleState* state) override {};
    void invalidatePositions(ParticleState* state) override {};
};
