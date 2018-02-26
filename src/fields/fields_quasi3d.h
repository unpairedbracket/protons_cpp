#pragma once

#include <string>

#include "fields.h"

struct Q3DField : FieldStructure {
    std::string filename;

    float* magx;
    float* magy;
    float* magz;

    int n_cells[3];

    Vector3 origin;
    Vector3 size;
    double xlim[2], ylim[2], zlim[2];

    double wavelength, frequency, B0, E0;
    double b_mult, e_mult;

    void initFields() override;
    void getFields(ParticleState* state) override;
    void invalidatePositions(ParticleState* state) override;
    void setWavelength(double wavelength);
};
