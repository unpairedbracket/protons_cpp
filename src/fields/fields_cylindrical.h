#pragma once

#include <string>

#include "fields.h"

struct CylindricalField : FieldStructure {
    std::string filename;

    double* magphi;
    double* eler;
    double* elez;

    double Bscale;
    double Escale;

    double origin_z;
    int nb[3];
    double rmax, zlim[2], sizez;

    void initFields() override;
    void getFields(ParticleState* state) override;
    void invalidatePositions(ParticleState* state) override;
};
