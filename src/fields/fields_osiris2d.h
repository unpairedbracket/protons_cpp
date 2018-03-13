#pragma once

#include <string>

#include "fields.h"

#include "../util/ndarray.h"

struct Channel {
    int* centres;
    int* starts;
    int* ends;

    void init(int length);
};

struct Osiris2DField : FieldStructure {
    std::string filename;

    ArrayND<2, float> b3;
    Channel* channels;
    int max_channels = 5;
    int channels_found = 0;
    int smooth_half = 50;

    double threshold_field = 0.03;
    int max_gap = 50;
    int min_width = 10;
    int max_fault = 10;
    double max_reduction_factor = 3;
    double reduction_factor = 1.1;
    int max_expansion = 100;

    int n_cells[2];

    Vector3 origin;
    Vector3 size;
    double xlim[2], ylim[2];
    double xlim3[2], ylim3[2], zlim3[2];

    double wavelength, frequency, B0, E0;
    double b_mult, e_mult;

    void initFields() override;
    void getFields(ParticleState* state) override;
    void invalidatePositions(ParticleState* state) override;
    void setWavelength(double wavelength);
    void findChannels();
};

