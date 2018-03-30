#pragma once

#include "interpolator.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include "../particles/structs.h"
#include "../sources/source.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
typedef K::FT Coord_type;
typedef K::Point_2 Point;

struct ScatteredLinearInterpolator : Interpolator {
    Delaunay tri;

    std::map<Point, Coord_type, K::Less_xy_2>
        pos_x_values, pos_y_values,
        vel_x_values, vel_y_values, vel_z_values,
        B_x_values, B_y_values, B_z_values,
        E_x_values, E_y_values, E_z_values;

    Vector3* sample_pos;
    Vector3* sample_vel;

    void setSamplePoints(ParticleState* samplePoints) override;
    void setSampleValues(ParticleState* sampleValues) override;

    void interpolate(ParticleState* probeState) override;
};

