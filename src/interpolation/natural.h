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

struct NaturalInterpolator : Interpolator {
    Delaunay tri;

    std::map<Point, Coord_type, K::Less_xy_2>
        pos_x_values, pos_y_values, pos_z_values,
        vel_x_values, vel_y_values, vel_z_values;

    Point* sample_points;

    void setSamplePoints(ParticleState* samplePoints) override;
    void setSampleValues(ParticleState* sampleValues) override;
    void initState(ParticleInfo* type) override;
    void interpolate() override;
};

