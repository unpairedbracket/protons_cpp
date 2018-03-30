#include "natural.h"

#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

typedef std::vector< std::pair< Point, Coord_type > > Points;

void NaturalInterpolator::setSamplePoints(ParticleState* samplePoints) {
    this->sample_pos = new Vector3[samplePoints->N];
    this->sample_vel = new Vector3[samplePoints->N];
    for(long j = 0; j < samplePoints->N; j++) {
        Point p(samplePoints->vel[j].x / samplePoints->vel[j].z, samplePoints->vel[j].y / samplePoints->vel[j].z);
        this->tri.insert(p);
        this->sample_pos[j] = samplePoints->pos[j];
        this->sample_vel[j] = samplePoints->vel[j];
    }
}


void NaturalInterpolator::setSampleValues(ParticleState* sampleValues) {
    for(long j = 0; j < sampleValues->N; j++) {
        Point p(this->sample_vel[j].x / this->sample_vel[j].z, this->sample_vel[j].y / this->sample_vel[j].z);
        Vector3 v_deflection = {
            sampleValues->vel[j].x - this->sample_vel[j].x,
            sampleValues->vel[j].y - this->sample_vel[j].y,
            sampleValues->vel[j].z - this->sample_vel[j].z
        };
        Vector3 r_deflection = {
            (sampleValues->pos[j].x - sampleValues->pos[j].z * sampleValues->vel[j].x / sampleValues->vel[j].z) 
          - (this->sample_pos[j].x  - this->sample_pos[j].z  * this->sample_vel[j].x  / this->sample_vel[j].z ),
            (sampleValues->pos[j].y - sampleValues->pos[j].z * sampleValues->vel[j].y / sampleValues->vel[j].z) 
          - (this->sample_pos[j].y  - this->sample_pos[j].z  * this->sample_vel[j].y  / this->sample_vel[j].z ),
          0
        };
        pos_x_values.insert(std::make_pair(p, r_deflection.x));
        pos_y_values.insert(std::make_pair(p, r_deflection.y));

        vel_x_values.insert(std::make_pair(p, v_deflection.x));
        vel_y_values.insert(std::make_pair(p, v_deflection.y));
        vel_z_values.insert(std::make_pair(p, v_deflection.z));

        B_x_values.insert(std::make_pair(p, sampleValues->intB[j].x));
        B_y_values.insert(std::make_pair(p, sampleValues->intB[j].y));
        B_z_values.insert(std::make_pair(p, sampleValues->intB[j].z));

        E_x_values.insert(std::make_pair(p, sampleValues->intE[j].x));
        E_y_values.insert(std::make_pair(p, sampleValues->intE[j].y));
        E_z_values.insert(std::make_pair(p, sampleValues->intE[j].z));
    }

    delete[] this->sample_pos;
    delete[] this->sample_vel;
}

void NaturalInterpolator::interpolate(ParticleState* probeState) {
    typedef CGAL::Data_access< std::map<Point, Coord_type, K::Less_xy_2> > Value_access;
    #pragma omp parallel for
    for(long j = 0; j < probeState->N; j++) {
        Point p(probeState->vel[j].x / probeState->vel[j].z, probeState->vel[j].y / probeState->vel[j].z);

        std::vector< std::pair<Point, Coord_type> > coords;
        auto result = CGAL::natural_neighbor_coordinates_2(
                this->tri, p, std::back_inserter(coords));
        if(result.third) {
            Coord_type norm = result.second;
            Coord_type d_pos_x = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(pos_x_values));
            Coord_type d_pos_y = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(pos_y_values));

            Coord_type d_vel_x = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(vel_x_values));
            Coord_type d_vel_y = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(vel_y_values));
            Coord_type d_vel_z = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(vel_z_values));

            Coord_type B_x = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(B_x_values));
            Coord_type B_y = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(B_y_values));
            Coord_type B_z = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(B_z_values));

            Coord_type E_x = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(E_x_values));
            Coord_type E_y = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(E_y_values));
            Coord_type E_z = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(E_z_values));

            probeState->pos[j].x -= p.x() * probeState->pos[j].z;
            probeState->pos[j].y -= p.y() * probeState->pos[j].z;
            probeState->pos[j].z = 0;

            probeState->pos[j].x += d_pos_x;
            probeState->pos[j].y += d_pos_y;

            probeState->vel[j].x += d_vel_x;
            probeState->vel[j].y += d_vel_y;
            probeState->vel[j].z += d_vel_z;

            probeState->intB[j].x += B_x;
            probeState->intB[j].y += B_y;
            probeState->intB[j].z += B_z;

            probeState->intE[j].x += E_x;
            probeState->intE[j].y += E_y;
            probeState->intE[j].z += E_z;
        }
    }
}
