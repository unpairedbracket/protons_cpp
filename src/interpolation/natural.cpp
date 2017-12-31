#include "natural.h"

#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

typedef std::vector< std::pair< Point, Coord_type > > Points;

void NaturalInterpolator::setSamplePoints(ParticleState* samplePoints) {
    this->sample_points = new Point[samplePoints->N];
    for(long j = 0; j < samplePoints->N; j++) {
        Point p(samplePoints->vel[j].x, samplePoints->vel[j].y);
        this->tri.insert(p);
        this->sample_points[j] = p;
    }
}

void NaturalInterpolator::setSampleValues(ParticleState* sampleValues) {
    for(long j = 0; j < sampleValues->N; j++) {
        pos_x_values.insert(std::make_pair(this->sample_points[j], sampleValues->pos[j].x));
        pos_y_values.insert(std::make_pair(this->sample_points[j], sampleValues->pos[j].y));
        pos_z_values.insert(std::make_pair(this->sample_points[j], sampleValues->pos[j].z));
        vel_x_values.insert(std::make_pair(this->sample_points[j], sampleValues->vel[j].x));
        vel_y_values.insert(std::make_pair(this->sample_points[j], sampleValues->vel[j].y));
        vel_z_values.insert(std::make_pair(this->sample_points[j], sampleValues->vel[j].z));
        //std::cout << "[" << this->sample_points[j] << "] -> [" << sampleValues->pos[j].x << ", " << sampleValues->pos[j].y << ", " << sampleValues->pos[j].z << "], [" << sampleValues->vel[j].x << ", " << sampleValues->vel[j].y << ", " << sampleValues->vel[j].z << "]" << std::endl;
    }

    delete[] this->sample_points;
}

void NaturalInterpolator::initState(ParticleInfo* type) {
    this->interpParticles = this->interpSource->genParticleState(type, this->interpParticles);
}

void NaturalInterpolator::interpolate() {
    typedef CGAL::Data_access< std::map<Point, Coord_type, K::Less_xy_2> > Value_access;
    #pragma omp parallel for
    for(long j = 0; j < this->interpParticles->N; j++) {
        Point p(this->interpParticles->vel[j].x, this->interpParticles->vel[j].y);
        //std::cout << std::endl << "Point: " << p << " -> ";

        std::vector< std::pair<Point, Coord_type> > coords;
        Coord_type norm = CGAL::natural_neighbor_coordinates_2(
                this->tri, p, std::back_inserter(coords)).second;
        //std::cout << "Norm: " << norm << std::endl << "Coords: ";
        //for(auto c : coords) {
            //std::cout << c.first << ", " << c.second << "; ";
        //}
        //std::cout << std::endl;
        Coord_type pos_x = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(pos_x_values));
        Coord_type pos_y = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(pos_y_values));
        Coord_type pos_z = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(pos_z_values));

        //std::cout << "[" << pos_x << ", " << pos_y << ", " << pos_z << "], ";

        Coord_type vel_x = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(vel_x_values));
        Coord_type vel_y = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(vel_y_values));
        Coord_type vel_z = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(vel_z_values));

        //std::cout << "[" << vel_x << ", " << vel_y << ", " << vel_z << "]";
        
        this->interpParticles->pos[j] = {pos_x, pos_y, pos_z};
        this->interpParticles->vel[j] = {vel_x, vel_y, vel_z};
    }
}
