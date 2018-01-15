#include "linear.h"

#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

typedef std::vector< std::pair< Point, Coord_type > > Points;
typedef K::Triangle_2 Triangle;
typedef Delaunay::Face_handle Face;
typedef K::Vector_2 Vector;
typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<K> Triangle_coordinates;

void LinearInterpolator::setSamplePoints(ParticleState* samplePoints) {
    this->sample_points = new Point[samplePoints->N];
    for(long j = 0; j < samplePoints->N; j++) {
        Point p(samplePoints->vel[j].x, samplePoints->vel[j].y);
        this->tri.insert(p);
        this->sample_points[j] = p;
    }
}

void LinearInterpolator::setSampleValues(ParticleState* sampleValues) {
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

void LinearInterpolator::interpolate() {
    typedef CGAL::Data_access< std::map<Point, Coord_type, K::Less_xy_2> > Value_access;
    #pragma omp parallel for
    for(long j = 0; j < this->interpParticles->N; j++) {
        Point p(this->interpParticles->vel[j].x, this->interpParticles->vel[j].y);
        //std::cout << std::endl << "Point: " << p << " -> ";

        std::vector< std::pair<Point, Coord_type> > coords;
        Face f = this->tri.inexact_locate(p);
        if(!this->tri.is_infinite(f)) {
            Triangle t = this->tri.triangle(f);

            Triangle_coordinates bary(t[0], t[1], t[2]);

            std::vector< std::pair<Point, Coord_type> > coords;
            std::vector< Coord_type > baryCoords;
            bary(p, baryCoords);

            coords.push_back(std::make_pair(t[0], baryCoords[0]));
            coords.push_back(std::make_pair(t[1], baryCoords[1]));
            coords.push_back(std::make_pair(t[2], baryCoords[2]));
            Coord_type norm = baryCoords[0] + baryCoords[1] + baryCoords[2];

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
}
