#include "fields_quasi3d.h"

#include <H5Cpp.h>
#include <cstdio>
#include <cmath>

#include "../util/math.h"
#include "../util/physical_constants.h"

void Q3DField::initFields() {
    using namespace H5;
    {
        H5File file(this->filename, H5F_ACC_RDONLY);

        double x0 = this->wavelength / (2*pi());

        DataSet xlim_data = file.openDataSet( "/x1_lim" );
        xlim_data.read( this->xlim, PredType::NATIVE_DOUBLE);
        xlim_data.close();
        this->xlim[0] *= x0;
        this->xlim[1] *= x0;
        this->xlim[0] -= this->origin.x;
        this->xlim[1] -= this->origin.x;

        DataSet ylim_data = file.openDataSet( "/x2_lim" );
        ylim_data.read( this->ylim, PredType::NATIVE_DOUBLE);
        ylim_data.close();
        this->ylim[0] *= x0;
        this->ylim[1] *= x0;
        this->ylim[0] -= this->origin.y;
        this->ylim[1] -= this->origin.y;

        DataSet zlim_data = file.openDataSet( "/x3_lim" );
        zlim_data.read( this->zlim, PredType::NATIVE_DOUBLE);
        zlim_data.close();
        this->zlim[0] *= x0;
        this->zlim[1] *= x0;
        this->zlim[0] -= this->origin.z;
        this->zlim[1] -= this->origin.z;

        this->size = {this->xlim[1] - this->xlim[0], this->ylim[1] - this->ylim[0], this->zlim[1] - this->zlim[0]};

        printf("Overall bounds [%f, %f, %f] -> [%f, %f, %f]\n", xlim[0], ylim[0], zlim[0],
                                                 xlim[1], ylim[1], zlim[1]);

        this->min_z = 0;
        double corner_z;
        for(int i = 0; i < 2; i++) {
            for(int j = 0; j < 2; j++) {
                for(int k = 0; k < 2; k++) {
                    corner_z = xlim[i] * zaxis.x + ylim[j] * zaxis.y + zlim[k] * zaxis.z;
                    if(corner_z < this->min_z) this->min_z = corner_z;
                }
            }
        }
        printf("Minimum z: %f\n", this->min_z);

        DataSet magx = file.openDataSet( "/b1" );

        /*
         * Get dataspace of grid variables.
         * TODO Assert on them all being equal?
         */
        DataSpace dataspace = magx.getSpace();

        /*
         * Get the dimension size of each dimension in the dataspace and
         * display them.
         */
        hsize_t mag_dims[3];
        dataspace.getSimpleExtentDims( mag_dims, NULL);
        int nxb = mag_dims[2], nyb = mag_dims[1], nzb = mag_dims[0];

        this->n_cells[0] = nxb;
        this->n_cells[1] = nyb;
        this->n_cells[2] = nzb;

        printf("Magnetic fields are %dx%dx%d", nxb, nyb, nzb);

        int total_cells = nxb*nyb*nzb;

        /*
         * Read data from the file into
         * memory and display the data.
         */
        this->magx = new float[total_cells];
        magx.read( this->magx, PredType::NATIVE_FLOAT);
        magx.close();

        DataSet magy = file.openDataSet( "/b2" );
        this->magy = new float[total_cells];
        magy.read( this->magy, PredType::NATIVE_FLOAT);
        magy.close();

        DataSet magz = file.openDataSet( "/b3" );
        this->magz = new float[total_cells];
        magz.read( this->magz, PredType::NATIVE_FLOAT);
        magz.close();

        file.close();

        #pragma omp parallel for
        for(long j = 0; j < total_cells; j++) {
            this->magx[j] *= this->B0 * this->b_mult;
            this->magy[j] *= this->B0 * this->b_mult;
            this->magz[j] *= this->B0 * this->b_mult;
        }
    }

    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        this->E[j] = { 0.0, 0.0, 0.0 };
        this->B[j] = { 0.0, 0.0, 0.0 };
    }
}

void Q3DField::getFields(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        if(state->running[j]) {
            if(state->pos[j].x >= this->xlim[0]) { //TODO etc.
                float dx = size.x / this->n_cells[0];
                float dy = size.y / this->n_cells[1];
                float dz = size.z / this->n_cells[2];

                float fracx = (state->pos[j].x - this->xlim[0]) / dx;
                float fracy = (state->pos[j].y - this->ylim[0]) / dy;
                float fracz = (state->pos[j].z - this->zlim[0]) / dz;

                int idx_m, idy_m, idz_m;

                bool valid = false;
                if(fracx > 0 && fracx < this->n_cells[0]) {
                    if(fracy > 0 && fracy < this->n_cells[1]) {
                        if(fracz > 0 && fracz < this->n_cells[2]) {
                            valid = true;
                            idx_m = floor(fracx);
                            idy_m = floor(fracy);
                            idz_m = floor(fracz);
                        }
                    }
                }
                float magx = 0, magy = 0, magz = 0;

                if(valid) {
                    int index = ( idz_m * this->n_cells[1] + idy_m ) * this->n_cells[0] + idx_m;
                    magx = this->magx[index];
                    magy = this->magy[index];
                    magz = this->magz[index];
                }


                this->B[j] = {magx, magy, magz};

        } else {
                this->B[j] = { 0.0, 0.0, 0.0 };
            }
        }
    }
}

void Q3DField::invalidatePositions(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            // Testing and tracing to min/max x planes
            double distance = state->pos[j].x - this->xlim[0];
            if(distance < 0) { // We're behind minimum x
                // Are we moving back towards the minimum x plane?
                if(state->vel[j].x < 0) {
                    // If no we're lost
                    state->running[j] = false;
                    #pragma omp atomic
                    state->N_running--;
                    continue;
                } else {
                    // If so let's just go there in one go
                    state->pos[j].x -= distance;
                    state->pos[j].y -= distance * state->vel[j].y / state->vel[j].x;
                    state->pos[j].z -= distance * state->vel[j].z / state->vel[j].x;
                }
            } else {
                // Same for the maximum x plane
                distance = state->pos[j].x - this->xlim[1];
                if(distance > 0) {
                    if(state->vel[j].x > 0) {
                        state->running[j] = false;
                        #pragma omp atomic
                        state->N_running--;
                        continue;
                    } else {
                        state->pos[j].x -= distance;
                        state->pos[j].y -= distance * state->vel[j].y / state->vel[j].x;
                        state->pos[j].z -= distance * state->vel[j].z / state->vel[j].x;
                    }
                }
            }
            // Testing and tracing to min/max y planes
            distance = state->pos[j].y - this->ylim[0];
            if(distance < 0) {
                if(state->vel[j].y < 0) {
                    state->running[j] = false;
                    #pragma omp atomic
                    state->N_running--;
                    continue;
                } else {
                    state->pos[j].x -= distance * state->vel[j].x / state->vel[j].y;
                    state->pos[j].y -= distance;
                    state->pos[j].z -= distance * state->vel[j].z / state->vel[j].y;
                }
            } else {
                distance = state->pos[j].y - this->ylim[1];
                if(distance > 0) {
                    if(state->vel[j].y > 0) {
                        state->running[j] = false;
                        #pragma omp atomic
                        state->N_running--;
                        continue;
                    } else {
                        state->pos[j].x -= distance * state->vel[j].x / state->vel[j].y;
                        state->pos[j].y -= distance;
                        state->pos[j].z -= distance * state->vel[j].z / state->vel[j].y;
                    }
                }
            }
            // Testing and tracing to min/max z planes
            distance = state->pos[j].z - this->zlim[0];
            if(distance < 0) {
                if(state->vel[j].z < 0) {
                    state->running[j] = false;
                    #pragma omp atomic
                    state->N_running--;
                    continue;
                } else {
                    state->pos[j].x -= distance * state->vel[j].x / state->vel[j].z;
                    state->pos[j].y -= distance * state->vel[j].y / state->vel[j].z;
                    state->pos[j].z -= distance;
                }
            } else {
                distance = state->pos[j].z - this->zlim[1];
                if(distance > 0) {
                    if(state->vel[j].z > 0) {
                        state->running[j] = false;
                        #pragma omp atomic
                        state->N_running--;
                        continue;
                    } else {
                        state->pos[j].x -= distance * state->vel[j].x / state->vel[j].z;
                        state->pos[j].y -= distance * state->vel[j].y / state->vel[j].z;
                        state->pos[j].z -= distance;
                    }
                }
            }
        }
    }
}

void Q3DField::setWavelength(double wavelength) {
    this->wavelength = wavelength;
    this->frequency = 2 * pi() * c / wavelength;
    this->B0 = this->frequency * m_e / e;
    this->E0 = this->B0 * c;
}

