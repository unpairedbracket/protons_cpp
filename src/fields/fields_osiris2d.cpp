#include "fields_osiris2d.h"

#include <H5Cpp.h>
#include <cstdio>
#include <cmath>

#include "../util/math.h"
#include "../util/ndarray.h"
#include "../util/h5files.h"
#include "../util/physical_constants.h"

void Osiris2DField::initFields() {
    using namespace H5;
    {
        H5File file(this->filename, H5F_ACC_RDWR);

        double x0 = this->wavelength / (2*pi());

        DataSet xlim_data = file.openDataSet( "/AXIS/AXIS1" );
        xlim_data.read( this->xlim, PredType::NATIVE_DOUBLE);
        xlim_data.close();
        this->xlim[0] *= x0;
        this->xlim[1] *= x0;
        this->xlim[0] -= this->origin.x;
        this->xlim[1] -= this->origin.x;

        DataSet ylim_data = file.openDataSet( "/AXIS/AXIS2" );
        ylim_data.read( this->ylim, PredType::NATIVE_DOUBLE);
        ylim_data.close();
        this->ylim[0] *= x0;
        this->ylim[1] *= x0;
        this->ylim[0] -= this->origin.y;
        this->ylim[1] -= this->origin.y;

        size = {xlim[1]-xlim[0] , ylim[1]-ylim[0], 0};

        printf("Overall bounds [%f, %f] -> [%f, %f]\n", xlim[0], ylim[0],
                                                 xlim[1], ylim[1]);

        DataSet mag = file.openDataSet( "/b3" );

        // OSIRIS HDF5 files are in y,x order, so index with (y, x)
        this->b3 = createArrayFromDataSet<2, float>(mag);
        this->n_cells[0] = this->b3.dimensions[1];
        this->n_cells[1] = this->b3.dimensions[0];
        mag.close();

        this->findChannels();

        Group channel = createOrOpenGroup( file, "/channels" );
        hsize_t channels_dims[2] = {(hsize_t)this->channels_found, (hsize_t)this->n_cells[0]};
        hsize_t channel_dims[1] = {(hsize_t)this->n_cells[0]};
        DataSpace channels_space(2, channels_dims);
        DataSpace single_channel(1, channel_dims);
        DataSet centres = createOrOpenDataSet( file, "/channels/centres", PredType::STD_I32LE, channels_space );
        DataSet starts = createOrOpenDataSet( file, "/channels/starts", PredType::STD_I32LE, channels_space );
        DataSet ends = createOrOpenDataSet( file, "/channels/ends", PredType::STD_I32LE, channels_space );
        hsize_t start[2] = {0, 0};
        hsize_t stride[2] = {1, 1};
        hsize_t count[2] = {1, (hsize_t)this->n_cells[0]};

        for(int i = 0; i < this->channels_found; i++) {
            start[0] = i;
            channels_space.selectHyperslab( H5S_SELECT_SET, count, start, stride, stride );
            centres.write(this->channels[i].centres, PredType::NATIVE_INT, single_channel, channels_space);
            starts.write(this->channels[i].starts, PredType::NATIVE_INT, single_channel, channels_space);
            ends.write(this->channels[i].ends, PredType::NATIVE_INT, single_channel, channels_space);
        }
        centres.close();
        starts.close();
        ends.close();

        file.close();
    }

    int total_cells = this->b3.dimensions[0] * this->b3.dimensions[1];

    #pragma omp parallel for
    for(long j = 0; j < total_cells; j++) {
        this->b3[j] *= this->B0 * this->b_mult;
    }

    this->min_z = 0;
    double corner_z;
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            for(int k = 0; k < 2; k++) {
                corner_z = xlim3[i] * zaxis.x + ylim3[j] * zaxis.y + zlim3[k] * zaxis.z;
                if(corner_z < this->min_z) this->min_z = corner_z;
            }
        }
    }
    printf("Minimum z: %f\n", this->min_z);


    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        this->E[j] = { 0.0, 0.0, 0.0 };
        this->B[j] = { 0.0, 0.0, 0.0 };
    }
}

void Osiris2DField::findChannels() {
    using std::max;
    using std::min;
    this->channels = new Channel[this->max_channels];
    for(int c = 0; c < this->max_channels; c++) {
        this->channels[c].init(this->n_cells[0]);
    }

    // Produce a Gaussian smoothing kernel
    int smoothing = 2*this->smooth_half + 1;
    double* kern = new double[smoothing];
    double sum = 0;
    for(int j = 0; j < smoothing; j++) {
        double x = 4.0 * ((double)(j - smooth_half)) / ((double)(smooth_half));
        kern[j] = exp(-(x*x));
        sum += kern[j];
    }
    for(int j = 0; j < smoothing; j++) {
        kern[j] /= sum;
    }

    int max_height = 0;
    int min_y = 0;
    int max_y = 0;
    this->channels_found = 0;

    double* current_row = new double[n_cells[1]];
    for(int j = 0; j < this->n_cells[0]; j++) {
        // Calculate a smoothed version of this row by convolution with kern
        for(int k = 0; k < this->n_cells[1]; k++) {
            // k + l - smooth_num needs to be >= 0 and < n_cells[1]
            // l >= 0 and < smoothing
            // l >= smooth_num - k and >= 0 and < n_cells[1] + smooth_num - k and < smoothing
            int left = max(0, smooth_half - k);
            int right = min(smoothing, (int)(n_cells[1] + smooth_half - k));
            current_row[k] = 0;
            for(int l = left; l < right; l++) {
                current_row[k] += this->b3(k+l-smooth_half, j) * kern[l];
            }
        }

        int found = 0;
        double threshold = this->threshold_field;
        int tries = 0;

        while(found < 1 && threshold >= this->threshold_field / this->max_reduction_factor) {
            int stage = 1;
            int start_pixel = -1;
            int gap_start = -1;
            int centre_pixel = -1;
            int gap_end = -1;
            int end_pixel = -1;

            for(int pixidx = 0; pixidx < this->n_cells[1]; pixidx++) {
                double pixel = current_row[pixidx];
                if(stage == 1) {
                    // Look for > threshold
                    if(pixel > threshold) {
                        start_pixel = pixidx;
                        stage = 2;
                    }
                } else if(stage == 2) {
                    // Look for start of gap
                    if(pixel < threshold) {
                        gap_start = pixidx;
                        stage = 3;
                    }
                } else if(stage == 3) {
                    // In gap, look for either jumping out of gap or following
                    // through to < threshold
                    if(pixel > threshold) {
                        if(pixidx > gap_start + this->max_fault) {
                            start_pixel = pixidx;
                        }
                        stage = 2;
                    } else if(-pixel > threshold) {
                        gap_end = pixidx;
                        stage = 4;
                    } else if(pixidx > gap_start + this->max_gap) {
                        stage = 1;
                    }
                } else if(stage == 4) {
                    if(-pixel < threshold) {
                        end_pixel = pixidx;
                        stage = 5;
                    }
                } else if(stage == 5) {
                    if(-pixel > threshold) {
                        stage = 4;
                    } else if(pixel > threshold || pixidx > end_pixel + this->max_fault || pixidx == this->n_cells[1] - 1) {
                        centre_pixel = (gap_start + gap_end) / 2;
                        int gap_width = gap_end - gap_start;
                        if(centre_pixel - start_pixel > this->min_width && end_pixel - centre_pixel > this->min_width && gap_width < this->max_gap) {
                            this->channels[found].centres[j] = centre_pixel;
                            this->channels[found].starts[j] = start_pixel;
                            this->channels[found].ends[j] = end_pixel;
                            found = found + 1;
                            int larger_width = max(end_pixel - centre_pixel, centre_pixel - start_pixel);
                            max_height = max(max_height, larger_width);
                            max_y = max(max_y, centre_pixel + larger_width);
                            min_y = min(min_y, centre_pixel - larger_width);
                            if(found >= this->max_channels) break;
                        }
                        stage = 1;
                    }
                }
            }
            if(found < 1) {
                threshold = threshold / this->reduction_factor;
                tries++;
            }
        }

        //Expand the limits of the channels a bit
        if(found > 0) {
            this->channels[0].starts[j] = max(0, this->channels[0].starts[j] - this->max_expansion);
            for(int k = 1; k < found; k++) {
                double midpoint = 0.5 * (this->channels[k-1].ends[j] + this->channels[k].starts[j]);
                this->channels[k-1].ends[j] = min((int)(floor(midpoint) - (midpoint == trunc(midpoint))), this->channels[k-1].ends[j] + this->max_expansion);
                this->channels[k].starts[j] = max((int)(ceil(midpoint) + (midpoint == trunc(midpoint))), this->channels[k-1].starts[j] - this->max_expansion);
            }
            this->channels[found-1].ends[j] = min(this->n_cells[1] - 1, this->channels[found-1].ends[j] + this->max_expansion);
        } else {
            int centre = this->n_cells[1] / 2;
            this->channels[0].starts[j] = max(0, centre - this->max_expansion);
            this->channels[0].ends[j] = min(this->n_cells[1] - 1, centre + this->max_expansion);
            this->channels[0].centres[j] = centre;
        }
        this->channels_found = max(this->channels_found, found);
    }
    printf("Found a total of %d channels.\n", this->channels_found);
    printf("Min/Max z: %d\n", max_height);
    printf("Min/Max y: %d, %d\n", min_y, max_y);

    float dx = this->size.x / this->n_cells[0];
    float dy = this->size.y / this->n_cells[1];

    float min_y_coord = this->ylim[0] + min_y * dy;
    float max_y_coord = this->ylim[0] + (max_y + 1) * dy;
    float max_absz = (max_height+1) * dy;

    this->xlim3[0] = this->xlim[0]; this->xlim3[1] = this->xlim[1];
    this->ylim3[0] = min_y_coord;   this->ylim3[1] = max_y_coord;
    this->zlim3[0] = -max_absz;     this->zlim3[1] = max_absz;
}

void Osiris2DField::getFields(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        if(state->running[j]) {
            float dx = size.x / this->n_cells[0];
            float dy = size.y / this->n_cells[1];

            float fracx = (state->pos[j].x - this->xlim[0]) / dx;
            float fracy = (state->pos[j].y - this->ylim[0]) / dy;

            int idx_m, idy_m;

            bool valid = false;
            if(fracx > 0.5 && fracx < this->n_cells[0] - 0.5) {
                if(fracy > 0.5 && fracy < this->n_cells[1] - 0.5) {
                    valid = true;
                    idx_m = floor(fracx);
                    idy_m = floor(fracy);
                }
            }

            this->B[j] = { 0, 0, 0 };

            if(valid) {
                double By_m = 0, By_p = 0;
                double Bz_m = 0, Bz_p = 0;
                for(int i = 0; i < max_channels; i++) {
                    int centre = this->channels[i].centres[idx_m];
                    int start = this->channels[i].starts[idx_m];
                    int end = this->channels[i].ends[idx_m];
                    if(centre >= 0){
                        double centre_y = this->ylim[0] + (0.5 + (double)centre) * dy;
                        double y_rel = state->pos[j].y - centre_y;
                        double z_rel = state->pos[j].z;
                        double r = sqrt(y_rel * y_rel + z_rel*z_rel);
                        int idr_p = centre + floor(r / dy);
                        int idr_m = centre - floor(r / dy);

                        double cos_phi = y_rel / r;
                        double sin_phi = z_rel / r;

                        double proportionRight = 0.5 * (1.0 + cos_phi);
                        double proportionLeft = 0.5 * (1.0 - cos_phi);

                        double b_phi = 0;
                        if(idr_p < end) {
                            b_phi += proportionRight * this->b3(idr_p, idx_m);
                        }
                        if(idr_m > start) {
                            b_phi -= proportionLeft * this->b3(idr_m, idx_m);
                        }
                        Bz_m += b_phi * cos_phi;
                        By_m -= b_phi * sin_phi;
                    }
                }
                int idx_p = idx_m + 1;
                for(int i = 0; i < max_channels; i++) {
                    int centre = this->channels[i].centres[idx_p];
                    int start = this->channels[i].starts[idx_p];
                    int end = this->channels[i].ends[idx_p];
                    if(centre >= 0){
                        double centre_y = this->ylim[0] + (0.5 + (double)centre) * dy;
                        double y_rel = state->pos[j].y - centre_y;
                        double z_rel = state->pos[j].z;
                        double r = sqrt(y_rel * y_rel + z_rel*z_rel);
                        int idr_p = centre + floor(r / dy);
                        int idr_m = centre - floor(r / dy);

                        double cos_phi = y_rel / r;
                        double sin_phi = z_rel / r;

                        double proportionRight = 0.5 * (1.0 + cos_phi);
                        double proportionLeft = 0.5 * (1.0 - cos_phi);

                        double b_phi = 0;
                        if(idr_p < end) {
                            b_phi += proportionRight * this->b3(idr_p, idx_p);
                        }
                        if(idr_m > start) {
                            b_phi -= proportionLeft * this->b3(idr_m, idx_p);
                        }
                        Bz_p += b_phi * cos_phi;
                        By_p -= b_phi * sin_phi;
                    }
                }
                this->B[j] = { 0.0, By_m + (fracx - (double)idx_m) * (By_p - By_m), Bz_m + (fracx - (double)idx_m) * (Bz_p - Bz_m) };
                //this->B[j] = { 0.0, By_m, Bz_m };
            }
        } else {
            this->B[j] = { 0.0, 0.0, 0.0 };
        }
    }
}

void Osiris2DField::invalidatePositions(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            // Testing and tracing to min/max x planes
            double distance = state->pos[j].x - this->xlim3[0];
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
                distance = state->pos[j].x - this->xlim3[1];
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
            distance = state->pos[j].y - this->ylim3[0];
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
                distance = state->pos[j].y - this->ylim3[1];
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
            distance = state->pos[j].z - this->zlim3[0];
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
                distance = state->pos[j].z - this->zlim3[1];
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

void Osiris2DField::setWavelength(double wavelength) {
    this->wavelength = wavelength;
    this->frequency = 2 * pi() * c / wavelength;
    this->B0 = this->frequency * m_e / e;
    this->E0 = this->B0 * c;
}

void Channel::init(int length) {
    printf("Channel length: %d\n", length);
    this->centres = new int[length];
    this->starts = new int[length];
    this->ends = new int[length];

    for(int j = 0; j < length; j++) {
        this->centres[j] = -1;
        this->starts[j] = -1;
        this->ends[j] = -1;
    }
}

