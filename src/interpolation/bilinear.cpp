#include "bilinear.h"

void BilinearInterpolator::setSamplePoints(ParticleState* samplePoints) {
    initParticleState(&sampleVelocities, samplePoints->particleInfo, samplePoints->N);
    for(long j = 0; j < sampleVelocities.N; j++) {
        sampleVelocities.pos[j] = samplePoints->pos[j];
        sampleVelocities.vel[j] = samplePoints->vel[j];
    }
}

void BilinearInterpolator::setSampleValues(ParticleState* sampleValues) {
    initParticleState(&sampleDeflections, sampleValues->particleInfo, sampleValues->N);
    printf("Sample Deflections N = %ld\n", sampleDeflections.N);
    for(long j = 0; j < sampleDeflections.N; j++) {
        int xpos = j / 80;
        int ypos = j % 80;
        sampleDeflections.vel[j] = {
            sampleValues->vel[j].x - sampleVelocities.vel[j].x,
            sampleValues->vel[j].y - sampleVelocities.vel[j].y,
            sampleValues->vel[j].z - sampleVelocities.vel[j].z
        };
        sampleDeflections.pos[j] = {
            (   sampleValues->pos[j].x -    sampleValues->pos[j].z *    sampleValues->vel[j].x /    sampleValues->vel[j].z)
          - (sampleVelocities.pos[j].x - sampleVelocities.pos[j].z * sampleVelocities.vel[j].x / sampleVelocities.vel[j].z),
            (   sampleValues->pos[j].y -    sampleValues->pos[j].z *    sampleValues->vel[j].y /    sampleValues->vel[j].z)
          - (sampleVelocities.pos[j].y - sampleVelocities.pos[j].z * sampleVelocities.vel[j].y / sampleVelocities.vel[j].z),
            0
        };
        sampleDeflections.intB[j] = sampleValues->intB[j];
        sampleDeflections.intE[j] = sampleValues->intE[j];
    }
}

void BilinearInterpolator::interpolate(ParticleState* probeState) {
    #pragma omp parallel for
    for(long j = 0; j < probeState->N; j++) {
        double x_interp = probeState->vel[j].x / probeState->vel[j].z;
        double y_interp = probeState->vel[j].y / probeState->vel[j].z;

        probeState->pos[j].x -= x_interp * probeState->pos[j].z;
        probeState->pos[j].y -= y_interp * probeState->pos[j].z;
        probeState->pos[j].z = 0;

        double x_frac = (x_interp / this->size[0] + 0.5) * this->n_cells[0];
        double y_frac = (y_interp / this->size[1] + 0.5) * this->n_cells[1];

        if(x_frac > 0.5 && x_frac < this->n_cells[0] - 0.5 && y_frac > 0.5 && y_frac < this->n_cells[1] - 0.5) {

            int xm = floor(x_frac - 0.5);
            int ym = floor(y_frac - 0.5);

            double xf = (x_frac-0.5) - xm;
            double yf = (y_frac-0.5) - ym;

            int idmm = xm * this->n_cells[1] + ym;
            int idmp = xm * this->n_cells[1] + ym + 1;
            int idpm = (xm+1) * this->n_cells[1] + ym;
            int idpp = (xm+1) * this->n_cells[1] + ym + 1;

            Vector3 d_pos = {
                (1-xf) * ((1-yf)*sampleDeflections.pos[idmm].x+yf*sampleDeflections.pos[idmp].x)
                 + xf  * ((1-yf)*sampleDeflections.pos[idpm].x+yf*sampleDeflections.pos[idpp].x),
                (1-xf) * ((1-yf)*sampleDeflections.pos[idmm].y+yf*sampleDeflections.pos[idmp].y)
                 + xf  * ((1-yf)*sampleDeflections.pos[idpm].y+yf*sampleDeflections.pos[idpp].y),
                (1-xf) * ((1-yf)*sampleDeflections.pos[idmm].z+yf*sampleDeflections.pos[idmp].z)
                 + xf  * ((1-yf)*sampleDeflections.pos[idpm].z+yf*sampleDeflections.pos[idpp].z)
            };

            Vector3 d_vel = {
                (1-xf) * ((1-yf)*sampleDeflections.vel[idmm].x+yf*sampleDeflections.vel[idmp].x)
                 + xf  * ((1-yf)*sampleDeflections.vel[idpm].x+yf*sampleDeflections.vel[idpp].x),
                (1-xf) * ((1-yf)*sampleDeflections.vel[idmm].y+yf*sampleDeflections.vel[idmp].y)
                 + xf  * ((1-yf)*sampleDeflections.vel[idpm].y+yf*sampleDeflections.vel[idpp].y),
                (1-xf) * ((1-yf)*sampleDeflections.vel[idmm].z+yf*sampleDeflections.vel[idmp].z)
                 + xf  * ((1-yf)*sampleDeflections.vel[idpm].z+yf*sampleDeflections.vel[idpp].z)
            };

            probeState->pos[j].x += d_pos.x;
            probeState->pos[j].y += d_pos.y;

            probeState->vel[j].x += d_vel.x;
            probeState->vel[j].y += d_vel.y;
            probeState->vel[j].z += d_vel.z;

            probeState->intB[j] = {
                (1-xf) * ((1-yf)*sampleDeflections.intB[idmm].x+yf*sampleDeflections.intB[idmp].x)
                 + xf  * ((1-yf)*sampleDeflections.intB[idpm].x+yf*sampleDeflections.intB[idpp].x),
                (1-xf) * ((1-yf)*sampleDeflections.intB[idmm].y+yf*sampleDeflections.intB[idmp].y)
                 + xf  * ((1-yf)*sampleDeflections.intB[idpm].y+yf*sampleDeflections.intB[idpp].y),
                (1-xf) * ((1-yf)*sampleDeflections.intB[idmm].z+yf*sampleDeflections.intB[idmp].z)
                 + xf  * ((1-yf)*sampleDeflections.intB[idpm].z+yf*sampleDeflections.intB[idpp].z)
            };

            probeState->intE[j] = {
                (1-xf) * ((1-yf)*sampleDeflections.intE[idmm].x+yf*sampleDeflections.intE[idmp].x)
                 + xf  * ((1-yf)*sampleDeflections.intE[idpm].x+yf*sampleDeflections.intE[idpp].x),
                (1-xf) * ((1-yf)*sampleDeflections.intE[idmm].y+yf*sampleDeflections.intE[idmp].y)
                 + xf  * ((1-yf)*sampleDeflections.intE[idpm].y+yf*sampleDeflections.intE[idpp].y),
                (1-xf) * ((1-yf)*sampleDeflections.intE[idmm].z+yf*sampleDeflections.intE[idmp].z)
                 + xf  * ((1-yf)*sampleDeflections.intE[idpm].z+yf*sampleDeflections.intE[idpp].z)
            };
        }
    }
}
