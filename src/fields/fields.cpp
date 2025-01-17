#include "fields.h"

#include <cmath>

void FieldStructure::initFieldArrays(long N) {
    this->N = N;

    this->E = new Vector3[N];
    this->B = new Vector3[N];

    if(zaxis.x == 0 && zaxis.y == 0) {
        xaxis = {1, 0, 0};
        yaxis = {0, 1, 0};
        zaxis = {0, 0, 1};
    } else {
        double len = sqrt(zaxis.x * zaxis.x
                        + zaxis.y * zaxis.y
                        + zaxis.z * zaxis.z);

        zaxis = {zaxis.x / len, zaxis.y / len, zaxis.z / len}; // z'
        xaxis = {-zaxis.y, zaxis.x, 0}; // x' = z x z'
        len = sqrt(xaxis.x * xaxis.x + xaxis.y * xaxis.y);
        xaxis.x /= len; xaxis.y /= len;

        yaxis = { // y' = z' x x'
                              - zaxis.z * xaxis.y,
            zaxis.z * xaxis.x,
            zaxis.x * xaxis.y - zaxis.y * xaxis.x
        };
    }

    Vector3 tmp = { // xaxis after phi rotation
        xaxis.x * cos(phi) - yaxis.x * sin(phi),
        xaxis.y * cos(phi) - yaxis.y * sin(phi),
        xaxis.z * cos(phi) - yaxis.z * sin(phi)
    };

    yaxis = {
        yaxis.x * cos(phi) + xaxis.x * sin(phi),
        yaxis.y * cos(phi) + xaxis.y * sin(phi),
        yaxis.z * cos(phi) + xaxis.z * sin(phi)
    };

    xaxis = tmp;

    tmp = { // zaxis after theta rotation
        zaxis.x * cos(theta) - yaxis.x * sin(theta),
        zaxis.y * cos(theta) - yaxis.y * sin(theta),
        zaxis.z * cos(theta) - yaxis.z * sin(theta)
    };
    
    yaxis = {
        yaxis.x * cos(theta) + zaxis.x * sin(theta),
        yaxis.y * cos(theta) + zaxis.y * sin(theta),
        yaxis.z * cos(theta) + zaxis.z * sin(theta)
    };

    zaxis = tmp;
}

void FieldStructure::orientBeam(ParticleState* state) {
    if(zaxis.x == 0 && zaxis.y == 0) return;

    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        state->pos[j] = {
            state->pos[j].x * xaxis.x + state->pos[j].y * yaxis.x + state->pos[j].z * zaxis.x,
            state->pos[j].x * xaxis.y + state->pos[j].y * yaxis.y + state->pos[j].z * zaxis.y,
            state->pos[j].x * xaxis.z + state->pos[j].y * yaxis.z + state->pos[j].z * zaxis.z
        };

        state->vel[j] = {
            state->vel[j].x * xaxis.x + state->vel[j].y * yaxis.x + state->vel[j].z * zaxis.x,
            state->vel[j].x * xaxis.y + state->vel[j].y * yaxis.y + state->vel[j].z * zaxis.y,
            state->vel[j].x * xaxis.z + state->vel[j].y * yaxis.z + state->vel[j].z * zaxis.z
        };

        state->intE[j] = {
            state->intE[j].x * xaxis.x + state->intE[j].y * yaxis.x + state->intE[j].z * zaxis.x,
            state->intE[j].x * xaxis.y + state->intE[j].y * yaxis.y + state->intE[j].z * zaxis.y,
            state->intE[j].x * xaxis.z + state->intE[j].y * yaxis.z + state->intE[j].z * zaxis.z
        };

        state->intB[j] = {
            state->intB[j].x * xaxis.x + state->intB[j].y * yaxis.x + state->intB[j].z * zaxis.x,
            state->intB[j].x * xaxis.y + state->intB[j].y * yaxis.y + state->intB[j].z * zaxis.y,
            state->intB[j].x * xaxis.z + state->intB[j].y * yaxis.z + state->intB[j].z * zaxis.z
        };
    }
}
void FieldStructure::deorientBeam(ParticleState* state) {
    if(zaxis.x == 0 && zaxis.y == 0) return;

    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        state->pos[j] = {
            state->pos[j].x * xaxis.x + state->pos[j].y * xaxis.y + state->pos[j].z * xaxis.z,
            state->pos[j].x * yaxis.x + state->pos[j].y * yaxis.y + state->pos[j].z * yaxis.z,
            state->pos[j].x * zaxis.x + state->pos[j].y * zaxis.y + state->pos[j].z * zaxis.z
        };

        state->vel[j] = {
            state->vel[j].x * xaxis.x + state->vel[j].y * xaxis.y + state->vel[j].z * xaxis.z,
            state->vel[j].x * yaxis.x + state->vel[j].y * yaxis.y + state->vel[j].z * yaxis.z,
            state->vel[j].x * zaxis.x + state->vel[j].y * zaxis.y + state->vel[j].z * zaxis.z
        };

        state->intE[j] = {
            state->intE[j].x * xaxis.x + state->intE[j].y * xaxis.y + state->intE[j].z * xaxis.z,
            state->intE[j].x * yaxis.x + state->intE[j].y * yaxis.y + state->intE[j].z * yaxis.z,
            state->intE[j].x * zaxis.x + state->intE[j].y * zaxis.y + state->intE[j].z * zaxis.z
        };

        state->intB[j] = {
            state->intB[j].x * xaxis.x + state->intB[j].y * xaxis.y + state->intB[j].z * xaxis.z,
            state->intB[j].x * yaxis.x + state->intB[j].y * yaxis.y + state->intB[j].z * yaxis.z,
            state->intB[j].x * zaxis.x + state->intB[j].y * zaxis.y + state->intB[j].z * zaxis.z
        };
    }
}
