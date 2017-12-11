#include "fields.h"

#include <cmath>

void FieldStructure::initFieldArrays(long N) {
    this->N = N;

    this->E = new Vector3[N];
    this->B = new Vector3[N];

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
    }
}
