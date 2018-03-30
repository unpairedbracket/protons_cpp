#include "invalidateStates.h"

#include "physical_constants.h"

void invalidateStates(ParticleState* particles) {
    #pragma omp parallel for
    for(long j = 0; j < particles->N; j++) {
        if(particles->running[j])
        {
            if(particles->vel[j].x*particles->vel[j].x + particles->vel[j].y*particles->vel[j].y + particles->vel[j].z*particles->vel[j].z > c*c
            ) {
                particles->running[j] = false;
                #pragma omp atomic
                particles->N_running--;
            }
        }
    }
}
