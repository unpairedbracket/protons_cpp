#include "fields_cylindrical.h"

#include <H5Cpp.h>
#include <cmath>

void CylindricalField::initFields() {
    using namespace H5;
    {
        H5File file(this->filename, H5F_ACC_RDONLY);


        DataSet rmax_data = file.openDataSet( "/rmax" );
        rmax_data.read( &(this->rmax), PredType::NATIVE_DOUBLE);
        rmax_data.close();

        DataSet zlim_data = file.openDataSet( "/zlim" );
        zlim_data.read( this->zlim, PredType::NATIVE_DOUBLE);
        zlim_data.close();
        this->zlim[0] -= this->origin_z;
        this->zlim[1] -= this->origin_z;
        this->sizez = (this->zlim[1] - this->zlim[0]);

        printf("Field dimensions: r = [0 -> %e]; z = [%e -> %e]", this->rmax, this->zlim[0], this->zlim[1]);

        DataSet magphi_data = file.openDataSet( "/Bphi" );

        DataSpace dataspace = magphi_data.getSpace();

        hsize_t dims[2];
        dataspace.getSimpleExtentDims( dims, NULL);
        int nzb = dims[1], nrb = dims[0];

        this->nb[0] = nrb;
        this->nb[1] = nzb;

        printf("Magnetic fields are %dx%d\n", nrb, nzb);

        /*
         * Read data from the file into
         * memory and display the data.
         */
        this->magphi = new double[nrb*nzb];
        magphi_data.read( this->magphi, PredType::NATIVE_DOUBLE);
        magphi_data.close();

        DataSet eler_data = file.openDataSet( "/Er" );
        this->eler = new double[nrb*nzb];
        eler_data.read( this->eler, PredType::NATIVE_DOUBLE);
        eler_data.close();

        DataSet elez_data = file.openDataSet( "/Ez" );
        this->elez = new double[nrb*nzb];
        elez_data.read( this->elez, PredType::NATIVE_DOUBLE);
        elez_data.close();

        file.close();
    }

    #pragma omp parallel for
    for(int i = 0; i < this->nb[0] * this->nb[1]; i++) {
        this->magphi[i] *= this->Bscale;
        this->eler[i] *= this->Escale;
        this->elez[i] *= this->Escale;
    }

    this->min_z = std::min(this->zlim[0] * zaxis.z, this->zlim[1] * zaxis.z) - this->rmax*sqrt(zaxis.x*zaxis.x + zaxis.y * zaxis.y);
    printf("Min z: %f", this->min_z);

    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        this->E[j] = { 0.0, 0.0, 0.0 };
        this->B[j] = { 0.0, 0.0, 0.0 };
    }
}

void CylindricalField::getFields(ParticleState* state) {
    double dr = this->rmax / this->nb[0];
    double dz = this->sizez / this->nb[1];

    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        if(state->running[j]) {
            double r = sqrt(state->pos[j].x*state->pos[j].x + state->pos[j].y * state->pos[j].y);
            double z = state->pos[j].z;

            if(r <= (this->rmax + dr/2) && z >= (this->zlim[0]-dz/2) && z <= (this->zlim[1]+dz/2)) {
                double fracr = r / dr;
                double fracz = (z - this->zlim[0]) / dz;

                int idr_m = floor(fracr-0.5);
                int idz_m = floor(fracz-0.5);

                double del_r = (fracr-0.5) - idr_m;
                double del_z = (fracz-0.5) - idr_m;
                double idel_r = 1-del_r;
                double idel_z = 1-del_z;

                double magphi_mz;
                double magphi_pz;

                double eler_mz;
                double eler_pz;

                double elez_mz;
                double elez_pz;

                int nz = nb[1];
                auto index = [nz](int idr, int idz) -> int { return nz * idr + idz; };

                if(idz_m < 0) {
                    magphi_mz = 0;
                    eler_mz = 0;
                    elez_mz = 0;
                    if(idr_m < 0) {
                        magphi_pz = magphi[index(0, 0)] * (2*del_r - 1);
                        eler_pz = eler[index(0, 0)] * (2*del_r - 1);
                        elez_pz = elez[index(0, 0)] * (2*del_r - 1);
                    } else if(idr_m + 1 > this->nb[0] - 1) {
                        magphi_pz = magphi[index(idr_m, 0)] * (1-del_r);
                        eler_pz = eler[index(idr_m, 0)] * (1-del_r);
                        elez_pz = elez[index(idr_m, 0)] * (1-del_r);
                    } else {
                        magphi_pz = magphi[index(idr_m, 0)] * (1-del_r) + magphi[index(idr_m+1, 0)] * del_r;
                        eler_pz = eler[index(idr_m, 0)] * (1-del_r) + eler[index(idr_m+1, 0)] * del_r;
                        elez_pz = elez[index(idr_m, 0)] * (1-del_r) + elez[index(idr_m+1, 0)] * del_r;
                    }
                } else if(idz_m + 1 > this->nb[1] - 1) {
                    magphi_pz = 0;
                    eler_pz = 0;
                    elez_pz = 0;
                    if(idr_m < 0) {
                        magphi_mz = magphi[index(0, this->nb[1]-1)] * (2*del_r - 1);
                        eler_mz = eler[index(0, this->nb[1]-1)] * (2*del_r - 1);
                        elez_mz = elez[index(0, this->nb[1]-1)] * (2*del_r - 1);
                    } else if(idr_m + 1 > this->nb[0] - 1) {
                        magphi_mz = magphi[index(idr_m, this->nb[1]-1)] * (1-del_r);
                        eler_mz = eler[index(idr_m, this->nb[1]-1)] * (1-del_r);
                        elez_mz = elez[index(idr_m, this->nb[1]-1)] * (1-del_r);
                    } else {
                        magphi_mz = magphi[index(idr_m, this->nb[1]-1)] * (1-del_r) + magphi[index(idr_m+1, this->nb[1]-1)] * del_r;
                        eler_mz = eler[index(idr_m, this->nb[1]-1)] * (1-del_r) + eler[index(idr_m+1, this->nb[1]-1)] * del_r;
                        elez_mz = elez[index(idr_m, this->nb[1]-1)] * (1-del_r) + elez[index(idr_m+1, this->nb[1]-1)] * del_r;
                    }
                } else {
                    if(idr_m < 0) {
                        magphi_mz = magphi[index(0, idz_m  )] * (2*del_r - 1);
                        magphi_pz = magphi[index(0, idz_m+1)] * (2*del_r - 1);
                        eler_mz = eler[index(0, idz_m  )] * (2*del_r - 1);
                        eler_pz = eler[index(0, idz_m+1)] * (2*del_r - 1);
                        elez_mz = elez[index(0, idz_m  )] * (2*del_r - 1);
                        elez_pz = elez[index(0, idz_m+1)] * (2*del_r - 1);
                    } else if(idr_m + 1 > this->nb[0] - 1) {
                        magphi_mz = magphi[index(idr_m, idz_m  )] * (1-del_r);
                        magphi_pz = magphi[index(idr_m, idz_m+1)] * (1-del_r);
                        eler_mz = eler[index(idr_m, idz_m  )] * (1-del_r);
                        eler_pz = eler[index(idr_m, idz_m+1)] * (1-del_r);
                        elez_mz = elez[index(idr_m, idz_m  )] * (1-del_r);
                        elez_pz = elez[index(idr_m, idz_m+1)] * (1-del_r);
                    } else {
                        magphi_mz = magphi[index(idr_m, idz_m  )] * (1-del_r) + magphi[index(idr_m+1, idz_m  )] * del_r;
                        magphi_pz = magphi[index(idr_m, idz_m+1)] * (1-del_r) + magphi[index(idr_m+1, idz_m+1)] * del_r;
                        eler_mz = eler[index(idr_m, idz_m  )] * (1-del_r) + eler[index(idr_m+1, idz_m  )] * del_r;
                        eler_pz = eler[index(idr_m, idz_m+1)] * (1-del_r) + eler[index(idr_m+1, idz_m+1)] * del_r;
                        elez_mz = elez[index(idr_m, idz_m  )] * (1-del_r) + elez[index(idr_m+1, idz_m  )] * del_r;
                        elez_pz = elez[index(idr_m, idz_m+1)] * (1-del_r) + elez[index(idr_m+1, idz_m+1)] * del_r;
                    }
                }

                double Bp = idel_z * magphi_mz + del_z * magphi_pz;
                double Er = idel_z * eler_mz + del_z * eler_pz;
                double Ez = idel_z * elez_mz + del_z * elez_pz;

                this->B[j] = { Bp * -state->pos[j].y/r, Bp * state->pos[j].x/r, 0 };
                this->E[j] = { Er * state->pos[j].x/r, Er * state->pos[j].y/r, Ez };
            } else {
                this->B[j] = { 0.0, 0.0, 0.0 };
                this->E[j] = { 0.0, 0.0, 0.0 };
            }
        }
    }
}

void CylindricalField::invalidatePositions(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            Vector3 pos = state->pos[j];
            Vector3 vel = state->vel[j];

            if(((pos.z < this->zlim[0] || pos.z > this->zlim[1]) && vel.z * pos.z > 0)
            || (pos.x*pos.x + pos.y*pos.y > this->rmax*this->rmax && vel.x * pos.x + vel.y * pos.y > - sqrt( (vel.x*vel.x + vel.y*vel.y) * (pos.x*pos.x + pos.y*pos.y - this->rmax*this->rmax)))
            ) {
                state->running[j] = false;
                #pragma omp atomic
                state->N_running--;
            }
        }
    }
}

