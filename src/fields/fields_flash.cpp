#include "fields_flash.h"

#include <H5Cpp.h>
#include <cassert>

void FlashField::initFields() {
    using namespace H5;
    try {
        Exception::dontPrint();

        H5File file(this->filename, H5F_ACC_RDONLY);

        DataSet bounding_boxes = file.openDataSet( "/bounding box" );

        /*
         * Get dataspace of the dataset.
         */
        DataSpace dataspace = bounding_boxes.getSpace();

        /*
         * Get the dimension size of each dimension in the dataspace and
         * display them.
         */
        hsize_t dims_out[3];
        dataspace.getSimpleExtentDims( dims_out, NULL);
        this->nblocks = dims_out[0];

        float* bounds = new float[nblocks*3*2];
        auto index = [](int block, int dir, int minmax) {return 6*block + 2 * dir + minmax;};

        /*
         * Read data from hyperslab in the file into the hyperslab in
         * memory and display the data.
         */
        bounding_boxes.read( bounds, PredType::NATIVE_FLOAT);

        this->bounds_min = new Vector3[nblocks];
        this->bounds_max = new Vector3[nblocks];

        printf("Bounding Boxes for FLASH Blocks:\n");
        for(int i = 0; i < nblocks; i++) {
            bounds_min[i] = {0.01*(bounds[index(i,0,0)]-this->origin.x), 0.01*(bounds[index(i,1,0)]-this->origin.y), 0.01*(bounds[index(i,2,0)]-this->origin.z)};
            bounds_max[i] = {0.01*(bounds[index(i,0,1)]-this->origin.x), 0.01*(bounds[index(i,1,1)]-this->origin.y), 0.01*(bounds[index(i,2,1)]-this->origin.z)};

            printf("[%f, %f, %f] -> [%f, %f, %f]\n", bounds_min[i].x, bounds_min[i].y, bounds_min[i].z,
                                                     bounds_max[i].x, bounds_max[i].y, bounds_max[i].z);
        }

        dataspace.close();
        bounding_boxes.close();

        DataSet magx = file.openDataSet( "/magx" );

        /*
         * Get dataspace of grid variables.
         * TODO Assert on them all being equal?
         */
        dataspace = magx.getSpace();

        /*
         * Get the dimension size of each dimension in the dataspace and
         * display them.
         */
        hsize_t mag_dims[4];
        dataspace.getSimpleExtentDims( mag_dims, NULL);
        assert(nblocks == mag_dims[0]);
        int nxb = mag_dims[1], nyb = mag_dims[2], nzb = mag_dims[3];

        this->magx = new float[nblocks*nxb*nyb*nzb];
        this->magy = new float[nblocks*nxb*nyb*nzb];
        this->magz = new float[nblocks*nxb*nyb*nzb];

        this->nb[0] = nxb;
        this->nb[1] = nyb;
        this->nb[2] = nzb;

        printf("Magnetic fields are %dx%dx%d on %d blocks\n", nxb, nyb, nzb, nblocks);

        /*
         * Read data from the file into
         * memory and display the data.
         */
        magx.read( this->magx, PredType::NATIVE_FLOAT);
        magx.close();

        DataSet magy = file.openDataSet( "/magy" );
        magy.read( this->magy, PredType::NATIVE_FLOAT);
        magy.close();

        DataSet magz = file.openDataSet( "/magz" );
        magz.read( this->magz, PredType::NATIVE_FLOAT);
        magz.close();

        file.close();
    } catch(FileIException error) {
        error.printError();
    } catch(DataSetIException error) {
        error.printError();
    }

    blockIn = new int[this->N];
    cellIn = new int[this->N];

    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        this->E[j] = { 0.0, 0.0, 0.0 };
        this->B[j] = { 0.0, 0.0, 0.0 };
        this->blockIn[j] = 0;
        this->cellIn[j] = 0;
    }
}

bool FlashField::findBlock(Vector3 pos, int &block) {
    int i = 0;
    while(i < this->nblocks) {
        if(pos.x >= this->bounds_min[block].x && pos.x < this->bounds_max[block].x
        && pos.y >= this->bounds_min[block].y && pos.y < this->bounds_max[block].y
        && pos.z >= this->bounds_min[block].z && pos.z < this->bounds_max[block].z) {
            return true;
        } else {
            i++;
            block++;
            if(block >= nblocks) {
                block = 0;
            }
        }
    }
    return false;
}

void FlashField::getFields(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        if(state->running[j]) {
            bool foundBlock = findBlock(state->pos[j], this->blockIn[j]);

            if(foundBlock) {
                int block = this->blockIn[j];

                float sizex = (this->bounds_max[block].x - this->bounds_min[block].x);
                float sizey = (this->bounds_max[block].y - this->bounds_min[block].y); 
                float sizez = (this->bounds_max[block].z - this->bounds_min[block].z); 

                float dx = sizex / this->nb[0];
                float dy = sizey / this->nb[1];
                float dz = sizez / this->nb[2];

                float fracx = (state->pos[j].x - this->bounds_min[block].x) / dx;
                float fracy = (state->pos[j].y - this->bounds_min[block].y) / dy;
                float fracz = (state->pos[j].z - this->bounds_min[block].z) / dz;

                int idx_m = fracx > 0.5 ? floor(fracx-0.5) : this->nb[0] - 1;
                int idx_p = idx_m < this->nb[0] - 1 ? idx_m + 1 : 0;
                int idy_m = fracy > 0.5 ? floor(fracy-0.5) : this->nb[1] - 1;
                int idy_p = idy_m < this->nb[1] - 1 ? idy_m + 1 : 0; 
                int idz_m = fracz > 0.5 ? floor(fracz-0.5) : this->nb[2] - 1;
                int idz_p = idz_m < this->nb[2] - 1 ? idz_m + 1 : 0;
/*
                idx_m = floor(fracx);
                idx_p = floor(fracx);
                idy_m = floor(fracy);
                idy_p = floor(fracy);
                idz_m = floor(fracz);
                idz_p = floor(fracz);
*/
                int block_corner = 0;

                float magx_mmm;
                float magx_mmp;
                float magx_mpm;
                float magx_mpp;
                float magx_pmm;
                float magx_pmp;
                float magx_ppm;
                float magx_ppp;

                float magy_mmm;
                float magy_mmp;
                float magy_mpm;
                float magy_mpp;
                float magy_pmm;
                float magy_pmp;
                float magy_ppm;
                float magy_ppp;

                float magz_mmm;
                float magz_mmp;
                float magz_mpm;
                float magz_mpp;
                float magz_pmm;
                float magz_pmp;
                float magz_ppm;
                float magz_ppp;

                int index;
                //dx = 0; dy = 0; dz = 0;

                if(findBlock({state->pos[j].x - 0.5*dx, state->pos[j].y - 0.5*dy, state->pos[j].z - 0.5*dz}, block_corner)) {
                    index = ((block_corner * nb[2] + idz_m) * nb[1] + idy_m) * nb[0] + idx_m;
                    magx_mmm = this->magx[index];
                    magy_mmm = this->magy[index];
                    magz_mmm = this->magz[index];
                } else {
                    magx_mmm = 0;
                    magy_mmm = 0; 
                    magz_mmm = 0; 
                }

                if(findBlock({state->pos[j].x - 0.5*dx, state->pos[j].y - 0.5*dy, state->pos[j].z + 0.5*dz}, block_corner)) {
                    index = ((block_corner * nb[2] + idz_p) * nb[1] + idy_m) * nb[0] + idx_m;
                    magx_mmp = this->magx[index];
                    magy_mmp = this->magy[index];
                    magz_mmp = this->magz[index];
                } else {
                    magx_mmp = 0;
                    magy_mmp = 0; 
                    magz_mmp = 0; 
                }

                if(findBlock({state->pos[j].x - 0.5*dx, state->pos[j].y + 0.5*dy, state->pos[j].z - 0.5*dz}, block_corner)) {
                    index = ((block_corner * nb[2] + idz_m) * nb[1] + idy_p) * nb[0] + idx_m;
                    magx_mpm = this->magx[index];
                    magy_mpm = this->magy[index];
                    magz_mpm = this->magz[index];
                } else {
                    magx_mpm = 0;
                    magy_mpm = 0; 
                    magz_mpm = 0; 
                }

                if(findBlock({state->pos[j].x - 0.5*dx, state->pos[j].y + 0.5*dy, state->pos[j].z + 0.5*dz}, block_corner)) {
                    index = ((block_corner * nb[2] + idz_p) * nb[1] + idy_p) * nb[0] + idx_m;
                    magx_mpp = this->magx[index];
                    magy_mpp = this->magy[index];
                    magz_mpp = this->magz[index];
                } else {
                    magx_mpp = 0;
                    magy_mpp = 0; 
                    magz_mpp = 0; 
                }

                if(findBlock({state->pos[j].x + 0.5*dx, state->pos[j].y - 0.5*dy, state->pos[j].z - 0.5*dz}, block_corner)) {
                    index = ((block_corner * nb[2] + idz_m) * nb[1] + idy_m) * nb[0] + idx_p;
                    magx_pmm = this->magx[index];
                    magy_pmm = this->magy[index];
                    magz_pmm = this->magz[index];
                } else {
                    magx_pmm = 0;
                    magy_pmm = 0; 
                    magz_pmm = 0; 
                }

                if(findBlock({state->pos[j].x + 0.5*dx, state->pos[j].y - 0.5*dy, state->pos[j].z + 0.5*dz}, block_corner)) {
                    index = ((block_corner * nb[2] + idz_p) * nb[1] + idy_m) * nb[0] + idx_p;
                    magx_pmp = this->magx[index];
                    magy_pmp = this->magy[index];
                    magz_pmp = this->magz[index];
                } else {
                    magx_pmp = 0;
                    magy_pmp = 0; 
                    magz_pmp = 0; 
                }

                if(findBlock({state->pos[j].x + 0.5*dx, state->pos[j].y + 0.5*dy, state->pos[j].z - 0.5*dz}, block_corner)) {
                    index = ((block_corner * nb[2] + idz_m) * nb[1] + idy_p) * nb[0] + idx_p;
                    magx_ppm = this->magx[index];
                    magy_ppm = this->magy[index];
                    magz_ppm = this->magz[index];
                } else {
                    magx_ppm = 0;
                    magy_ppm = 0; 
                    magz_ppm = 0; 
                }

                if(findBlock({state->pos[j].x + 0.5*dx, state->pos[j].y + 0.5*dy, state->pos[j].z + 0.5*dz}, block_corner)) {
                    index = ((block_corner * nb[2] + idz_p) * nb[1] + idy_p) * nb[0] + idx_p;
                    magx_ppp = this->magx[index];
                    magy_ppp = this->magy[index];
                    magz_ppp = this->magz[index];
                } else {
                    magx_ppp = 0;
                    magy_ppp = 0; 
                    magz_ppp = 0; 
                }

                float del_x = (fracx-0.5) - floor(fracx-0.5);
                float del_y = (fracy-0.5) - floor(fracy-0.5);
                float del_z = (fracz-0.5) - floor(fracz-0.5);
                float idel_x = 1-del_x;
                float idel_y = 1-del_y;
                float idel_z = 1-del_z;

                float magx = idel_x * (idel_y * (idel_z * magx_mmm + del_z * magx_mmp) + del_y * (idel_z * magx_mpm + del_z * magx_mpp))
                            + del_x * (idel_y * (idel_z * magx_pmm + del_z * magx_pmp) + del_y * (idel_z * magx_ppm + del_z * magx_ppp));
                float magy = idel_x * (idel_y * (idel_z * magy_mmm + del_z * magy_mmp) + del_y * (idel_z * magy_mpm + del_z * magy_mpp))
                            + del_x * (idel_y * (idel_z * magy_pmm + del_z * magy_pmp) + del_y * (idel_z * magy_ppm + del_z * magy_ppp));
                float magz = idel_x * (idel_y * (idel_z * magz_mmm + del_z * magz_mmp) + del_y * (idel_z * magz_mpm + del_z * magz_mpp))
                            + del_x * (idel_y * (idel_z * magz_pmm + del_z * magz_pmp) + del_y * (idel_z * magz_ppm + del_z * magz_ppp));

                this->B[j] = {magx*1E-4, magy*1E-4, magz*1E-4};

            } else {
                this->B[j] = { 0.0, 0.0, 0.0 };
            }
        }
    }

}

