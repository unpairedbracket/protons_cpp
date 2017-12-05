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
            bounds_min[i] = {0.01*(bounds[index(i,0,0)]), 0.01*(bounds[index(i,1,0)]-0.05), 0.01*(bounds[index(i,2,0)]+0.05)};
            bounds_max[i] = {0.01*(bounds[index(i,0,1)]), 0.01*(bounds[index(i,1,1)]-0.05), 0.01*(bounds[index(i,2,1)]+0.05)};

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

void FlashField::getFields(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        if(state->running[j]) {
            bool foundBlock = false;
            int block;
            int i = 0;
            while(!foundBlock && i < this->nblocks - 1) {
                if(state->pos[j].x >= this->bounds_min[this->blockIn[j]].x && state->pos[j].x < this->bounds_max[this->blockIn[j]].x
                && state->pos[j].y >= this->bounds_min[this->blockIn[j]].y && state->pos[j].y < this->bounds_max[this->blockIn[j]].y
                && state->pos[j].z >= this->bounds_min[this->blockIn[j]].z && state->pos[j].z < this->bounds_max[this->blockIn[j]].z) {
                    foundBlock = true;
                    block = this->blockIn[j];
                } else {
                    i++;
                    (this->blockIn[j])++;
                    if(this->blockIn[j] >= nblocks) {
                        this->blockIn[j] = 0;
                    }
                }
            }

            if(foundBlock) {
                float fracx = (state->pos[j].x - this->bounds_min[block].x) / (this->bounds_max[block].x - this->bounds_min[block].x);
                float fracy = (state->pos[j].y - this->bounds_min[block].y) / (this->bounds_max[block].y - this->bounds_min[block].y);
                float fracz = (state->pos[j].z - this->bounds_min[block].z) / (this->bounds_max[block].z - this->bounds_min[block].z);

                int idx = floor(fracx * this->nb[0]);
                int idy = floor(fracy * this->nb[1]);
                int idz = floor(fracz * this->nb[2]);

                int index = ((block * nb[0] + idx) * nb[1] + idy) * nb[2] + idz;

                this->B[j] = {this->magx[index]*1E-4, this->magy[index]*1E-4, this->magz[index]*1E-4};

            } else {
                this->B[j] = { 0.0, 0.0, 0.0 };
            }
        }
    }

}

