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
        int nblocks = dims_out[0];

        this->bounds = new float[nblocks*3*2];
        auto index = [](int block, int dir, int minmax) {return 6*block + 2 * dir + minmax;};

        /*
         * Read data from hyperslab in the file into the hyperslab in
         * memory and display the data.
         */
        bounding_boxes.read( this->bounds, PredType::NATIVE_FLOAT);

        printf("Bounding Boxes for FLASH Blocks:\n");
        for(int i = 0; i < nblocks; i++) {
            printf("[%f, %f, %f] -> [%f, %f, %f]\n", bounds[index(i,0,0)], bounds[index(i,1,0)], bounds[index(i,2,0)],
                                                     bounds[index(i,0,1)], bounds[index(i,1,1)], bounds[index(i,2,1)]);
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
        int NDIM = dataspace.getSimpleExtentDims( mag_dims, NULL);
        assert(nblocks == mag_dims[0]);
        int nxb = mag_dims[1], nyb = mag_dims[2], nzb = mag_dims[3];

        this->magx = new float[nblocks*nxb*nyb*nzb];
        this->magy = new float[nblocks*nxb*nyb*nzb];
        this->magz = new float[nblocks*nxb*nyb*nzb];

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
    #pragma omp parallel for
    for(long j = 0; j < this->N; j++) {
        this->E[j] = { 0.0, 0.0, 0.0 };
        this->B[j] = { 0.0, 0.0, 0.0 };
    }
}

void FlashField::getFields(ParticleState* state) {
}

