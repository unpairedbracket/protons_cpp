#include "main.h"

#include <cstdio>

#include <chrono>
#include <fstream>
#include <iostream>

#include <omp.h>

#include <H5Cpp.h>

#include "reconstruction/invert.h"

int main(int argc, char *argv[]) {
    int pixels[2];
    int num_pixels;
    double detectorSize[2];

    // Note fluence and fluence_undeviated are arrays, fluence_expected is a scalar value
    double *fluence, *fluence_undeviated, fluence_expected;
    double *Phi, *X, *Y, *XX, *YY, *XY;

    using namespace H5;
    {
        H5File file(argv[1], H5F_ACC_RDONLY);

        DataSet fluence_data = file.openDataSet( "/fluence" );
        DataSpace dataspace = fluence_data.getSpace();

        hsize_t dims[2];
        dataspace.getSimpleExtentDims( dims, NULL);
        pixels[0] = dims[1]; pixels[1] = dims[0];
        num_pixels = pixels[0]*pixels[1];

        fluence = new double[num_pixels];
        fluence_data.read(fluence, PredType::NATIVE_DOUBLE);
        dataspace.close();
        fluence_data.close();

        DataSet undeviated_data = file.openDataSet( "/fluence_0" );
        fluence_undeviated = new double[num_pixels];
        undeviated_data.read(fluence_undeviated, PredType::NATIVE_DOUBLE);
        undeviated_data.close();

        DataSet expected_data = file.openDataSet( "/fluence_expected" );
        expected_data.read(&fluence_expected, PredType::NATIVE_DOUBLE);
        expected_data.close();

        DataSet size_data = file.openDataSet( "/detector_size" );
        size_data.read(&detectorSize, PredType::NATIVE_DOUBLE);
        size_data.close();

        file.close();
    }

    Phi = new double[num_pixels];
    X = new double[num_pixels];
    Y = new double[num_pixels];

    XX = new double[num_pixels];
    YY = new double[num_pixels];
    XY = new double[num_pixels];

    for(int i = 0; i < num_pixels; i++) {
        // Fluence that's actually zero will blow up to -inf in the logarithm - avoid
        fluence[i] = std::max(fluence_expected / 1000, fluence[i] - fluence_undeviated[i] + fluence_expected);
    }

    double dt = 1e-10;
    if(argc > 2) {
        dt = atof(argv[2]);
        if(dt == 0) { dt = 1e-10; }
    }

    invert_fluences(fluence, fluence_expected, pixels, detectorSize, dt, Phi, X, Y, XX, XY, YY);


    using namespace H5;
    {
        H5File file(argv[1], H5F_ACC_RDWR);

        hsize_t dims[2];
        dims[0] = pixels[1];
        dims[1] = pixels[0];

        DataSet dataset = file.openDataSet("/potential");
        dataset.write(Phi, PredType::NATIVE_DOUBLE);
        dataset.close();

        dataset = file.openDataSet("/X");
        dataset.write(X, PredType::NATIVE_DOUBLE);
        dataset.close();

        dataset = file.openDataSet("/Y");
        dataset.write(Y, PredType::NATIVE_DOUBLE);
        dataset.close();

        dataset = file.openDataSet("/XX");
        dataset.write(XX, PredType::NATIVE_DOUBLE);
        dataset.close();

        dataset = file.openDataSet("/XY");
        dataset.write(XY, PredType::NATIVE_DOUBLE);
        dataset.close();

        dataset = file.openDataSet("/YY");
        dataset.write(YY, PredType::NATIVE_DOUBLE);
        dataset.close();

        file.close();
    }
}

