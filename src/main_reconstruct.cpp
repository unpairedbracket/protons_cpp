#include "main.h"

#include <cmath>

#include <H5Cpp.h>

#include "reconstruction/invert.h"
#include "util/ndarray.h"

using namespace H5;

typedef ArrayND<2, double> Array;

int main(int argc, char *argv[]) {
    int pixels[2];
    double detectorSize[2];

    // Note fluence and fluence_undeviated are arrays, fluence_expected is a scalar value
    Array fluence, fluence_undeviated;
    double fluence_expected;
    Array Phi, X, Y, XX, YY, XY;

    H5File file(argv[1], H5F_ACC_RDWR);

    DataSet fluence_data = file.openDataSet( "/fluence" );
    DataSpace dataspace = fluence_data.getSpace();

    hsize_t dims[2];
    dataspace.getSimpleExtentDims( dims, NULL);
    pixels[0] = dims[0]; pixels[1] = dims[1];
    int num_pixels = pixels[0]*pixels[1];

    fluence.assignMemory(pixels);
    fluence_data.read(fluence.memory, PredType::NATIVE_DOUBLE);
    dataspace.close();
    fluence_data.close();

    DataSet undeviated_data = file.openDataSet( "/fluence_0" );
    fluence_undeviated.assignMemory(pixels);

    undeviated_data.read(fluence_undeviated.memory, PredType::NATIVE_DOUBLE);
    undeviated_data.close();

    DataSet expected_data = file.openDataSet( "/fluence_expected" );
    expected_data.read(&fluence_expected, PredType::NATIVE_DOUBLE);
    expected_data.close();

    DataSet size_data = file.openDataSet( "/detector_size" );
    size_data.read(detectorSize, PredType::NATIVE_DOUBLE);
    size_data.close();

    printf("[%d, %d] pixels, size %e by %e metres\n", pixels[0], pixels[1], detectorSize[0], detectorSize[1]);

    Phi.assignMemory(pixels);
    X.assignMemory(pixels);
    Y.assignMemory(pixels);

    XX.assignMemory(pixels);
    YY.assignMemory(pixels);
    XY.assignMemory(pixels);

    for(int i = 0; i < num_pixels; i++) {
        // Fluence that's actually zero will blow up to -inf in the logarithm - avoid
        fluence[i] = std::max(fluence_expected / 1000, fluence[i] - fluence_undeviated[i] + fluence_expected);
    }

    double dt = 1e-10;
    if(argc > 2) {
        dt = atof(argv[2]);
        if(dt == 0) { dt = 1e-10; }
    }

    int N = 10000;
    if(argc > 3) {
        N = atoi(argv[3]);
    }

    initialise_inversion_arrays(pixels, detectorSize, Phi, X, Y, XX, XY, YY);

    while(N > 0) {
        DataSet dataset = file.openDataSet("/potential");
        dataset.write(Phi.memory, PredType::NATIVE_DOUBLE);
        dataset.close();

        dataset = file.openDataSet("/X");
        dataset.write(X.memory, PredType::NATIVE_DOUBLE);
        dataset.close();

        dataset = file.openDataSet("/Y");
        dataset.write(Y.memory, PredType::NATIVE_DOUBLE);
        dataset.close();

        dataset = file.openDataSet("/XX");
        dataset.write(XX.memory, PredType::NATIVE_DOUBLE);
        dataset.close();

        dataset = file.openDataSet("/XY");
        dataset.write(XY.memory, PredType::NATIVE_DOUBLE);
        dataset.close();

        dataset = file.openDataSet("/YY");
        dataset.write(YY.memory, PredType::NATIVE_DOUBLE);
        dataset.close();

        file.flush(H5F_SCOPE_LOCAL);

        int Nrun = std::min(N, 100);

        invert_fluences(fluence, fluence_expected, pixels, detectorSize, dt, Nrun, Phi, X, Y, XX, XY, YY);

        N -= Nrun;
    }
    file.close();
}

