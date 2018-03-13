#pragma once

#include <cassert>

#include <H5Cpp.h>

using namespace H5;

template<typename T>
DataType type = PredType::NATIVE_CHAR;

/*
 * Templates for getting the DataType we want for a given dataset
 * There should be one for all of them but I will just do them as I need them
 */

template<>
DataType type<double> = PredType::NATIVE_DOUBLE;

template<>
DataType type<float> = PredType::NATIVE_FLOAT;

template<int N, typename T>
ArrayND<N, T> createArrayFromDataSet(DataSet data) {
    DataSpace dataspace = data.getSpace();

    /*
     * Get the dimension size of each dimension in the dataspace and
     * display them.
     */
    hsize_t N_actual = dataspace.getSimpleExtentNdims();
    assert(N == N_actual);
    hsize_t dims_hsize[N];
    int dims[N];

    dataspace.getSimpleExtentDims( dims_hsize, NULL);

    for(int i = 0; i<N; i++) {
        dims[i] = dims_hsize[i];
    }

    /*
     * Read data from the file into
     * memory and display the data.
     */
    ArrayND<N, T> result;
    result.assignMemory(dims);
    data.read( result.memory, type<T>);
    return result;
}

DataSet createOrOpenDataSet(H5Location &location, const char* name, DataType type, DataSpace fileSpace) {
    //Exception::dontPrint();
    if(location.exists(name)) {
        return location.openDataSet(name);
    } else {
        return location.createDataSet(name, type, fileSpace);
    }
}

Group createOrOpenGroup(H5File &location, const char* name) {
    //Exception::dontPrint();
    if(location.exists(name)) {
        return location.openGroup(name);
    } else {
        return location.createGroup(name);
    }
}
