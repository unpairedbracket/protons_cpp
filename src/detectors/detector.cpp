#include "detector.h"

#include <cmath>
#include <iostream>
#include <fstream>

#include <H5Cpp.h>

void ParticleDetector::init(ParticleInfo* particle, double distance) {
    this->particleInfo = particle;
    this->distance = distance;
}

void ParticleDetector::finalPush(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        state->pos[j].x += (distance - state->pos[j].z) * state->vel[j].x / state->vel[j].z;
        state->pos[j].y += (distance - state->pos[j].z) * state->vel[j].y / state->vel[j].z;
        state->pos[j].z += (distance - state->pos[j].z);
    }
}

void DetectorTextFile::detectUndeviated(ParticleState* state) {
    //This will be overridden by detect() but I don't care about that at this point
    this->state = state;
}

void DetectorTextFile::detect(ParticleState* state) {
    this->state = state;
}

void DetectorTextFile::output() {
    std::ofstream posfile;
    posfile.open("pos.txt");
    for(long i = 0; i < state->N; i++) {
        posfile << state->pos[i].x << "," << state->pos[i].y << "," << state->pos[i].z << std::endl;
    }
    posfile.close();

    std::ofstream velfile;
    velfile.open("vel.txt");
    for(long i = 0; i < state->N; i++) {
        velfile << state->vel[i].x << "," << state->vel[i].y << "," << state->vel[i].z << std::endl;
    }
    velfile.close();
}

void DetectorHDF5::detectUndeviated(ParticleState* state) {
    //This will be overridden by detect() but I don't care about that at this point
    this->state = state;
}

void DetectorHDF5::detect(ParticleState* state) {
    this->state = state;
}

void DetectorHDF5::output() {
    using namespace H5;
    try {
        Exception::dontPrint();

        H5File file("output.h5", H5F_ACC_TRUNC);

        hsize_t dims[2];               // dataset dimensions
        dims[0] = state->N;
        dims[1] = 3;
        DataSpace dataspace(2, dims);

        // Create the dataset.      
        DataSet dataset_pos = file.createDataSet("Positions",  PredType::IEEE_F64BE, dataspace);
        DataSet dataset_vel = file.createDataSet("Velocities",  PredType::IEEE_F64BE, dataspace);

        // Write the data to the dataset using default memory space, file
        // space, and transfer properties.
        dataset_pos.write(state->pos, PredType::NATIVE_DOUBLE);
        dataset_vel.write(state->vel, PredType::NATIVE_DOUBLE);

        dataspace.close();
        dataset_pos.close();
        dataset_vel.close();
        file.close();
    } catch(FileIException error) {
    	error.printError();
    } catch(DataSetIException error) {
        error.printError();
    }
}

void DetectorFluence::detect(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        double xfrac = state->pos[j].x / this->detectorSize[0] + 0.5;
        double yfrac = state->pos[j].y / this->detectorSize[1] + 0.5;
        if(xfrac >= 0 && xfrac <= 1 && yfrac >= 0 && yfrac <= 1) {
            int xcell = floor(this->detectorPixels[0] * xfrac);
            int ycell = floor(this->detectorPixels[1] * yfrac);
            int idx = ycell * this->detectorPixels[0] + xcell;
            #pragma omp atomic
            detectorArray[idx]++;
        }
    }
}

void DetectorFluence::detectUndeviated(ParticleState* state) {
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        double d = this->distance - state->pos[j].z;
        double x = state->pos[j].x + d * state->vel[j].x / state->vel[j].z;
        double y = state->pos[j].y + d * state->vel[j].y / state->vel[j].z;
        double xfrac = x / this->detectorSize[0] + 0.5;
        double yfrac = y / this->detectorSize[1] + 0.5;
        if(xfrac >= 0 && xfrac <= 1 && yfrac >= 0 && yfrac <= 1) {
            int xcell = floor(this->detectorPixels[0] * xfrac);
            int ycell = floor(this->detectorPixels[1] * yfrac);
            int idx = ycell * this->detectorPixels[0] + xcell;
            #pragma omp atomic
            nullDetectorArray[idx]++;
        }
    }
}

void DetectorFluence::output() {
    using namespace H5;
    try {
        Exception::dontPrint();

        H5File file("fluence.h5", H5F_ACC_TRUNC);

        hsize_t dims[2];               // dataset dimensions
        dims[0] = this->detectorPixels[1];
        dims[1] = this->detectorPixels[0];
        DataSpace dataspace(2, dims);

        // Create the dataset.
        DataSet dataset = file.createDataSet("fluence",  PredType::IEEE_F64BE, dataspace);
        dataset.write(this->detectorArray, PredType::NATIVE_DOUBLE);
        dataset.close();

        dataset = file.createDataSet("fluence_0",  PredType::IEEE_F64BE, dataspace);
        dataset.write(this->nullDetectorArray, PredType::NATIVE_DOUBLE);

        dataspace.close();
        dataset.close();
        file.close();
    } catch(FileIException error) {
        error.printError();
    } catch(DataSetIException error) {
        error.printError();
    }
}
