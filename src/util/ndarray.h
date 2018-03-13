#pragma once

#include <cstdarg>

template<int N, typename T>
struct ArrayND {
    T* memory;
    int dimensions[N];

    inline int getIndex(int inds[N] ) {
        int index = 0;
        for(int i = 0; i < N; i++) {
            index = index * this->dimensions[i] + inds[i];
        }
        return index;
    }

    inline int getIndex( int idx, va_list indices ) {
        int index = idx;
        for(int i = 1; i < N; i++) {
            int idxi = va_arg(indices, int);
            index = index * this->dimensions[i] + idxi;
        }
        return index;
    }

    inline int getIndex( int idx, ... ) {
        va_list indices;
        va_start(indices, idx);
        int index = getIndex(idx, indices);
        va_end(indices);
        return index;
    }

    void assignMemory(int dims[N]) {
        int numel = 1;
        for(int i = 0; i < N; i++) {
            int Ni = dims[i];
            this->dimensions[i] = Ni;
            numel *= Ni;
        }
        this->memory = new T[numel];
    }

    void assignMemory( int Nx, ... ) {
        va_list N_cells;
        va_start(N_cells, Nx);
        int numel = Nx;
        for(int i = 1; i < N; i++) {
            int Ni = va_arg(N_cells, int);
            this->dimensions[i] = Ni;
            numel *= Ni;
        }
        this->memory = new T[numel];
        va_end(N_cells);
    }

    void freeMemory() {
        delete[] this->memory;
    }

    T& operator()( int inds[N] ) {
        return this->memory[getIndex( inds )];
    }

    T& operator()( int idx, ... ) {
        va_list indices;
        va_start(indices, idx);
        T& result = this->memory[getIndex( idx, indices )];
        va_end(indices);
        return result;
    }

    T& operator[](int linear_index) {
        return this->memory[linear_index];
    }
};

template<typename T>
struct ArrayND<2, T> {
    T* memory;
    int dimensions[2];

    inline int getIndex(int xidx, int yidx) {
        return xidx * this->dimensions[1] + yidx;
    }

    void assignMemory(int dims[2]) {
        this->dimensions[0] = dims[0];
        this->dimensions[1] = dims[1];
        this->memory = new T[dims[0] * dims[1]];
    }

    void assignMemory(int Nx, int Ny) {
        this->dimensions = { Nx, Ny };
        this->memory = new T[Nx * Ny];
    }

    void freeMemory() {
        delete[] this->memory;
    }

    T& operator()(int xidx, int yidx) {
        return this->memory[getIndex(xidx, yidx)];
    }

    T& operator[](int linear_index) {
        return this->memory[linear_index];
    }
};

