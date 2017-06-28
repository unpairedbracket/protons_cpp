#pragma once

struct FieldStructure {
    double *Ex, *Ey, *Ez;
    double *Bx, *By, *Bz;
};

FieldStructure* makeFieldStructure(long N);

void getFields(FieldStructure* fields, double* x, double* y, double* z, long N);
