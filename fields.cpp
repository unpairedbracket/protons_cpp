#include "fields.h"

FieldStructure* makeFieldStructure(long N) {
    FieldStructure* fields = new FieldStructure();
    fields->Ex = new double[N];
    fields->Ey = new double[N];
    fields->Ez = new double[N];

    fields->Bx = new double[N];
    fields->By = new double[N];
    fields->Bz = new double[N];

    return fields;
}

void getFields(FieldStructure* fields, double* x, double* y, double* z, long N) {
    #pragma omp parallel for
    for(long j = 0; j < N; j++) {
        fields->Ex[j] = 0.0;
        fields->Ey[j] = 0.0;
        fields->Ez[j] = +0.1;

        fields->Bx[j] = 0.0;
        fields->By[j] = 0.0;
        fields->Bz[j] = 0.0;
    }
}

