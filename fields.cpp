#include "fields.h"

FieldStructure* makeFieldStructure(long N) {
    FieldStructure* fields = new FieldStructure();
    fields->N = N;
    fields->Ex = new double[N];
    fields->Ey = new double[N];
    fields->Ez = new double[N];

    fields->Bx = new double[N];
    fields->By = new double[N];
    fields->Bz = new double[N];

    return fields;
}

