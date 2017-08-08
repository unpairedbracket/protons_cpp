#include "fields.h"

void initFieldArrays(FieldStructure* field, long N) {
    field->N = N;

    field->E = new Vector3[N];
    field->B = new Vector3[N];
}
