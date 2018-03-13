#pragma once

#include "../util/ndarray.h"

typedef ArrayND<2, double> Array;

void initialise_inversion_arrays(int pixels[2], double size[2], Array potential_out, Array X_out, Array Y_out,
        Array XX_out, Array XY_out, Array YY_out);
void invert_fluences(Array fluence_adjusted, double fluence_expected_per_area, int pixels[2], double size[2], double dt, int N,
        Array potential_out, Array X_out, Array Y_out, Array XX_out, Array XY_out, Array YY_out);
