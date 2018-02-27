#pragma once

void initialise_inversion_arrays(int pixels[2], double size[2], double* potential_out, double* X_out, double* Y_out,
        double* XX_out, double* XY_out, double* YY_out);
void invert_fluences(double* fluence_adjusted, double fluence_expected_per_area, int pixels[2], double size[2], double dt, int N,
        double* potential_out, double* X_out, double* Y_out, double* XX_out, double* XY_out, double* YY_out);
