#pragma once

void invert_fluences(double* fluence_adjusted, double fluence_expected_per_area, int pixels[2], double size[2], double dt,
        double* potential_out, double* X_out, double* Y_out, double* XX_out, double* XY_out, double* YY_out);
