#include "invert.h"

#include <cmath>
#include <cstdio>

double getCoord(int index, double size, double steps) {
    double delta = size / steps;
    double xmin = - (size/2);
    return xmin + ((double)index + 0.5) * delta;
}

void initialise_inversion_arrays(int pixels[2], double size[2], Array potential_out,
        Array X_out, Array Y_out, Array XX_out, Array XY_out, Array YY_out) {

    int n_pixels = pixels[0] * pixels[1];

    // Setup initial potential and gradients
    for(int index = 0; index < n_pixels; index++) {
        int ix = index / pixels[1];
        int iy = index % pixels[1];

        double x = getCoord(ix, size[0], pixels[0]);
        double y = getCoord(iy, size[1], pixels[1]);
        potential_out[index] = 0.5 * ( x*x + y*y );
        X_out[index] = x;
        Y_out[index] = y;

        XX_out[index] = 1;
        XY_out[index] = 0;
        YY_out[index] = 1;
    }

}

void invert_fluences(Array fluence_adjusted, double fluence_expected, int pixels[2], double size[2], double dt, int N,
        Array potential_out, Array X_out, Array Y_out,
        Array XX_out, Array XY_out, Array YY_out) {
    using std::abs;

    int n_pixels = pixels[0] * pixels[1];
    double dx = size[0] / pixels[0];
    double dy = size[1] / pixels[1];

    //Out-parameters
    Array potential = potential_out;
    Array dpdx = X_out;
    Array dpdy = Y_out;

    //Internal things
    Array d2pdxdx = XX_out;
    Array d2pdxdy = XY_out;
    Array d2pdydy = YY_out;

    auto interp = [&fluence_adjusted, size, pixels, dx, dy, fluence_expected] (double x, double y) -> double {
        if(abs(x) >= size[0]/2 || abs(y) >= size[1]/2) {
            return fluence_expected;
        }

        double fracx = (x + size[0]/2) / dx;
        double fracy = (y + size[1]/2) / dy;

        int idx_m = floor(fracx-0.5);
        int idy_m = floor(fracy-0.5);

        double del_x = (fracx-0.5) - idx_m;
        double del_y = (fracy-0.5) - idy_m;

        double mm, mp, pm, pp;

        if(idx_m < 0 || idy_m < 0) {
            mm = fluence_expected;
        } else {
            mm = fluence_adjusted(idx_m, idy_m);
        }

        if(idx_m < 0 || idy_m + 1 > pixels[1] - 1) {
            mp = fluence_expected;
        } else {
            mp = fluence_adjusted(idx_m, idy_m+1);
        }

        if(idx_m + 1 > pixels[0] - 1 || idy_m < 0) {
            pm = fluence_expected;
        } else {
            pm = fluence_adjusted(idx_m+1, idy_m);
        }

        if(idx_m + 1 > pixels[0] - 1 || idy_m + 1 > pixels[1] - 1) {
            pp = fluence_expected;
        } else {
            pp = fluence_adjusted(idx_m+1, idy_m+1);
        }

        double fluence_my;
        double fluence_py;

        fluence_my = mm + (pm - mm) * del_x;
        fluence_py = mp + (pp - mp) * del_x;

        return fluence_my + (fluence_py - fluence_my) * del_y;
    };

    for(int i = 0; i < N; i++) {
        double sumsq = 0;

        for(int index = 0; index < n_pixels; index++) {
            double det = d2pdxdx[index] * d2pdydy[index] - d2pdxdy[index] * d2pdxdy[index];

            double fluence = interp(dpdx[index], dpdy[index]);

            double Fdt = log(fluence*abs(det)/fluence_expected) * dt;

            //Reduction
            sumsq += Fdt*Fdt / (potential[index] * potential[index]);

            potential[index] += Fdt;
        }

        printf("%d: rms proportional diff in potential: %e\n", i, sqrt(sumsq / (pixels[0] * pixels[1])));

        for(int index = 0; index < n_pixels; index++) {
            int ix = index / pixels[1];
            int iy = index % pixels[1];

            if(ix == 0 || ix == pixels[0] - 1) {
                dpdx[index] = getCoord(ix, size[0], pixels[0]);
            } else {
                dpdx[index] = (potential(ix+1, iy) - potential(ix-1, iy)) / (2 * dx);
            }

            if(iy == 0 || iy == pixels[1] - 1) {
                dpdy[index] = getCoord(iy, size[1], pixels[1]);
            } else {
                dpdy[index] = (potential(ix, iy+1) - potential(ix, iy-1)) / (2 * dy);
            }

            int im = ix-1;
            int i0 = ix;
            int ip = ix+1;

            if(ix == 0) {
                im++;
                i0++;
                ip++;
            } else if(ix == pixels[0] - 1) {
                im--;
                i0--;
                ip--;
            }
            d2pdxdx[index] = (potential(ip, iy) - 2 * potential(i0, iy) + potential(im, iy)) / (dx * dx);

            im = iy-1;
            i0 = iy;
            ip = iy+1;

            if(iy == 0) {
                im++;
                i0++;
                ip++;
            } else if(iy == pixels[1] - 1) {
                im--;
                i0--;
                ip--;
            }
            d2pdydy[index] = (potential(ix, ip) - 2 * potential(ix, i0) + potential(ix, im)) / (dy * dy);

            if(ix == 0 || ix == pixels[0] - 1 || iy == 0 || iy == pixels[1] - 1) {
                d2pdxdy[index] = 0;
            } else {
                double pp = potential(ix+1, iy+1);
                double pm = potential(ix+1, iy-1);
                double mp = potential(ix-1, iy+1);
                double mm = potential(ix-1, iy-1);
                d2pdxdy[index] = ((pp - mp)
                                - (pm - mm)) / (4 * dx *dy);
            }
        }
    }
}
