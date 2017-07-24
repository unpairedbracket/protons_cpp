#include "integrator.h"

const double integrator::a2[] = {1.0/5.0};
const double integrator::a3[] = {3.0/40.0, 9.0/40.0};
const double integrator::a4[] = {44.0/45.0, -56.0/15.0, 32.0/9.0};
const double integrator::a5[] = {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0};
const double integrator::a6[] = {9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0};
const double integrator::a7[] = {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0};

const double integrator::b4[] = {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0};
const double integrator::b5[] = {5179.0/57600.0, 0.0, 7571.0/16695.0, 393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 1.0/40.0};

integrator* makeIntegrator(long N) {
    integrator* method = new integrator();
    method->accelx = new double[N];
    method->accely = new double[N];
    method->accelz = new double[N];
    method->dt = new double[N];
    return method;
}

void step(
        integrator* method,
        ParticleState* state, FieldStructure* field,
        void (*accelFunc)(ParticleState*, FieldStructure*, double*, double*, double*)
    ) {

    accelFunc(state, field, method->accelx, method->accely, method->accelz);
    
    #pragma omp parallel for
    for(long j = 0; j < state->N; j++) {
        if(state->running[j])
        {
            state->posX[j] += state->velX[j] * method->dt[j];
            state->posY[j] += state->velY[j] * method->dt[j];
            state->posZ[j] += state->velZ[j] * method->dt[j];

            state->velX[j] += method->accelx[j] * method->dt[j];
            state->velY[j] += method->accely[j] * method->dt[j];
            state->velZ[j] += method->accelz[j] * method->dt[j];
        }
    }
}
//function [posOut, velOut, accelOut, dtOut] = RKDPIntegrate(posIn, velIn, accelIn, dt, func)
//    global verbose
    //%% Integration params
    //% Coefficients for calculating the various midpoints
    //a = [
           //0             0             0             0             0             0             0
           //1/5           0             0             0             0             0             0
           //3/40          9/40          0             0             0             0             0
           //44/45        -56/15         32/9          0             0             0             0
           //19372/6561   -25360/2187    64448/6561   -212/729       0             0             0
           //9017/3168    -355/33        46732/5247    49/176       -5103/18656    0             0
           //35/384        0             500/1113      125/192      -2187/6784     11/84         0
    //];
//
    //% Coefficients for 4th-order method
    //b4 = a(7, :);
    //% Coefficients for 5th-order method
    //b5 = [ 5179/57600    0             7571/16695    393/640      -92097/339200  187/2100      1/40 ];
    //% Coefficients for error between 4th- and 5th-order
    //E = b4 - b5;
    //
    //%% Adaptation params
    //% How large a proportion of the velocity are we allowed to be wrong by?
    //rtol = 1e-3;
    //% err ~ O(h^5) so we take our step size down by (err/tol)^(1/5) if we're wrong
    //pow = 1/5;
    //% How much are we allowed to increase step size by per non-failed step?
    //maxLengthen = 1.1;
    //% How much are we allowed to decrease step size by on first failure?
    //maxFirstShorten = 10;
    //% How much do we decrease step size by on subsequent failures?
    //maxOtherShorten = 2;
    //% How short is a step allowed to be
    //dt_min = 5E-16; % seconds
 //
    //
    //% Preallocate output vars
    //posOut = zeros(size(posIn));
    //velOut = zeros(size(velIn));
    //accelOut = zeros(size(accelIn));
    //dtOut = dt;
//
    //% Assign the first set of inputs
    //pos1 = posIn;
    //vel1 = velIn;
    //accel1 = accelIn;
    //
    //% All of our input protons are running to begin with
    //idxrunning = true(size(pos1(:, 1)));
    //
    //nofailures = true;
    //
    //while size(pos1, 1) > 0
    //
        //vel2 = vel1 + dt .* ( a(2,1) * accel1 ); % v2 = v1 + h a21 kv1
        //pos2 = pos1 + dt .* ( a(2,1) *   vel1 ); % x2 = x1 + h a21 kx1
//
        //accel2 = func(pos2, vel2); % kv2 = f(x2, v2)
//
        //vel3 = vel1 + dt .* ( a(3,1) * accel1 + a(3,2) * accel2 ); % v3 = v1 + h ( a31 kv1 + a32 kv2 )
        //pos3 = pos1 + dt .* ( a(3,1) *   vel1 + a(3,2) *   vel2 ); % x3 = x1 + h ( a31 kx1 + a32 kx2 )
//
        //accel3 = func(pos3, vel3); % kv3 = f(x3, v3)
//
        //vel4 = vel1 + dt .* ( a(4,1) * accel1 + a(4,2) * accel2 + a(4,3) * accel3 ); % v4 = v1 + h ( a41 kv1 + a42 kv2 + a43 kv3 )
        //pos4 = pos1 + dt .* ( a(4,1) *   vel1 + a(4,2) *   vel2 + a(4,3) *   vel3 ); % x4 = x1 + h ( a41 kx1 + a42 kx2 + a43 kx3 )
//
        //accel4 = func(pos4, vel4); % kv4 = f(x4, v4)
//
        //vel5 = vel1 + dt .* ( a(5,1) * accel1 + a(5,2) * accel2 + a(5,3) * accel3 + a(5,4) * accel4 ); % v5 = v1 + h ( a51 kv1 + a52 kv2 + a53 kv3 + a54 kv4 )
        //pos5 = pos1 + dt .* ( a(5,1) *   vel1 + a(5,2) *   vel2 + a(5,3) *   vel3 + a(5,4) *   vel4 ); % x5 = x1 + h ( a51 kx1 + a52 kx2 + a53 kx3 + a54 kx4 )
//
        //accel5 = func(pos5, vel5); % kv5 = f(x5, v5)
//
        //vel6 = vel1 + dt .* ( a(6,1) * accel1 + a(6,2) * accel2 + a(6,3) * accel3 + a(6,4) * accel4 + a(6,5) * accel5 ); % v6 = v1 + h ( a61 kv1 + a62 kv2 + a63 kv3 + a64 kv4 + a65 kv5 )
        //pos6 = pos1 + dt .* ( a(6,1) *   vel1 + a(6,2) *   vel2 + a(6,3) *   vel3 + a(6,4) *   vel4 + a(6,5) *   vel5 ); % x6 = x1 + h ( a61 kx1 + a62 kx2 + a63 kx3 + a64 kx4 + a65 kx5 )
//
        //accel6 = func(pos6, vel6); % kv6 = f(x6, v6)
//
        //vel7 = vel1 + dt .* ( a(7,1) * accel1 + a(7,2) * accel2 + a(7,3) * accel3 + a(7,4) * accel4 + a(7,5) * accel5 + a(7,6) * accel6 ); % v7 = v1 + h ( a71 kv1 + a72 kv2 + a73 kv3 + a74 kv4 + a75 kv5 + a76 kv6 )
        //pos7 = pos1 + dt .* ( a(7,1) *   vel1 + a(7,2) *   vel2 + a(7,3) *   vel3 + a(7,4) *   vel4 + a(7,5) *   vel5 + a(7,6) *   vel6 ); % x7 = x1 + h ( a71 kx1 + a72 kx2 + a73 kx3 + a74 kx4 + a75 kx5 + a76 kx6 )
//
        //accel7 = func(pos7, vel7); % kv7 = f(x7, v7)
//
        //% Absolute difference in acceleration between 4th- and 5th-order.
        //dErr = E(1) * accel1 + E(2) * accel2 + E(3) * accel3 + E(4) * accel4 + E(5) * accel5 + E(6) * accel6 + E(7) * accel7;
        //
        //% Proportion error in |vel| over the timestep, per proton.
        //err = dt .* sqrt(dot(dErr, dErr, 2) ./ dot(vel1, vel1, 2));
        //if verbose
            //display(sprintf('        Errors: min: %.1e, max: %.1e, mean: %.1e', min(err), max(err), mean(err)));
        //end
//
        //idxfail = err > rtol; % Indices of failed particles, which we'll want to run again
        //
        //if nofailures
            //temp = 1.25 * (err / rtol).^pow;
            //temp(~idxfail & temp < 1/maxLengthen) = 1/maxLengthen;
            //temp( idxfail & temp > maxFirstShorten ) = maxFirstShorten;
            //dt = dt./temp;
        //else
            //dt(idxfail) = dt(idxfail) / maxOtherShorten;
        //end
        //
        //if sum(idxfail) == 0 % All our particles are fine! Yay!
            //velOut(idxrunning, :) = vel7;
            //posOut(idxrunning, :) = pos7;
            //accelOut(idxrunning, :) = accel7;
            //dtOut(idxrunning) = dt;
            //break;
        //end
        //
        //if verbose
            //display(sprintf('            Failed for %d particles, reducing stepsize to try again', sum(idxfail)));
        //end
        //
        //idxsuccess = ~idxfail;
        //
        //idxsuccessfull = false(size(posIn(:, 1)));
        //
        //idxsuccessfull(idxrunning) = idxsuccess;
        //
        //velOut(idxsuccessfull, :) = vel7(idxsuccess, :);
        //posOut(idxsuccessfull, :) = pos7(idxsuccess, :);
        //accelOut(idxsuccessfull, :) = accel7(idxsuccess, :);
        //dtOut(idxsuccessfull) = dt(idxsuccess);
        //
        //vel1 = vel1(idxfail, :);
        //pos1 = pos1(idxfail, :);
        //accel1 = accel1(idxfail, :);
        //dt = dt(idxfail);
        //
        //idxrunning(idxsuccessfull) = false;
        //
        //idxkill = dt < dt_min;
        //
        //if sum(idxkill) > 0
            //display(sprintf('            %d particles failed to converge at miminum stepsize; killing them.', sum(idxkill)));
//
            //idxkillfull = false(size(posIn(:, 1)));
            //idxkillfull(idxrunning) = idxkill;
//
            //vel1 = vel1(~idxkill, :);
            //pos1 = pos1(~idxkill, :);
            //accel1 = accel1(~idxkill, :);
            //dt = dt(~idxkill);
//
            //idxrunning(idxkillfull) = false;
            //
            //if sum(idxrunning) <= 0
                //break;
            //end
        //end
        //
        //nofailures = false;
    //end
//end
