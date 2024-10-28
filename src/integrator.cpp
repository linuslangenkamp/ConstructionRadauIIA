/**
 * GDOPT - General Dynamic Optimizer
 * Copyright (C) 2024  Linus Langenkamp
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/

/*
Implementation of Radau IIA collocation schemes. Contains schemes for s = 1, 2, ..., 36 steps
with corresponding orders p = 2s - 1. The coefficients are calculated to 53 exact decimal places via
https://github.com/linuslangenkamp/ConstructionRadauIIA. The theoretical values
a_{s,s} = 1/s^2 and ainv_{s,s} = 1/2 (1 + s^2) hold.
*/

#include "integrator.h"

#include "util.h"
// TODO: remove A
// Outsource calculations if possible
// Note: A^-1 = D = diff Lagrange, differentiation matrix

Integrator::Integrator(const std::vector<double>& c, const std::vector<double>& c0,
    const std::vector<double>& cBisection, const std::vector<double>& c0Bisection, const std::vector<double>& b, 
    const std::vector<std::vector<double>>& derMat, const std::vector<std::vector<double>>& derMat2,
    const std::vector<std::vector<double>>& lagrangeBisectedFactorsC, const std::vector<std::vector<double>>& lagrangeBisectionC0, 
    int steps) 
    : c(c), c0(c0), cBisection(cBisection), c0Bisection(c0Bisection), b(b), derMat(derMat), derMat2(derMat2),
    lagrangeBisectedFactorsC(lagrangeBisectedFactorsC), lagrangeBisectionC0(lagrangeBisectionC0), steps(steps)
{}

// First and second derivative of the Lagrange interpolating polynomial on the nominal interval [0, 1]
// Will be used for detecting discontinuities, corners, sections that are steep or have a huge curvature
// use this for every interval, but the 0-th control interval

double Integrator::integrate(std::vector<double>& f) {
    // input: vector of f(c_j) excluding c_0 = 0
    // integrates according to Radau scheme
    double integral = 0;
    for (int j = 0; j < sz(b); j++) {
        integral += b[j] * f[j];
    }
    return integral;
}


std::vector<double> Integrator::evalLagrangeDiff(std::vector<double>& coefficients) {
    // evaluates the diff of the lagrange polynomial of given coefficients at c0 = 0, c1, c2, ...
    std::vector<double> lagrangeDiff;
    for (int i = 0; i < sz(derMat); i++) {
        double diffC = 0;
        for (int j = 0; j < sz(coefficients); j++) {
            diffC += derMat[i][j] * coefficients[j];
        }
        lagrangeDiff.push_back(diffC);
    }
    return lagrangeDiff;
}

std::vector<double> Integrator::evalLagrangeDiff2(std::vector<double>& coefficients) {
    // evaluates the 2nd diff of the lagrange polynomial of given coefficients at c0 = 0, c1, c2, ...
    std::vector<double> lagrangeDiff;
    for (int i = 0; i < sz(derMat2); i++) {
        double diffC = 0;
        for (int j = 0; j < sz(coefficients); j++) {
            diffC += derMat2[i][j] * coefficients[j];
        }
        lagrangeDiff.push_back(diffC);
    }
    return lagrangeDiff;
}

double Integrator::evalLagrange(std::vector<double> grid, std::vector<double>& f, double x) {
    // evaluates the interpolating polynomial p with p(grid[i]) = f[i], at x -> returns p(x)
    // runtime O(n^2) with n = #gridpoints
    double val = 0;
    for (int j = 0; j < sz(grid); j++) {
        double basisFactor = 1;
        for (int m = 0; m < sz(grid); m++) {
            if (m != j) {
                basisFactor *= (x - grid[m]) / (grid[j] - grid[m]);
            }
        }
        val += basisFactor * f[j];
    }
    return val;
}

// output values at c_0/2, c_1/2, ..., c_m/2 = 1/2, 1/2 + c_0/2, 1/2 + c_1/2, ..., 1
std::vector<double> Integrator::interpolateFirstControl(std::vector<double>& uValues) {
    std::vector<double> vals;
    for (auto coeffs : lagrangeBisectedFactorsC) {
        double sum = 0;
        for (int k = 0; k < steps; k++) {
            sum += uValues[k] * coeffs[k];
        }
        vals.push_back(sum);
    }
    return vals;
}

// output values of given interpolating polynomial
// at  c_0/2, c_1/2, ..., c_m/2 = 1/2, 1/2 + c_0/2, 1/2 + c_1/2, ..., 1
std::vector<double> Integrator::evalInterpolationNewNodes(std::vector<double>& values) {
    std::vector<double> vals;
    for (int j = 1; j < sz(lagrangeBisectionC0); j++) {
        double sum = 0;
        for (int k = 0; k < steps + 1; k++) {
            sum += values[k] * lagrangeBisectionC0[j][k];
        }
        vals.push_back(sum);
    }
    return vals;
}

std::vector<double> Integrator::evalLinearSplineNewNodes(std::vector<double>& values) {
    // evaluates the values input array at the new nodes c1/2, ..., cm/2=1/2, c1/2 + 1/2, ...
    // via a linear spline on the entire interval
    std::vector<double> newVals;
    int idx = 0;

    for (double cn : cBisection) {
        while (cn > c0[idx + 1] && idx < sz(c0) - 2) {
            idx++;
        }

        // cn must be inside c0[idx], c0[idx+1]
        if (c0[idx] <= cn && cn <= c0[idx + 1]) {
            // linear interpolation
            const double slope = (values[idx + 1] - values[idx]) / (c0[idx + 1] - c0[idx]);
            const double interpolatedVal = values[idx] + slope * (cn - c0[idx]);
            newVals.push_back(interpolatedVal);
        }
    }

    // remove val at t=c0
    if (!newVals.empty()) {
        newVals.erase(newVals.begin());
    }

    return newVals;
}

Integrator Integrator::radauIIA(IntegratorSteps steps) {

    switch (steps) {
        
    default:
        return {{1.0},
{0.0,1.0},
{0.5,1.0},
{0.0,0.5,1.0},
{1.0},
{{-1.0,1.0},
{-1.0,1.0}},
{{-2.0,2.},
{-2.0,2.}},
{{1.},
{1.}},
{{1.0,0.0},
{0.5,0.5},
{0.0,1.0}},
1};
    }
}
