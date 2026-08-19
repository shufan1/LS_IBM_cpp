#include "LSVelocityExtrapolation.h"
#include "../Utilities/Interp.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

// Ported through LSVelocityExtrapolation.m:137 (`q_mlout = q_mlout'`).
// Step 5 below and the scatter back onto LS are still TODO -- see the note
// at the bottom of the function.
void LSVelocityExtrapolation(LS &ls, const StateVar &stateVar, const Domain &domain,
                              const Variables &variables)
{// For every node in the narrow band |psi| < variables.LSgamma, on BOTH
// sides of the interface (solveHJEq needs a velocity on both to advect
// psi):
//
//   1. Project onto the interface along the normal: x_I = x -/+ n*psi
//
//   2. Place two probes further out into the FLUID, at sqrt(2)*dx and
//      2*sqrt(2)*dx from x_I along n. Both probes are on the fluid side
//      regardless of which side the node itself is on.
//   3. Bilinearly interpolate all Np species at both probes
//      (Utilities/Interp.h; MATLAB uses biLinearInterpolation3d).
//   4. Solve for the wall concentration with the kinetics applied:
//         phi_wall = inv_A * (4*phi_probe1 - phi_probe2)
//      which is the second-order one-sided derivative stencil set equal to
//      the reaction rate, i.e. dphi/dn = A*phi_wall -- diffusive flux
//      arriving at the surface equals the reaction consuming it. One 3x3
//      solve per band node.
//   5. Convert species 1 (0-based; MATLAB's q_mlout(:,2)) molar flux to a
//      recession speed via the mineral's molar volume, then scale by the
//      normal: ls.u = speed*nx, ls.v = speed*ny.

    const int Nx = ls.psi.nx();  // imax+1
    const int Ny = ls.psi.ny();  // jmax+1
    const int Np = variables.Np;

    // ---- Outputs, zeroed up front (LSVelocityExtrapolation.m:29) ----
    // Nodes outside the narrow band keep these zeros; only band nodes get
    // written. MATLAB also zeroes LS.q_out here (line 11), which is not
    // ported -- see struct LS.
    ls.u = Field2D(Nx, Ny);
    ls.v = Field2D(Nx, Ny);

    const Field2D &psi = ls.psi;
    const Field2D &nx = ls.nx;
    const Field2D &ny = ls.ny;
    const std::vector<double> &x = domain.xp;
    const std::vector<double> &y = domain.yp;

    const double gamma = variables.LSgamma;
    const double dx = *std::min_element(domain.dxp.begin(), domain.dxp.end());
    const double dy = *std::min_element(domain.dyp.begin(), domain.dyp.end());

    // ---- i-sweep range (LSVelocityExtrapolation.m:47-56) ----
    // Cases 3/4 (rough fracture) clip the i-range to the inlet/outlet
    // buffer using DOMAIN.x_0. Neither is ported -- Domain has no x_0 --
    // so reject them rather than silently sweeping the wrong range.
    if (ls.caseId == 3 || ls.caseId == 4) {
        throw std::runtime_error(
            "LSVelocityExtrapolation: LSCase 3/4 (rough fracture) clips the i-range with "
            "DOMAIN.x_0, which is not ported -- only caseId==1 (grain) is supported");
    }
    const int iBegin = 1, iEnd = Nx - 2;  // MATLAB I_s=2 .. I_f=Nx-1, 1-based
    const int jBegin = 1, jEnd = Ny - 2;  // MATLAB j=2 .. Ny-1

    // ---- Batch storage ----
    // Flat, laid out [k*Np + s] -- node-major, species-minor. That is
    // already MATLAB's *transposed* (N_band x Np) orientation, so the
    // `phi_save1 = phi_save1'` and `q_mlout = q_mlout'` transposes at
    // lines 120-121 and 137 cost nothing here.
    //
    // The precount is a reserve hint only. MATLAB sizes its batch arrays
    // to N_band counted over the WHOLE array while the sweep covers only
    // iBegin..iEnd, then runs the matrix products across the unused
    // trailing columns anyway and discards them. Here the batch holds
    // exactly the nodes actually found.
    size_t bandHint = 0;
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            if (std::abs(psi(i, j)) < gamma) ++bandHint;

    std::vector<int> bandI, bandJ;
    std::vector<double> phiProbe1, phiProbe2;
    bandI.reserve(bandHint);
    bandJ.reserve(bandHint);
    phiProbe1.reserve(bandHint * static_cast<size_t>(Np));
    phiProbe2.reserve(bandHint * static_cast<size_t>(Np));

    const double SQRT2 = std::sqrt(2.0);

    // ---- Steps 1-3: narrow-band sweep (lines 60-118) ----
    for (int i = iBegin; i <= iEnd; ++i) {
        for (int j = jBegin; j <= jEnd; ++j) {
            if (std::abs(psi(i, j)) >= gamma) continue;

            // 1. MATLAB splits this into psi>=0 and psi<0 branches
            // (`x - n*psi` and `x + n*abs(psi)`), but the two are the same
            // expression -- abs(psi) == -psi when psi<0 -- so it is
            // written once here.
            const double x_I = x[i] - nx(i, j) * psi(i, j);
            const double y_I = y[j] - ny(i, j) * psi(i, j);

            // 2. x steps by dx and y by dy, as MATLAB does. Identical on
            // the uniform grid; kept distinct for fidelity.
            const double x_pr = x_I + SQRT2 * dx * nx(i, j);
            const double y_pr = y_I + SQRT2 * dy * ny(i, j);
            const double x_dpr = x_pr + SQRT2 * dx * nx(i, j);
            const double y_dpr = y_pr + SQRT2 * dy * ny(i, j);
            // (MATLAB also sets `d = sqrt(2)*dy` in both branches and
            // never reads it again -- dropped.)

            bandI.push_back(i);
            bandJ.push_back(j);

            // 3. biLinearInterpolation3d does all Np species in one call;
            // stateVar.phi is one Field2D per species, so this is a loop
            // over the same bilinear interpolation.
            for (int s = 0; s < Np; ++s) {
                phiProbe1.push_back(bilinearInterp(x, y, stateVar.phi[s], x_pr, y_pr));
                phiProbe2.push_back(bilinearInterp(x, y, stateVar.phi[s], x_dpr, y_dpr));
            }
        }
    }

    const size_t nBand = bandI.size();

    // ---- Step 4: wall concentration (line 135) ----
    //   phi_wall = inv_A * (4*phi_probe1 - phi_probe2)
    // inv_A is inv(3*I - 2*l*A) with l = sqrt(2)*dx, precomputed once in
    // defineReactivity(), so this is a matrix-vector product per node --
    // no per-node solve.
    std::vector<double> phiWall(nBand * static_cast<size_t>(Np), 0.0);
    for (size_t k = 0; k < nBand; ++k) {
        for (int r = 0; r < Np; ++r) {
            double acc = 0.0;
            for (int c = 0; c < Np; ++c) {
                acc += variables.inv_A[r * Np + c] *
                       (4.0 * phiProbe1[k * Np + c] - phiProbe2[k * Np + c]);
            }
            phiWall[k * Np + r] = acc;
        }
    }

    // ---- Molar flux at the wall (lines 136-137) ----
    //   q_mlout = -A * phi_wall,  then transposed to (N_band x Np)
    // -- already in that orientation given the [k*Np + s] layout.
    //           node 0          node 1          node 2
    // qWall =  │s0│s1│s2│      │s0│s1│s2│      │s0│s1│s2│  ...
    std::vector<double> qWall(nBand * static_cast<size_t>(Np), 0.0);
    for (size_t k = 0; k < nBand; ++k) {
        for (int r = 0; r < Np; ++r) {
            double acc = 0.0;
            for (int c = 0; c < Np; ++c) acc += variables.A[r * Np + c] * phiWall[k * Np + c];
            qWall[k * Np + r] = -acc;
        }
    }

    // ---- Step 5: recession speed, then scatter (lines 139, 143-146) ----
    //   u_I_out = q_mlout(:,2)/Pe * c_real * 36.9e-6 / 1.261
    //
    // Only one species' flux moves the interface -- the mineral being
    // dissolved -- so this reads a single column out of qWall's
    // (nBand x Np) layout and ignores the other Np-1. Index 1 (0-based)
    // is MATLAB's q_mlout(:,2). Its molar flux times the mineral's molar
    // volume is a volume of solid removed per unit area per unit time,
    // which is a velocity; c_real and 1/Pe carry that back into the
    // nondimensional units psi lives in.
    const int kSpeciesRecession = 1; // the second species is the one being disolved

    // MATLAB's literal 36.9e-6 is VARIABLES.molarVol (36.9 cm^3/mol,
    // calcite) expressed in m^3/mol.
    const double molarVol_m3_per_mol = variables.molarVol * 1e-6;

    // Undo defineReactivity()'s `A = inv(Pd)*A` for this species. The
    // scaled A is right for the transport BC -- each species' equation is
    // nondimensionalized by its own diffusivity -- but a recession speed
    // needs the TRUE molar flux, which is diffusivity-independent.
    //
    // MATLAB writes this as a bare `/1.261` literal, which is 1/0.793
    // rounded to four figures, i.e. 1/Pd(2,2). Multiplying by Pd here
    // instead makes the cancellation exact and keeps it correct if Pd
    // ever changes. It also costs bit-exactness against MATLAB: 1/1.261
    // = 0.7930214 vs 0.793, a 2.7e-5 relative difference in the interface
    // speed. See the note in the header.
    const double undoPdScaling = variables.Pd[kSpeciesRecession];

    for (size_t k = 0; k < nBand; ++k) {
        const double speed = qWall[k * Np + kSpeciesRecession] / variables.Pe * variables.c_real *
                             molarVol_m3_per_mol * undoPdScaling;
        const int i = bandI[k];
        const int j = bandJ[k];
        ls.u(i, j) = speed * nx(i, j);
        ls.v(i, j) = speed * ny(i, j);
    }

}
