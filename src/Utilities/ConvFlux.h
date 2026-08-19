#pragma once
#include "../VariableNonDim.h"
#include "../ControlVar.h"

// Result of computing one equation's convective face fluxes -- the
// moving-fluid counterpart to DiffFluxCoeffs (Utilities/DiffuFlux.h).
struct ConvFluxCoeffs {
    // Face mass fluxes (velocity interpolated to each face, times face
    // area).
    Field2D Fe, Fw, Fn, Fs;

    // Upwind direction switches: 1 where the corresponding flux is
    // positive, 0 otherwise (0 also where the flux hasn't been computed
    // at all -- Field2D's zero default already gives the right answer
    // there, since an unset flux reads as 0, which this maps to 0 too).
    Field2D alphae, alphaw, alphan, alphas;

    // Net flux imbalance: Fe - Fw + Fn - Fs.
    Field2D dF;

    // accumulation coefficient: cell volume / dt, set to zero everywhere
    // when the equation is being solved as steady.
    Field2D A0_p;
};

// Computes U-momentum's convective face fluxes from the current velocity
// field. Interior range matches computeDiffFluxU's (i in [1,imax-2], j
// in [1,jmax-1], 0-indexed). If bc.BC_e_p==1, Fw/Fn/Fs additionally
// extend that same formula to the outlet row (i=imax-1); Fe instead
// falls back to a direct (non-interpolated) value there, since there's
// no further east neighbor to interpolate against.
ConvFluxCoeffs computeConvFluxU(const StateVar &stateVar, const ControlVar &controlVar,
                                 const Domain &domain, const Variables &variables, const BC &bc);

// Computes V-momentum's convective face fluxes. Interior range matches
// computeDiffFluxV's (i in [1,imax-1], j in [1,jmax-2], 0-indexed). No
// outlet-row extension exists for V.
ConvFluxCoeffs computeConvFluxV(const StateVar &stateVar, const ControlVar &controlVar,
                                 const Domain &domain, const Variables &variables, const BC &bc);

// Scalar transport (ConvFlux.m's flag==0 branch), sized (imax+1,
// jmax+1) like computeDiffFluxPhi. Computed once per
// SolveTransportADRE() call (velocity is frozen for the whole transport
// solve), not re-derived per QUICK iteration. Unlike U/V's own
// convective flux, there's no CoEW/CoNS interpolation and no BC-
// dependent special case -- U/V already live exactly at a P-cell's
// faces, so each face flux is just that face's raw U/V value times the
// face length, always, everywhere in the interior.
ConvFluxCoeffs computeConvFluxPhi(const StateVar &stateVar, const ControlVar &controlVar,
                                   const Domain &domain, const Variables &variables, const BC &bc);
