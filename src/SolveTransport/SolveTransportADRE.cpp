#include "SolveTransportADRE.h"

// TODO: port SolveTransport/SolveTransportADRE.m -- one unsteady ADRE
// step for phi (COEFFPHIADRE/RHSPHIADRE assembly + KSP solve, per
// species). No-op stub for now: stateVar.phi is left unchanged.
void solveTransportADRE(ControlVar &controlVar, const Domain &domain,
                        const Variables &variables, StateVar &stateVar,
                        const IBM &ibm, const BC &bc, const LS &ls) {
    (void)controlVar; (void)domain; (void)variables; (void)stateVar;
    (void)ibm; (void)bc; (void)ls;
}
