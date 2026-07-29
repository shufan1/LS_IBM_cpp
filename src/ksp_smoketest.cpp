// Minimal, standalone test: does PETSc's own bicgstab KSP type resolve
// at all on this cluster/build, independent of anything else in this
// project? Isolates whether "Unable to find requested KSP type
// bicgstab" is a genuine environment/library issue or something
// specific to how the real program calls PETSc.
#include <petsc.h>
#include <cstdio>

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    Mat A;
    MatCreateSeqAIJ(PETSC_COMM_SELF, 3, 3, 3, nullptr, &A);
    for (int i = 0; i < 3; ++i) MatSetValue(A, i, i, 1.0, INSERT_VALUES);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    KSP ksp;
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, A, A);
    printf("about to call KSPSetType(bicgstab)...\n");
    fflush(stdout);
    PetscErrorCode ierr = KSPSetType(ksp, KSPBCGS);
    printf("KSPSetType returned %d\n", (int)ierr);

    KSPDestroy(&ksp);
    MatDestroy(&A);
    PetscFinalize();
    return 0;
}
