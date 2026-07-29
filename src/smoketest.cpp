// Smoke test: confirm MPI + PETSc are usable on this cluster before any
// real solver code is written. No other project files are included.
#include <petsc.h>
#include <mpi.h>
#include <cstdio>
#include <cstring>
#include <set>
#include <string>
#include <vector>
#include <iostream>

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char hostname[MPI_MAX_PROCESSOR_NAME];
    int hostlen;
    MPI_Get_processor_name(hostname, &hostlen);

    // Gather every rank's hostname to rank 0 so we can confirm ranks actually
    // landed on more than one physical node, not just print per-rank claims.
    char padded_hostname[MPI_MAX_PROCESSOR_NAME] = {0};
    std::strncpy(padded_hostname, hostname, MPI_MAX_PROCESSOR_NAME - 1);
    std::vector<char> all_hostnames;
    if (rank == 0) all_hostnames.resize((size_t)size * MPI_MAX_PROCESSOR_NAME);
    MPI_Gather(padded_hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
               rank == 0 ? all_hostnames.data() : nullptr, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
               0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::set<std::string> distinct_nodes;
        for (int r = 0; r < size; ++r)
            distinct_nodes.insert(std::string(&all_hostnames[(size_t)r * MPI_MAX_PROCESSOR_NAME]));
        printf("=== node spread check: %d rank(s) across %zu distinct node(s) ===\n",
               size, distinct_nodes.size());
        for (const auto &node : distinct_nodes) printf("  node: %s\n", node.c_str());
        if (distinct_nodes.size() > 1)
            printf("=== CONFIRMED: MPI_COMM_WORLD spans multiple physical nodes ===\n");
        else
            printf("=== WARNING: all ranks landed on a single node -- not a multi-node test ===\n");
    }

    PetscInt major, minor, subminor;
    PetscGetVersionNumber(&major, &minor, &subminor, nullptr);

    // Trivial PETSc parallel op: distributed vector, set to rank id, take norm.
    Vec x;
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, size * 10);
    VecSetFromOptions(x);
    VecSet(x, (PetscScalar)(rank + 1));
    PetscReal norm;
    VecNorm(x, NORM_2, &norm);

    PetscPrintf(PETSC_COMM_WORLD,
                "PETSc %d.%d.%d OK | MPI ranks=%d | global Vec norm=%g\n",
                major, minor, subminor, size, (double)norm);
    printf("rank %d/%d on host %s\n", rank, size, hostname);

    VecDestroy(&x);
    PetscFinalize();

    std::cout << "main.cpp done running" << std::endl;
    return 0;
}
