#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>
#include <mpi.h>
#if TEST_P4EST
#include <sc.h>
#endif
#if TEST_PETSC
#include <petscsys.h>
#endif

int main(int argc, char *argv[])
{
#if TEST_P4EST
	sc_set_log_defaults(NULL, NULL, SC_LP_SILENT);
#endif
	// global setup...
#if TEST_PETSC
	PetscInitialize(nullptr, nullptr, nullptr, nullptr);
#else
	MPI_Init(nullptr, nullptr);
#endif

	int result = Catch::Session().run(argc, argv);

	// abort if failure, some tests can hang otherwise
	if (result > 0) {
		MPI_Abort(MPI_COMM_WORLD, result);
	}

	// global clean-up...
#if TEST_PETSC
	PetscFinalize();
#else
	MPI_Finalize();
#endif
	return result;
}
