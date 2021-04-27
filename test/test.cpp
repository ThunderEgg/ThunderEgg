#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>
#include <mpi.h>
#if TEST_P4EST
#include <sc.h>
#endif

int main(int argc, char *argv[])
{
#if TEST_P4EST
	sc_set_log_defaults(NULL, NULL, SC_LP_SILENT);
#endif
	// global setup...
	MPI_Init(nullptr, nullptr);

	int result = Catch::Session().run(argc, argv);

	// abort if failure, some tests can hang otherwise
	if (result > 0) {
		MPI_Abort(MPI_COMM_WORLD, result);
	}

	// global clean-up...
	MPI_Finalize();

	return result;
}
