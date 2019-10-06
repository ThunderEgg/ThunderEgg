#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include <mpi.h>

int main(int argc, char *argv[])
{
	// global setup...
	MPI_Init(nullptr, nullptr);

	int result = Catch::Session().run(argc, argv);

	// global clean-up...
	MPI_Finalize();

	return result;
}