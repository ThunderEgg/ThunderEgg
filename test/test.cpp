#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include <petsc.h>

int main(int argc, char *argv[])
{
	// global setup...
	//PetscInitialize(nullptr, nullptr, nullptr, nullptr);

	int result = Catch::Session().run(argc, argv);

	// global clean-up...
	//PetscFinalize();

	return result;
}