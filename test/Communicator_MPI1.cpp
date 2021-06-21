#include <ThunderEgg/Communicator.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Default Constructor getMPIComm throws", "[Communicator]")
{
	Communicator comm;
	CHECK_THROWS_AS(comm.getMPIComm(), RuntimeError);
}
TEST_CASE("Default Constructor copy constructor", "[Communicator]")
{
	Communicator comm;
	Communicator comm_copy(comm);
	CHECK_THROWS_AS(comm_copy.getMPIComm(), RuntimeError);
}
TEST_CASE("Default Constructor copy assignment", "[Communicator]")
{
	Communicator comm;
	Communicator comm_copy;
	comm_copy = comm;
	CHECK_THROWS_AS(comm_copy.getMPIComm(), RuntimeError);
}
TEST_CASE("Default Constructor getRank throws", "[Communicator]")
{
	Communicator comm;
	CHECK_THROWS_AS(comm.getRank(), RuntimeError);
}
TEST_CASE("Default Constructor getSize throws", "[Communicator]")
{
	Communicator comm;
	CHECK_THROWS_AS(comm.getSize(), RuntimeError);
}
TEST_CASE("Comm Constructor getMPIComm", "[Communicator]")
{
	MPI_Comm     world = MPI_COMM_WORLD;
	Communicator comm(world);
	int          result;
	int          err = MPI_Comm_compare(comm.getMPIComm(), world, &result);
	REQUIRE(err == MPI_SUCCESS);
	CHECK(result == MPI_CONGRUENT);
}
TEST_CASE("Comm Constructor copy constructor", "[Communicator]")
{
	MPI_Comm     world = MPI_COMM_WORLD;
	Communicator comm(world);
	Communicator comm_copy(comm);
	int          result;
	int          err = MPI_Comm_compare(comm.getMPIComm(), comm_copy.getMPIComm(), &result);
	REQUIRE(err == MPI_SUCCESS);
	CHECK(result == MPI_CONGRUENT);
}
TEST_CASE("Comm Constructor copy assignment", "[Communicator]")
{
	MPI_Comm     world = MPI_COMM_WORLD;
	Communicator comm(world);
	Communicator comm_copy;
	comm_copy = comm;
	int result;
	int err = MPI_Comm_compare(comm.getMPIComm(), comm_copy.getMPIComm(), &result);
	REQUIRE(err == MPI_SUCCESS);
	CHECK(result == MPI_CONGRUENT);
}
TEST_CASE("Comm Constructor getRank", "[Communicator]")
{
	Communicator comm(MPI_COMM_WORLD);
	CHECK(comm.getRank() == 0);
}
TEST_CASE("Default Constructor getSize", "[Communicator]")
{
	Communicator comm(MPI_COMM_WORLD);
	CHECK(comm.getSize() == 1);
}