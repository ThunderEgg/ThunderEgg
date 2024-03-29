/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2021      Scott Aiton
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/
#include <ThunderEgg/Communicator.h>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Default Constructor getMPIComm throws")
{
  Communicator comm;
  CHECK_THROWS_AS(comm.getMPIComm(), RuntimeError);
}
TEST_CASE("Default Constructor copy constructor")
{
  Communicator comm;
  Communicator comm_copy(comm);
  CHECK_THROWS_AS(comm_copy.getMPIComm(), RuntimeError);
}
TEST_CASE("Default Constructor copy assignment")
{
  Communicator comm;
  Communicator comm_copy;
  comm_copy = comm;
  CHECK_THROWS_AS(comm_copy.getMPIComm(), RuntimeError);
}
TEST_CASE("Default Constructor getRank throws")
{
  Communicator comm;
  CHECK_THROWS_AS(comm.getRank(), RuntimeError);
}
TEST_CASE("Default Constructor getSize throws")
{
  Communicator comm;
  CHECK_THROWS_AS(comm.getSize(), RuntimeError);
}
TEST_CASE("Comm Constructor getMPIComm")
{
  MPI_Comm world = MPI_COMM_WORLD;
  Communicator comm(world);
  int result;
  int err = MPI_Comm_compare(comm.getMPIComm(), world, &result);
  REQUIRE_EQ(err, MPI_SUCCESS);
  CHECK_EQ(result, MPI_CONGRUENT);
}
TEST_CASE("Comm Constructor copy constructor")
{
  MPI_Comm world = MPI_COMM_WORLD;
  Communicator comm(world);
  Communicator comm_copy(comm);
  int result;
  int err = MPI_Comm_compare(comm.getMPIComm(), comm_copy.getMPIComm(), &result);
  REQUIRE_EQ(err, MPI_SUCCESS);
  CHECK_EQ(result, MPI_IDENT);
}
TEST_CASE("Comm Constructor copy assignment")
{
  MPI_Comm world = MPI_COMM_WORLD;
  Communicator comm(world);
  Communicator comm_copy;
  comm_copy = comm;
  int result;
  int err = MPI_Comm_compare(comm.getMPIComm(), comm_copy.getMPIComm(), &result);
  REQUIRE_EQ(err, MPI_SUCCESS);
  CHECK_EQ(result, MPI_IDENT);
}
TEST_CASE("Comm Constructor move constructor")
{
  MPI_Comm world = MPI_COMM_WORLD;
  Communicator comm(world);
  Communicator moved_comm(std::move(comm));
  int result;
  int err = MPI_Comm_compare(world, moved_comm.getMPIComm(), &result);
  REQUIRE_EQ(err, MPI_SUCCESS);
  CHECK_EQ(result, MPI_CONGRUENT);
  //CHECK_THROWS_AS(comm.getMPIComm(), RuntimeError);
}
TEST_CASE("Comm Constructor move assignment")
{
  MPI_Comm world = MPI_COMM_WORLD;
  Communicator comm(world);
  Communicator moved_comm;
  moved_comm = std::move(comm);
  int result;
  int err = MPI_Comm_compare(world, moved_comm.getMPIComm(), &result);
  REQUIRE_EQ(err, MPI_SUCCESS);
  CHECK_EQ(result, MPI_CONGRUENT);
  //CHECK_THROWS_AS(comm.getMPIComm(), RuntimeError);
}
TEST_CASE("Comm Constructor getRank")
{
  Communicator comm(MPI_COMM_WORLD);
  CHECK_EQ(comm.getRank(), 0);
}
TEST_CASE("Default Constructor getSize")
{
  Communicator comm(MPI_COMM_WORLD);
  CHECK_EQ(comm.getSize(), 1);
}
