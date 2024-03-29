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

#include <ThunderEgg/Vector.h>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Vector<1> getNumGhostCells default constructor")
{
  Communicator comm(MPI_COMM_WORLD);

  Vector<1> vec;

  CHECK_EQ(vec.getNumGhostCells(), 0);
}
TEST_CASE("Vector<1> getMPIComm default constructor")
{
  Communicator comm(MPI_COMM_WORLD);

  Vector<1> vec;

  CHECK_THROWS_AS(vec.getCommunicator().getMPIComm(), RuntimeError);
}
TEST_CASE("Vector<1> getNumLocalPatches default constructor")
{
  Communicator comm(MPI_COMM_WORLD);

  Vector<1> vec;

  CHECK_EQ(vec.getNumLocalPatches(), 0);
}
TEST_CASE("Vector<1> getNumComponents default constructor")
{
  Communicator comm(MPI_COMM_WORLD);

  Vector<1> vec;

  CHECK_EQ(vec.getNumComponents(), 0);
}
TEST_CASE("Vector<1> getNumLocalCells default constructor")
{
  Communicator comm(MPI_COMM_WORLD);

  Vector<1> vec;

  CHECK_EQ(vec.getNumLocalCells(), 0);
}
TEST_CASE("Vector<2> getNumGhostCells default constructor")
{
  Communicator comm(MPI_COMM_WORLD);

  Vector<2> vec;

  CHECK_EQ(vec.getNumGhostCells(), 0);
}
TEST_CASE("Vector<2> getMPIComm default constructor")
{
  Vector<2> vec;

  CHECK_THROWS_AS(vec.getCommunicator().getMPIComm(), RuntimeError);
}
TEST_CASE("Vector<2> getNumLocalPatches default constructor")
{
  Vector<2> vec;

  CHECK_EQ(vec.getNumLocalPatches(), 0);
}
TEST_CASE("Vector<2> getNumComponents default constructor")
{
  Vector<2> vec;

  CHECK_EQ(vec.getNumComponents(), 0);
}
TEST_CASE("Vector<2> getNumLocalCells default constructor")
{
  Vector<2> vec;

  CHECK_EQ(vec.getNumLocalCells(), 0);
}
TEST_CASE("Vector<3> getNumGhostCells default constructor")
{
  Vector<3> vec;

  CHECK_EQ(vec.getNumGhostCells(), 0);
}
TEST_CASE("Vector<3> getMPIComm default constructor")
{
  Vector<3> vec;

  CHECK_THROWS_AS(vec.getCommunicator().getMPIComm(), RuntimeError);
}
TEST_CASE("Vector<3> getNumLocalPatches default constructor")
{
  Vector<3> vec;

  CHECK_EQ(vec.getNumLocalPatches(), 0);
}
TEST_CASE("Vector<3> getNumComponents default constructor")
{
  Vector<3> vec;

  CHECK_EQ(vec.getNumComponents(), 0);
}
TEST_CASE("Vector<3> getNumLocalCells default constructor")
{
  Vector<3> vec;

  CHECK_EQ(vec.getNumLocalCells(), 0);
}

TEST_CASE("Vector<1> getComponentView.h default constructor")
{
  Vector<1> vec;

  if (ENABLE_DEBUG) {
    CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, 0), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, 0), RuntimeError);
  }
}
TEST_CASE("Vector<1> getComponentView const default constructor")
{
  Vector<1> vec_non_const;
  const Vector<1>& vec = vec_non_const;

  if (ENABLE_DEBUG) {
    CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, 0), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, 0), RuntimeError);
  }
}
TEST_CASE("Vector<2> getComponentView default constructor")
{
  Vector<2> vec;

  if (ENABLE_DEBUG) {
    CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, 0), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, 0), RuntimeError);
  }
}
TEST_CASE("Vector<2> getComponentView const default constructor")
{
  Vector<2> vec_non_const;
  const Vector<2>& vec = vec_non_const;

  if (ENABLE_DEBUG) {
    CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, 0), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, 0), RuntimeError);
  }
}
TEST_CASE("Vector<3> getComponentView default constructor")
{
  Vector<3> vec;

  if (ENABLE_DEBUG) {
    CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, 0), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, 0), RuntimeError);
  }
}
TEST_CASE("Vector<3> getComponentView const default constructor")
{
  Vector<3> vec_non_const;
  const Vector<3>& vec = vec_non_const;

  if (ENABLE_DEBUG) {
    CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, 0), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
    CHECK_THROWS_AS(vec.getComponentView(0, 0), RuntimeError);
  }
}
TEST_CASE("Vector<1> getPatchView default constructor")
{
  Vector<1> vec;

  if (ENABLE_DEBUG) {
    CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
    CHECK_THROWS_AS(vec.getPatchView(0), RuntimeError);
  }
}
TEST_CASE("Vector<1> getPatchView const default constructor")
{
  Vector<1> vec_non_const;
  const Vector<1>& vec = vec_non_const;

  if (ENABLE_DEBUG) {
    CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
    CHECK_THROWS_AS(vec.getPatchView(0), RuntimeError);
  }
}
TEST_CASE("Vector<2> getPatchView default constructor")
{
  Vector<2> vec;

  if (ENABLE_DEBUG) {
    CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
    CHECK_THROWS_AS(vec.getPatchView(0), RuntimeError);
  }
}
TEST_CASE("Vector<2> getPatchView const default constructor")
{
  Vector<2> vec_non_const;
  const Vector<2>& vec = vec_non_const;

  if (ENABLE_DEBUG) {
    CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
    CHECK_THROWS_AS(vec.getPatchView(0), RuntimeError);
  }
}
TEST_CASE("Vector<3> getPatchView default constructor")
{
  Vector<3> vec;

  if (ENABLE_DEBUG) {
    CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
    CHECK_THROWS_AS(vec.getPatchView(0), RuntimeError);
  }
}
TEST_CASE("Vector<3> getPatchView const default constructor")
{
  Vector<3> vec_non_const;
  const Vector<3>& vec = vec_non_const;

  if (ENABLE_DEBUG) {
    CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
    CHECK_THROWS_AS(vec.getPatchView(0), RuntimeError);
  }
}
