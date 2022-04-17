/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
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
#include <ThunderEgg/PatchInfo.h>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("CoarseNbrInfo Serialization/Deserialization")
{
  CoarseNbrInfo<2> info;
  info.id = 5;
  info.rank = 1;
  info.orth_on_coarse = Orthant<2>::nw();
  // serialize and then deserialize
  char* buff = new char[info.serialize(nullptr)];
  info.serialize(buff);
  CoarseNbrInfo<2> out;
  out.deserialize(buff);
  delete[] buff;
  REQUIRE(out.id == 5);
  REQUIRE(out.rank == 1);
  REQUIRE(out.orth_on_coarse == Orthant<2>::nw());
}
