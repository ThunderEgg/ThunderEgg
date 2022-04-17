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
using namespace ThunderEgg::tpl;

TEST_CASE("FineNbrInfo Serialization/Deserialization")
{
  FineNbrInfo<2> info;
  info.ids[0] = 1;
  info.ids[1] = 2;
  info.ids[2] = 3;
  info.ids[3] = 4;
  info.ranks[0] = 9;
  info.ranks[1] = 8;
  info.ranks[2] = 7;
  info.ranks[3] = 6;
  // serialize and then deserialize
  char* buff = new char[info.serialize(nullptr)];
  info.serialize(buff);
  FineNbrInfo<2> out;
  out.deserialize(buff);
  delete[] buff;
  REQUIRE(out.ids[0] == 1);
  REQUIRE(out.ids[1] == 2);
  REQUIRE(out.ids[2] == 3);
  REQUIRE(out.ids[3] == 4);
  REQUIRE(out.ranks[0] == 9);
  REQUIRE(out.ranks[1] == 8);
  REQUIRE(out.ranks[2] == 7);
  REQUIRE(out.ranks[3] == 6);
}
