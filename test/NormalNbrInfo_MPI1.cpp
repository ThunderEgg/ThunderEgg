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

TEST_CASE("NormalNbrInfo getNbrType works")
{
  NbrInfo<2>* info = new NormalNbrInfo<2>();
  REQUIRE_EQ(info->getNbrType(), NbrType::Normal);
  delete info;
}

TEST_CASE("NormalNbrInfo Serialization/Deserialization")
{
  NormalNbrInfo<2> info;
  info.id = 5;
  info.rank = 1;
  // serialize and then deserialize
  char* buff = new char[info.serialize(nullptr)];
  info.serialize(buff);
  NormalNbrInfo<2> out;
  out.deserialize(buff);
  delete[] buff;
  REQUIRE_EQ(out.id, 5);
  REQUIRE_EQ(out.rank, 1);
}
