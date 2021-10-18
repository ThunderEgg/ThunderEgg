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

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::tpl;

TEST_CASE("NormalNbrInfo getNbrType works", "[NormalNbrInfo]")
{
  NbrInfo<3>* info = new NormalNbrInfo<3>();
  REQUIRE(info->getNbrType() == NbrType::Normal);
  delete info;
}

TEST_CASE("NormalNbrInfo Serialization/Deserialization", "[NormalNbrInfo]")
{
  NormalNbrInfo<3> info;
  info.id = 5;
  info.rank = 1;
  // serialize and then deserialize
  char* buff = new char[info.serialize(nullptr)];
  info.serialize(buff);
  NormalNbrInfo<3> out;
  out.deserialize(buff);
  delete[] buff;
  REQUIRE(out.id == 5);
  REQUIRE(out.rank == 1);
}
TEST_CASE("NormalNbrInfo to_json", "[NormalNbrInfo]")
{
  NormalNbrInfo<3> info;
  info.id = GENERATE(1, 2, 3);
  info.rank = GENERATE(0, 1, 2);

  nlohmann::json j = info;

  CHECK(j["type"] == "NORMAL");
  REQUIRE(j["ids"].is_array());
  CHECK(j["ids"].size() == 1);
  CHECK(j["ids"][0] == info.id);
  REQUIRE(j["ranks"].is_array());
  CHECK(j["ranks"].size() == 1);
  CHECK(j["ranks"][0] == info.rank);
}
TEST_CASE("NormalNbrInfo from_json", "[NormalNbrInfo]")
{
  int id = GENERATE(1, 2, 3);
  int rank = GENERATE(0, 1, 2);

  nlohmann::json j;
  j["type"] = "NORMAL";
  j["ids"] = { id };
  j["ranks"] = { rank };

  NormalNbrInfo<3> info = j.get<NormalNbrInfo<3>>();
  CHECK(info.id == id);
  CHECK(info.rank == rank);
}