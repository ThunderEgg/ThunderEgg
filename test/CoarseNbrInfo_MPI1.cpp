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

TEST_CASE("CoarseNbrInfo Serialization/Deserialization", "[CoarseNbrInfo]")
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
TEST_CASE("CoarseNbrInfo to_json", "[CoarseNbrInfo]")
{
  CoarseNbrInfo<2> info;
  info.id = GENERATE(1, 2, 3);
  info.rank = GENERATE(0, 1, 2);
  info.orth_on_coarse = GENERATE(Orthant<2>::sw(), Orthant<3>::se(), Orthant<2>::nw());

  ThunderEgg::tpl::nlohmann::json j = info;

  CHECK(j["type"] == "COARSE");
  REQUIRE(j["ids"].is_array());
  CHECK(j["ids"].size() == 1);
  CHECK(j["ids"][0] == info.id);
  REQUIRE(j["ranks"].is_array());
  CHECK(j["ranks"].size() == 1);
  CHECK(j["ranks"][0] == info.rank);
  CHECK(j["orth_on_coarse"].get<Orthant<2>>() == info.orth_on_coarse);
}
TEST_CASE("CoarseNbrInfo from_json", "[CoarseNbrInfo]")
{
  int id = GENERATE(1, 2, 3);
  int rank = GENERATE(0, 1, 2);
  Orthant<2> orth_on_coarse = GENERATE(Orthant<2>::sw(), Orthant<3>::se(), Orthant<2>::nw());

  ThunderEgg::tpl::nlohmann::json j;
  j["type"] = "COARSE";
  j["ids"] = { id };
  j["ranks"] = { rank };
  j["orth_on_coarse"] = orth_on_coarse;

  CoarseNbrInfo<2> info = j.get<CoarseNbrInfo<2>>();
  CHECK(info.id == id);
  CHECK(info.rank == rank);
  CHECK(info.orth_on_coarse == orth_on_coarse);
}