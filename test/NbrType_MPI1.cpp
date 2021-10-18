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
#include <ThunderEgg/NbrType.h>

#include <catch2/catch_test_macros.hpp>

using namespace ThunderEgg;
using namespace ThunderEgg::tpl;
TEST_CASE("NbrType to_json", "[NbrType]")
{
  nlohmann::json j;
  j["normal"] = NbrType::Normal;
  j["coarse"] = NbrType::Coarse;
  j["fine"] = NbrType::Fine;
  CHECK(j["normal"] == "NORMAL");
  CHECK(j["coarse"] == "COARSE");
  CHECK(j["fine"] == "FINE");
}
TEST_CASE("NbrType from_json", "[NbrType]")
{
  nlohmann::json j;
  j["normal"] = "NORMAL";
  j["coarse"] = "COARSE";
  j["fine"] = "FINE";
  CHECK(j["normal"].get<NbrType>() == NbrType::Normal);
  CHECK(j["coarse"].get<NbrType>() == NbrType::Coarse);
  CHECK(j["fine"].get<NbrType>() == NbrType::Fine);
}