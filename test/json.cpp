/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2019-2022 Scott Aiton
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

#include <ThunderEgg/Face.h>
#include <ThunderEgg/Orthant.h>
#include <ThunderEgg/tpl/json.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace ThunderEgg;
using namespace ThunderEgg::tpl;

TEST_CASE("Test from_json for Edge", "[Edge][Face]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["bs"] = "BS";
  j["tn"] = "TN";
  j["bn"] = "BN";
  j["ts"] = "TS";
  j["bw"] = "BW";
  j["te"] = "TE";
  j["be"] = "BE";
  j["tw"] = "TW";
  j["sw"] = "SW";
  j["ne"] = "NE";
  j["se"] = "SE";
  j["nw"] = "NW";
  CHECK(j["null"].get<Edge>() == Edge::null());
  CHECK(j["bs"].get<Edge>() == Edge::bs());
  CHECK(j["tn"].get<Edge>() == Edge::tn());
  CHECK(j["bn"].get<Edge>() == Edge::bn());
  CHECK(j["ts"].get<Edge>() == Edge::ts());
  CHECK(j["bw"].get<Edge>() == Edge::bw());
  CHECK(j["te"].get<Edge>() == Edge::te());
  CHECK(j["be"].get<Edge>() == Edge::be());
  CHECK(j["tw"].get<Edge>() == Edge::tw());
  CHECK(j["sw"].get<Edge>() == Edge::sw());
  CHECK(j["ne"].get<Edge>() == Edge::ne());
  CHECK(j["nw"].get<Edge>() == Edge::nw());
}
TEST_CASE("Test to_json for Edge", "[Edge][Face]")
{
  nlohmann::json j;
  j["null"] = Edge::null();
  j["bs"] = Edge::bs();
  j["tn"] = Edge::tn();
  j["bn"] = Edge::bn();
  j["ts"] = Edge::ts();
  j["bw"] = Edge::bw();
  j["te"] = Edge::te();
  j["be"] = Edge::be();
  j["tw"] = Edge::tw();
  j["sw"] = Edge::sw();
  j["ne"] = Edge::ne();
  j["se"] = Edge::se();
  j["nw"] = Edge::nw();
  CHECK(j["null"] == nullptr);
  CHECK(j["bs"] == "BS");
  CHECK(j["tn"] == "TN");
  CHECK(j["bn"] == "BN");
  CHECK(j["ts"] == "TS");
  CHECK(j["bw"] == "BW");
  CHECK(j["te"] == "TE");
  CHECK(j["be"] == "BE");
  CHECK(j["tw"] == "TW");
  CHECK(j["sw"] == "SW");
  CHECK(j["ne"] == "NE");
  CHECK(j["se"] == "SE");
  CHECK(j["nw"] == "NW");
}
TEST_CASE("Test to_json for Side<1>", "[Side][Face]")
{
  nlohmann::json j;
  j["null"] = Side<1>::null();
  j["west"] = Side<1>::west();
  j["east"] = Side<1>::east();
  CHECK(j["null"] == nullptr);
  CHECK(j["west"] == "WEST");
  CHECK(j["east"] == "EAST");
}
TEST_CASE("Test to_json for Side<2>", "[Side][Face]")
{
  nlohmann::json j;
  j["null"] = Side<2>::null();
  j["west"] = Side<2>::west();
  j["east"] = Side<2>::east();
  j["south"] = Side<2>::south();
  j["north"] = Side<2>::north();
  CHECK(j["null"] == nullptr);
  CHECK(j["west"] == "WEST");
  CHECK(j["east"] == "EAST");
  CHECK(j["south"] == "SOUTH");
  CHECK(j["north"] == "NORTH");
}
TEST_CASE("Test to_json for Side<3>", "[Side][Face]")
{
  nlohmann::json j;
  j["null"] = Side<3>::null();
  j["west"] = Side<3>::west();
  j["east"] = Side<3>::east();
  j["south"] = Side<3>::south();
  j["north"] = Side<3>::north();
  j["bottom"] = Side<3>::bottom();
  j["top"] = Side<3>::top();
  CHECK(j["null"] == nullptr);
  CHECK(j["west"] == "WEST");
  CHECK(j["east"] == "EAST");
  CHECK(j["south"] == "SOUTH");
  CHECK(j["north"] == "NORTH");
  CHECK(j["bottom"] == "BOTTOM");
  CHECK(j["top"] == "TOP");
}
TEST_CASE("Test from_json for Side<1>", "[Side][Face]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["west"] = "WEST";
  j["east"] = "EAST";
  CHECK(j["null"].get<Side<1>>() == Side<1>::null());
  CHECK(j["west"].get<Side<1>>() == Side<1>::west());
  CHECK(j["east"].get<Side<1>>() == Side<1>::east());
}
TEST_CASE("Test from_json for Side<2>", "[Side][Face]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["west"] = "WEST";
  j["east"] = "EAST";
  j["south"] = "SOUTH";
  j["north"] = "NORTH";
  CHECK(j["null"].get<Side<2>>() == Side<2>::null());
  CHECK(j["west"].get<Side<2>>() == Side<2>::west());
  CHECK(j["east"].get<Side<3>>() == Side<3>::east());
  CHECK(j["south"].get<Side<2>>() == Side<2>::south());
  CHECK(j["north"].get<Side<2>>() == Side<2>::north());
}
TEST_CASE("Test from_json for Side<3>", "[Side][Face]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["west"] = "WEST";
  j["east"] = "EAST";
  j["south"] = "SOUTH";
  j["north"] = "NORTH";
  j["bottom"] = "BOTTOM";
  j["top"] = "TOP";
  CHECK(j["null"].get<Side<3>>() == Side<3>::null());
  CHECK(j["west"].get<Side<3>>() == Side<3>::west());
  CHECK(j["east"].get<Side<3>>() == Side<3>::east());
  CHECK(j["south"].get<Side<3>>() == Side<3>::south());
  CHECK(j["north"].get<Side<3>>() == Side<3>::north());
  CHECK(j["bottom"].get<Side<3>>() == Side<3>::bottom());
  CHECK(j["top"].get<Side<3>>() == Side<3>::top());
}
TEST_CASE("Test from_json for Orthant<0>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  CHECK(j["null"].get<Orthant<0>>() == Orthant<0>::null());
}
TEST_CASE("Test from_json for Orthant<1>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["lower"] = "LOWER";
  j["upper"] = "UPPER";
  CHECK(j["null"].get<Orthant<1>>() == Orthant<1>::null());
  CHECK(j["lower"].get<Orthant<1>>() == Orthant<1>::lower());
  CHECK(j["upper"].get<Orthant<1>>() == Orthant<1>::upper());
}
TEST_CASE("Test from_json for Orthant<2>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["sw"] = "SW";
  j["se"] = "SE";
  j["nw"] = "NW";
  j["ne"] = "NE";
  CHECK(j["null"].get<Orthant<2>>() == Orthant<2>::null());
  CHECK(j["sw"].get<Orthant<2>>() == Orthant<2>::sw());
  CHECK(j["se"].get<Orthant<2>>() == Orthant<2>::se());
  CHECK(j["nw"].get<Orthant<2>>() == Orthant<2>::nw());
  CHECK(j["ne"].get<Orthant<2>>() == Orthant<2>::ne());
}
TEST_CASE("Test from_json for Orthant<3>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["bsw"] = "BSW";
  j["bse"] = "BSE";
  j["bnw"] = "BNW";
  j["bne"] = "BNE";
  j["tsw"] = "TSW";
  j["tse"] = "TSE";
  j["tnw"] = "TNW";
  j["tne"] = "TNE";
  CHECK(j["null"].get<Orthant<3>>() == Orthant<3>::null());
  CHECK(j["bsw"].get<Orthant<3>>() == Orthant<3>::bsw());
  CHECK(j["bse"].get<Orthant<3>>() == Orthant<3>::bse());
  CHECK(j["bnw"].get<Orthant<3>>() == Orthant<3>::bnw());
  CHECK(j["bne"].get<Orthant<3>>() == Orthant<3>::bne());
  CHECK(j["tsw"].get<Orthant<3>>() == Orthant<3>::tsw());
  CHECK(j["tse"].get<Orthant<3>>() == Orthant<3>::tse());
  CHECK(j["tnw"].get<Orthant<3>>() == Orthant<3>::tnw());
  CHECK(j["tne"].get<Orthant<3>>() == Orthant<3>::tne());
}
TEST_CASE("Test to_json for Orthant<0>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = Orthant<0>::null();
  CHECK(j["null"] == nullptr);
}
TEST_CASE("Test to_json for Orthant<1>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = Orthant<1>::null();
  j["lower"] = Orthant<1>::lower();
  j["upper"] = Orthant<1>::upper();
  CHECK(j["null"] == nullptr);
  CHECK(j["lower"] == "LOWER");
  CHECK(j["upper"] == "UPPER");
}
TEST_CASE("Test to_json for Orthant<2>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = Orthant<2>::null();
  j["sw"] = Orthant<2>::sw();
  j["se"] = Orthant<2>::se();
  j["nw"] = Orthant<2>::nw();
  j["ne"] = Orthant<2>::ne();
  CHECK(j["null"] == nullptr);
  CHECK(j["sw"] == "SW");
  CHECK(j["se"] == "SE");
  CHECK(j["nw"] == "NW");
  CHECK(j["ne"] == "NE");
}
TEST_CASE("Test to_json for Orthant<3>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = Orthant<3>::null();
  j["bsw"] = Orthant<3>::bsw();
  j["bse"] = Orthant<3>::bse();
  j["bnw"] = Orthant<3>::bnw();
  j["bne"] = Orthant<3>::bne();
  j["tsw"] = Orthant<3>::tsw();
  j["tse"] = Orthant<3>::tse();
  j["tnw"] = Orthant<3>::tnw();
  j["tne"] = Orthant<3>::tne();
  CHECK(j["null"] == nullptr);
  CHECK(j["bsw"] == "BSW");
  CHECK(j["bse"] == "BSE");
  CHECK(j["bnw"] == "BNW");
  CHECK(j["bne"] == "BNE");
  CHECK(j["tsw"] == "TSW");
  CHECK(j["tse"] == "TSE");
  CHECK(j["tnw"] == "TNW");
  CHECK(j["tne"] == "TNE");
}
TEST_CASE("Test from_json for Corner<2>", "[Corner][Face]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["sw"] = "SW";
  j["se"] = "SE";
  j["nw"] = "NW";
  j["ne"] = "NE";
  CHECK(j["null"].get<Corner<2>>() == Corner<2>::null());
  CHECK(j["sw"].get<Corner<2>>() == Corner<2>::sw());
  CHECK(j["se"].get<Corner<2>>() == Corner<2>::se());
  CHECK(j["nw"].get<Corner<2>>() == Corner<2>::nw());
  CHECK(j["ne"].get<Corner<2>>() == Corner<2>::ne());
}
TEST_CASE("Test from_json for Corner<3>", "[Corner][Face]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["bsw"] = "BSW";
  j["bse"] = "BSE";
  j["bnw"] = "BNW";
  j["bne"] = "BNE";
  j["tsw"] = "TSW";
  j["tse"] = "TSE";
  j["tnw"] = "TNW";
  j["tne"] = "TNE";
  CHECK(j["null"].get<Corner<3>>() == Corner<3>::null());
  CHECK(j["bsw"].get<Corner<3>>() == Corner<3>::bsw());
  CHECK(j["bse"].get<Corner<3>>() == Corner<3>::bse());
  CHECK(j["bnw"].get<Corner<3>>() == Corner<3>::bnw());
  CHECK(j["bne"].get<Corner<3>>() == Corner<3>::bne());
  CHECK(j["tsw"].get<Corner<3>>() == Corner<3>::tsw());
  CHECK(j["tse"].get<Corner<3>>() == Corner<3>::tse());
  CHECK(j["tnw"].get<Corner<3>>() == Corner<3>::tnw());
  CHECK(j["tne"].get<Corner<3>>() == Corner<3>::tne());
}
TEST_CASE("Test to_json for Corner<2>", "[Corner][Face]")
{
  nlohmann::json j;
  j["null"] = Corner<2>::null();
  j["sw"] = Corner<2>::sw();
  j["se"] = Corner<2>::se();
  j["nw"] = Corner<2>::nw();
  j["ne"] = Corner<2>::ne();
  CHECK(j["null"] == nullptr);
  CHECK(j["sw"] == "SW");
  CHECK(j["se"] == "SE");
  CHECK(j["nw"] == "NW");
  CHECK(j["ne"] == "NE");
}
TEST_CASE("Test to_json for Corner<3>", "[Corner][Face]")
{
  nlohmann::json j;
  j["null"] = Corner<3>::null();
  j["bsw"] = Corner<3>::bsw();
  j["bse"] = Corner<3>::bse();
  j["bnw"] = Corner<3>::bnw();
  j["bne"] = Corner<3>::bne();
  j["tsw"] = Corner<3>::tsw();
  j["tse"] = Corner<3>::tse();
  j["tnw"] = Corner<3>::tnw();
  j["tne"] = Corner<3>::tne();
  CHECK(j["null"] == nullptr);
  CHECK(j["bsw"] == "BSW");
  CHECK(j["bse"] == "BSE");
  CHECK(j["bnw"] == "BNW");
  CHECK(j["bne"] == "BNE");
  CHECK(j["tsw"] == "TSW");
  CHECK(j["tse"] == "TSE");
  CHECK(j["tnw"] == "TNW");
  CHECK(j["tne"] == "TNE");
}