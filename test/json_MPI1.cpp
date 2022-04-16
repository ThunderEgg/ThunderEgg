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

#include <ThunderEgg/Domain.h>
#include <ThunderEgg/tpl/json.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

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
TEST_CASE("PatchInfo from_json with children", "[PatchInfo]")
{
  nlohmann::json j;
  j["id"] = 9;
  j["rank"] = 3;
  j["refine_level"] = 329;
  j["parent_id"] = 2;
  j["parent_rank"] = 3;
  j["orth_on_parent"] = "TNW";
  j["starts"] = { 1, 2, 3 };
  j["lengths"] = { 10, 20, 30 };
  j["child_ids"] = { 1, 2, 3, 4, 5, 6, 7, 8 };
  j["child_ranks"] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  j["nbrs"] = { NormalNbrInfo<2>(1),
                CoarseNbrInfo<2>(2, Orthant<2>::nw()),
                FineNbrInfo<2>({ 3, 4, 5, 6 }) };
  j["nbrs"][0]["side"] = "NORTH";
  j["nbrs"][1]["side"] = "EAST";
  j["nbrs"][2]["side"] = "SOUTH";

  PatchInfo<3> d = j.get<PatchInfo<3>>();
  CHECK(d.id == 9);
  CHECK(d.rank == 3);
  CHECK(d.refine_level == 329);
  CHECK(d.parent_id == 2);
  CHECK(d.parent_rank == 3);
  CHECK(d.orth_on_parent == Orthant<3>::tnw());
  CHECK(d.starts[0] == 1);
  CHECK(d.starts[1] == 2);
  CHECK(d.starts[2] == 3);
  CHECK(d.spacings[0] == 10);
  CHECK(d.spacings[1] == 20);
  CHECK(d.spacings[2] == 30);
  CHECK(d.ns[0] == 1);
  CHECK(d.ns[1] == 1);
  CHECK(d.ns[2] == 1);
  CHECK(d.child_ids[0] == 1);
  CHECK(d.child_ids[1] == 2);
  CHECK(d.child_ids[2] == 3);
  CHECK(d.child_ids[3] == 4);
  CHECK(d.child_ids[4] == 5);
  CHECK(d.child_ids[5] == 6);
  CHECK(d.child_ids[6] == 7);
  CHECK(d.child_ids[7] == 8);
  CHECK(d.child_ranks[0] == 0);
  CHECK(d.child_ranks[1] == 1);
  CHECK(d.child_ranks[2] == 2);
  CHECK(d.child_ranks[3] == 3);
  CHECK(d.child_ranks[4] == 4);
  CHECK(d.child_ranks[5] == 5);
  CHECK(d.child_ranks[6] == 6);
  CHECK(d.child_ranks[7] == 7);
  CHECK_FALSE(d.hasNbr(Side<3>::west()));
  CHECK(d.hasNbr(Side<3>::east()));
  CHECK(d.getNbrType(Side<3>::east()) == NbrType::Coarse);
  CHECK(d.hasNbr(Side<3>::south()));
  CHECK(d.getNbrType(Side<3>::south()) == NbrType::Fine);
  CHECK(d.hasNbr(Side<3>::north()));
  CHECK(d.getNbrType(Side<3>::north()) == NbrType::Normal);
  CHECK_FALSE(d.hasNbr(Side<3>::bottom()));
  CHECK_FALSE(d.hasNbr(Side<3>::top()));
}
TEST_CASE("PatchInfo to_json no children", "[PatchInfo]")
{
  PatchInfo<3> d;
  d.id = 9;
  d.rank = 0;
  d.parent_id = 2;
  d.parent_rank = 3;
  d.orth_on_parent = Orthant<3>::tnw();
  d.starts = { 1, 2, 3 };
  d.spacings = { 0.1, 0.2, 0.3 };
  d.ns = { 10, 20, 30 };
  d.setNbrInfo(Side<3>::north(), new NormalNbrInfo<2>(1));
  d.setNbrInfo(Side<3>::east(), new CoarseNbrInfo<2>(2, Orthant<2>::nw()));
  d.setNbrInfo(Side<3>::south(), new FineNbrInfo<2>({ 3, 4, 5, 6 }));
  d.setNbrInfo(Corner<3>::bsw(), new NormalNbrInfo<0>(1));
  d.setNbrInfo(Corner<3>::tse(), new CoarseNbrInfo<0>(2, Orthant<0>(0)));
  d.setNbrInfo(Corner<3>::bnw(), new FineNbrInfo<0>({ 1 }));
  d.setNbrInfo(Edge::sw(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Edge::bn(), new CoarseNbrInfo<1>(2, Orthant<1>::lower()));
  d.setNbrInfo(Edge::tw(), new FineNbrInfo<1>({ 1, 2 }));

  nlohmann::json j = d;

  CHECK(j["id"] == d.id);
  CHECK(j["parent_id"] == d.parent_id);
  CHECK(j["parent_rank"] == d.parent_rank);
  CHECK(j["orth_on_parent"] == "TNW");
  CHECK(j["rank"] == d.rank);
  CHECK(j["child_ids"] == nullptr);
  CHECK(j["child_ranks"] == nullptr);

  REQUIRE(j["starts"].is_array());
  REQUIRE(j["starts"].size() == 3);
  CHECK(j["starts"][0] == d.starts[0]);
  CHECK(j["starts"][1] == d.starts[1]);
  CHECK(j["starts"][2] == d.starts[2]);

  REQUIRE(j["lengths"].is_array());
  REQUIRE(j["lengths"].size() == 3);
  CHECK(j["lengths"][0] == d.spacings[0] * d.ns[0]);
  CHECK(j["lengths"][1] == d.spacings[1] * d.ns[1]);
  CHECK(j["lengths"][2] == d.spacings[2] * d.ns[2]);

  REQUIRE(j["nbrs"].is_array());
  REQUIRE(j["nbrs"].size() == 3);

  CHECK(j["nbrs"][0]["type"] == "COARSE");
  CHECK(j["nbrs"][0]["side"] == "EAST");

  CHECK(j["nbrs"][1]["type"] == "FINE");
  CHECK(j["nbrs"][1]["side"] == "SOUTH");

  CHECK(j["nbrs"][2]["type"] == "NORMAL");
  CHECK(j["nbrs"][2]["side"] == "NORTH");

  REQUIRE(j["corner_nbrs"].is_array());
  REQUIRE(j["corner_nbrs"].size() == 3);

  CHECK(j["corner_nbrs"][0]["type"] == "NORMAL");
  CHECK(j["corner_nbrs"][0]["corner"] == "BSW");

  CHECK(j["corner_nbrs"][1]["type"] == "FINE");
  CHECK(j["corner_nbrs"][1]["corner"] == "BNW");

  CHECK(j["corner_nbrs"][2]["type"] == "COARSE");
  CHECK(j["corner_nbrs"][2]["corner"] == "TSE");

  REQUIRE(j["edge_nbrs"].is_array());
  REQUIRE(j["edge_nbrs"].size() == 3);

  CHECK(j["edge_nbrs"][0]["type"] == "COARSE");
  CHECK(j["edge_nbrs"][0]["edge"] == "BN");

  CHECK(j["edge_nbrs"][1]["type"] == "FINE");
  CHECK(j["edge_nbrs"][1]["edge"] == "TW");

  CHECK(j["edge_nbrs"][2]["type"] == "NORMAL");
  CHECK(j["edge_nbrs"][2]["edge"] == "SW");
}
TEST_CASE("PatchInfo to_json no children no neighbors", "[PatchInfo]")
{
  PatchInfo<3> d;
  d.id = 9;
  d.rank = 0;
  d.parent_id = 2;
  d.parent_rank = 3;
  d.refine_level = 329;
  d.starts = { 1, 2, 3 };
  d.spacings = { 0.1, 0.2, 0.3 };
  d.ns = { 10, 20, 30 };

  nlohmann::json j = d;

  CHECK(j["id"] == d.id);
  CHECK(j["parent_id"] == d.parent_id);
  CHECK(j["parent_rank"] == d.parent_rank);
  CHECK(j["rank"] == d.rank);
  CHECK(j["refine_level"] == 329);
  CHECK(j["child_ids"] == nullptr);
  CHECK(j["child_ranks"] == nullptr);
  CHECK(j["orth_on_parent"] == nullptr);

  REQUIRE(j["starts"].is_array());
  REQUIRE(j["starts"].size() == 3);
  CHECK(j["starts"][0] == d.starts[0]);
  CHECK(j["starts"][1] == d.starts[1]);
  CHECK(j["starts"][2] == d.starts[2]);

  REQUIRE(j["lengths"].is_array());
  REQUIRE(j["lengths"].size() == 3);
  CHECK(j["lengths"][0] == d.spacings[0] * d.ns[0]);
  CHECK(j["lengths"][1] == d.spacings[1] * d.ns[1]);
  CHECK(j["lengths"][2] == d.spacings[2] * d.ns[2]);

  REQUIRE(j["nbrs"].is_array());
  REQUIRE(j["nbrs"].size() == 0);
}
TEST_CASE("PatchInfo to_json with children", "[PatchInfo]")
{
  PatchInfo<3> d;
  d.id = 9;
  d.rank = 0;
  d.parent_id = 2;
  d.parent_rank = 3;
  d.refine_level = 329;
  d.starts = { 1, 2, 3 };
  d.spacings = { 0.1, 0.2, 0.3 };
  d.ns = { 10, 20, 30 };
  d.child_ids = { 3, 4, 5, 6, 7, 8, 9, 10 };
  d.child_ranks = { 1, 2, 3, 4, 5, 6, 7, 8 };
  d.setNbrInfo(Side<3>::north(), new NormalNbrInfo<2>(1));
  d.setNbrInfo(Side<3>::east(), new CoarseNbrInfo<2>(2, Orthant<2>::nw()));
  d.setNbrInfo(Side<3>::south(), new FineNbrInfo<2>({ 3, 4, 5, 6 }));

  nlohmann::json j = d;

  CHECK(j["id"] == d.id);
  CHECK(j["parent_id"] == d.parent_id);
  CHECK(j["parent_rank"] == d.parent_rank);
  CHECK(j["rank"] == d.rank);
  CHECK(j["refine_level"] == 329);

  REQUIRE(j["child_ids"].is_array());
  REQUIRE(j["child_ids"].size() == 8);
  CHECK(j["child_ids"][0] == d.child_ids[0]);
  CHECK(j["child_ids"][1] == d.child_ids[1]);
  CHECK(j["child_ids"][2] == d.child_ids[2]);
  CHECK(j["child_ids"][3] == d.child_ids[3]);
  CHECK(j["child_ids"][4] == d.child_ids[4]);
  CHECK(j["child_ids"][5] == d.child_ids[5]);
  CHECK(j["child_ids"][6] == d.child_ids[6]);
  CHECK(j["child_ids"][7] == d.child_ids[7]);

  REQUIRE(j["child_ranks"].is_array());
  REQUIRE(j["child_ranks"].size() == 8);
  CHECK(j["child_ranks"][0] == d.child_ranks[0]);
  CHECK(j["child_ranks"][1] == d.child_ranks[1]);
  CHECK(j["child_ranks"][2] == d.child_ranks[2]);
  CHECK(j["child_ranks"][3] == d.child_ranks[3]);
  CHECK(j["child_ranks"][4] == d.child_ranks[4]);
  CHECK(j["child_ranks"][5] == d.child_ranks[5]);
  CHECK(j["child_ranks"][6] == d.child_ranks[6]);
  CHECK(j["child_ranks"][7] == d.child_ranks[7]);

  REQUIRE(j["starts"].is_array());
  REQUIRE(j["starts"].size() == 3);
  CHECK(j["starts"][0] == d.starts[0]);
  CHECK(j["starts"][1] == d.starts[1]);
  CHECK(j["starts"][2] == d.starts[2]);

  REQUIRE(j["lengths"].is_array());
  REQUIRE(j["lengths"].size() == 3);
  CHECK(j["lengths"][0] == d.spacings[0] * d.ns[0]);
  CHECK(j["lengths"][1] == d.spacings[1] * d.ns[1]);
  CHECK(j["lengths"][2] == d.spacings[2] * d.ns[2]);

  REQUIRE(j["nbrs"].is_array());
  REQUIRE(j["nbrs"].size() == 3);

  CHECK(j["nbrs"][0]["type"] == "COARSE");
  CHECK(j["nbrs"][0]["side"] == "EAST");

  CHECK(j["nbrs"][1]["type"] == "FINE");
  CHECK(j["nbrs"][1]["side"] == "SOUTH");

  CHECK(j["nbrs"][2]["type"] == "NORMAL");
  CHECK(j["nbrs"][2]["side"] == "NORTH");
}
TEST_CASE("PatchInfo from_json no children", "[PatchInfo]")
{
  nlohmann::json j;
  j["id"] = 9;
  j["rank"] = 3;
  j["refine_level"] = 329;
  j["parent_id"] = 2;
  j["parent_rank"] = 3;
  j["starts"] = { 1, 2, 3 };
  j["lengths"] = { 10, 20, 30 };
  j["nbrs"] = { NormalNbrInfo<2>(1),
                CoarseNbrInfo<2>(2, Orthant<2>::nw()),
                FineNbrInfo<2>({ 3, 4, 5, 6 }) };
  j["nbrs"][0]["side"] = "NORTH";
  j["nbrs"][1]["side"] = "EAST";
  j["nbrs"][2]["side"] = "SOUTH";
  j["corner_nbrs"] = { NormalNbrInfo<0>(1),
                       CoarseNbrInfo<0>(2, Orthant<0>(0)),
                       FineNbrInfo<0>({ 1 }) };
  j["corner_nbrs"][0]["corner"] = "BSW";
  j["corner_nbrs"][1]["corner"] = "TSE";
  j["corner_nbrs"][2]["corner"] = "BNW";
  j["edge_nbrs"] = { NormalNbrInfo<1>(1),
                     CoarseNbrInfo<1>(2, Orthant<1>::lower()),
                     FineNbrInfo<1>({ 1, 2 }) };
  j["edge_nbrs"][0]["edge"] = "SW";
  j["edge_nbrs"][1]["edge"] = "BN";
  j["edge_nbrs"][2]["edge"] = "TW";

  PatchInfo<3> d = j.get<PatchInfo<3>>();
  CHECK(d.id == 9);
  CHECK(d.rank == 3);
  CHECK(d.refine_level == 329);
  CHECK(d.parent_id == 2);
  CHECK(d.parent_rank == 3);
  CHECK(d.orth_on_parent == Orthant<3>::null());
  CHECK(d.starts[0] == 1);
  CHECK(d.starts[1] == 2);
  CHECK(d.starts[2] == 3);
  CHECK(d.spacings[0] == 10);
  CHECK(d.spacings[1] == 20);
  CHECK(d.spacings[2] == 30);
  CHECK(d.ns[0] == 1);
  CHECK(d.ns[1] == 1);
  CHECK(d.ns[2] == 1);
  CHECK_FALSE(d.hasNbr(Side<3>::west()));
  CHECK(d.hasNbr(Side<3>::east()));
  CHECK(d.getNbrType(Side<3>::east()) == NbrType::Coarse);
  CHECK(d.hasNbr(Side<3>::south()));
  CHECK(d.getNbrType(Side<3>::south()) == NbrType::Fine);
  CHECK(d.hasNbr(Side<3>::north()));
  CHECK(d.getNbrType(Side<3>::north()) == NbrType::Normal);
  CHECK_FALSE(d.hasNbr(Side<3>::bottom()));
  CHECK_FALSE(d.hasNbr(Side<3>::top()));

  CHECK(d.hasNbr(Corner<3>::bsw()));
  CHECK(d.getNbrType(Corner<3>::bsw()) == NbrType::Normal);
  CHECK_FALSE(d.hasNbr(Corner<3>::bse()));
  CHECK(d.hasNbr(Corner<3>::bnw()));
  CHECK(d.getNbrType(Corner<3>::bnw()) == NbrType::Fine);
  CHECK_FALSE(d.hasNbr(Corner<3>::bne()));
  CHECK_FALSE(d.hasNbr(Corner<3>::tsw()));
  CHECK(d.hasNbr(Corner<3>::tse()));
  CHECK(d.getNbrType(Corner<3>::tse()) == NbrType::Coarse);
  CHECK_FALSE(d.hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(d.hasNbr(Corner<3>::tne()));

  CHECK_FALSE(d.hasNbr(Edge::bs()));
  CHECK_FALSE(d.hasNbr(Edge::tn()));
  CHECK(d.hasNbr(Edge::bn()));
  CHECK(d.getNbrType(Edge::bn()) == NbrType::Coarse);
  CHECK_FALSE(d.hasNbr(Edge::ts()));
  CHECK_FALSE(d.hasNbr(Edge::bw()));
  CHECK_FALSE(d.hasNbr(Edge::te()));
  CHECK_FALSE(d.hasNbr(Edge::be()));
  CHECK(d.hasNbr(Edge::tw()));
  CHECK(d.getNbrType(Edge::tw()) == NbrType::Fine);
  CHECK(d.hasNbr(Edge::sw()));
  CHECK(d.getNbrType(Edge::sw()) == NbrType::Normal);
  CHECK_FALSE(d.hasNbr(Edge::ne()));
  CHECK_FALSE(d.hasNbr(Edge::se()));
  CHECK_FALSE(d.hasNbr(Edge::nw()));
}
TEST_CASE("FineNbrInfo to_json", "[FineNbrInfo]")
{
  FineNbrInfo<2> info;
  info.ids[0] = GENERATE(1, 2);
  info.ids[1] = GENERATE(1, 2);
  info.ids[2] = GENERATE(1, 2);
  info.ids[3] = GENERATE(1, 2);
  info.ranks[0] = GENERATE(0, 1);
  info.ranks[1] = GENERATE(0, 1);
  info.ranks[2] = GENERATE(0, 1);
  info.ranks[3] = GENERATE(0, 1);

  nlohmann::json j = info;

  CHECK(j["type"] == "FINE");
  REQUIRE(j["ids"].is_array());
  REQUIRE(j["ids"].size() == 4);
  CHECK(j["ids"][0] == info.ids[0]);
  CHECK(j["ids"][1] == info.ids[1]);
  CHECK(j["ids"][2] == info.ids[2]);
  CHECK(j["ids"][3] == info.ids[3]);
  REQUIRE(j["ranks"].is_array());
  REQUIRE(j["ranks"].size() == 4);
  CHECK(j["ranks"][0] == info.ranks[0]);
  CHECK(j["ranks"][1] == info.ranks[1]);
  CHECK(j["ranks"][2] == info.ranks[2]);
  CHECK(j["ranks"][3] == info.ranks[3]);
}
TEST_CASE("FineNbrInfo from_json", "[FineNbrInfo]")
{
  int id1 = GENERATE(1, 2);
  int id2 = GENERATE(1, 2);
  int id3 = GENERATE(1, 2);
  int id4 = GENERATE(1, 2);
  int rank1 = GENERATE(0, 1);
  int rank2 = GENERATE(0, 1);
  int rank3 = GENERATE(0, 1);
  int rank4 = GENERATE(0, 1);

  nlohmann::json j;
  j["type"] = "NORMAL";
  j["ids"] = { id1, id2, id3, id4 };
  j["ranks"] = { rank1, rank2, rank3, rank4 };

  FineNbrInfo<2> info = j.get<FineNbrInfo<2>>();
  CHECK(info.ids[0] == id1);
  CHECK(info.ids[1] == id2);
  CHECK(info.ids[2] == id3);
  CHECK(info.ids[3] == id4);
  CHECK(info.ranks[0] == rank1);
  CHECK(info.ranks[1] == rank2);
  CHECK(info.ranks[2] == rank3);
  CHECK(info.ranks[3] == rank4);
}
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
TEST_CASE("NormalNbrInfo to_json", "[NormalNbrInfo]")
{
  NormalNbrInfo<2> info;
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

  NormalNbrInfo<2> info = j.get<NormalNbrInfo<2>>();
  CHECK(info.id == id);
  CHECK(info.rank == rank);
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