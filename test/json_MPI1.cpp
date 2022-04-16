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

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::tpl;

TEST_CASE("Test from_json for Edge")
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
  CHECK_EQ(j["null"].get<Edge>(), Edge::null());
  CHECK_EQ(j["bs"].get<Edge>(), Edge::bs());
  CHECK_EQ(j["tn"].get<Edge>(), Edge::tn());
  CHECK_EQ(j["bn"].get<Edge>(), Edge::bn());
  CHECK_EQ(j["ts"].get<Edge>(), Edge::ts());
  CHECK_EQ(j["bw"].get<Edge>(), Edge::bw());
  CHECK_EQ(j["te"].get<Edge>(), Edge::te());
  CHECK_EQ(j["be"].get<Edge>(), Edge::be());
  CHECK_EQ(j["tw"].get<Edge>(), Edge::tw());
  CHECK_EQ(j["sw"].get<Edge>(), Edge::sw());
  CHECK_EQ(j["ne"].get<Edge>(), Edge::ne());
  CHECK_EQ(j["nw"].get<Edge>(), Edge::nw());
}

TEST_CASE("Test to_json for Edge")
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
  CHECK_EQ(j["null"], nullptr);
  CHECK_EQ(j["bs"], "BS");
  CHECK_EQ(j["tn"], "TN");
  CHECK_EQ(j["bn"], "BN");
  CHECK_EQ(j["ts"], "TS");
  CHECK_EQ(j["bw"], "BW");
  CHECK_EQ(j["te"], "TE");
  CHECK_EQ(j["be"], "BE");
  CHECK_EQ(j["tw"], "TW");
  CHECK_EQ(j["sw"], "SW");
  CHECK_EQ(j["ne"], "NE");
  CHECK_EQ(j["se"], "SE");
  CHECK_EQ(j["nw"], "NW");
}

TEST_CASE("Test to_json for Side<1>")
{
  nlohmann::json j;
  j["null"] = Side<1>::null();
  j["west"] = Side<1>::west();
  j["east"] = Side<1>::east();
  CHECK_EQ(j["null"], nullptr);
  CHECK_EQ(j["west"], "WEST");
  CHECK_EQ(j["east"], "EAST");
}

TEST_CASE("Test to_json for Side<2>")
{
  nlohmann::json j;
  j["null"] = Side<2>::null();
  j["west"] = Side<2>::west();
  j["east"] = Side<2>::east();
  j["south"] = Side<2>::south();
  j["north"] = Side<2>::north();
  CHECK_EQ(j["null"], nullptr);
  CHECK_EQ(j["west"], "WEST");
  CHECK_EQ(j["east"], "EAST");
  CHECK_EQ(j["south"], "SOUTH");
  CHECK_EQ(j["north"], "NORTH");
}

TEST_CASE("Test to_json for Side<3>")
{
  nlohmann::json j;
  j["null"] = Side<3>::null();
  j["west"] = Side<3>::west();
  j["east"] = Side<3>::east();
  j["south"] = Side<3>::south();
  j["north"] = Side<3>::north();
  j["bottom"] = Side<3>::bottom();
  j["top"] = Side<3>::top();
  CHECK_EQ(j["null"], nullptr);
  CHECK_EQ(j["west"], "WEST");
  CHECK_EQ(j["east"], "EAST");
  CHECK_EQ(j["south"], "SOUTH");
  CHECK_EQ(j["north"], "NORTH");
  CHECK_EQ(j["bottom"], "BOTTOM");
  CHECK_EQ(j["top"], "TOP");
}

TEST_CASE("Test from_json for Side<1>")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["west"] = "WEST";
  j["east"] = "EAST";
  CHECK_EQ(j["null"].get<Side<1>>(), Side<1>::null());
  CHECK_EQ(j["west"].get<Side<1>>(), Side<1>::west());
  CHECK_EQ(j["east"].get<Side<1>>(), Side<1>::east());
}

TEST_CASE("Test from_json for Side<2>")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["west"] = "WEST";
  j["east"] = "EAST";
  j["south"] = "SOUTH";
  j["north"] = "NORTH";
  CHECK_EQ(j["null"].get<Side<2>>(), Side<2>::null());
  CHECK_EQ(j["west"].get<Side<2>>(), Side<2>::west());
  CHECK_EQ(j["east"].get<Side<3>>(), Side<3>::east());
  CHECK_EQ(j["south"].get<Side<2>>(), Side<2>::south());
  CHECK_EQ(j["north"].get<Side<2>>(), Side<2>::north());
}

TEST_CASE("Test from_json for Side<3>")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["west"] = "WEST";
  j["east"] = "EAST";
  j["south"] = "SOUTH";
  j["north"] = "NORTH";
  j["bottom"] = "BOTTOM";
  j["top"] = "TOP";
  CHECK_EQ(j["null"].get<Side<3>>(), Side<3>::null());
  CHECK_EQ(j["west"].get<Side<3>>(), Side<3>::west());
  CHECK_EQ(j["east"].get<Side<3>>(), Side<3>::east());
  CHECK_EQ(j["south"].get<Side<3>>(), Side<3>::south());
  CHECK_EQ(j["north"].get<Side<3>>(), Side<3>::north());
  CHECK_EQ(j["bottom"].get<Side<3>>(), Side<3>::bottom());
  CHECK_EQ(j["top"].get<Side<3>>(), Side<3>::top());
}

TEST_CASE("Test from_json for Orthant<0>")
{
  nlohmann::json j;
  j["null"] = nullptr;
  CHECK_EQ(j["null"].get<Orthant<0>>(), Orthant<0>::null());
}

TEST_CASE("Test from_json for Orthant<1>")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["lower"] = "LOWER";
  j["upper"] = "UPPER";
  CHECK_EQ(j["null"].get<Orthant<1>>(), Orthant<1>::null());
  CHECK_EQ(j["lower"].get<Orthant<1>>(), Orthant<1>::lower());
  CHECK_EQ(j["upper"].get<Orthant<1>>(), Orthant<1>::upper());
}

TEST_CASE("Test from_json for Orthant<2>")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["sw"] = "SW";
  j["se"] = "SE";
  j["nw"] = "NW";
  j["ne"] = "NE";
  CHECK_EQ(j["null"].get<Orthant<2>>(), Orthant<2>::null());
  CHECK_EQ(j["sw"].get<Orthant<2>>(), Orthant<2>::sw());
  CHECK_EQ(j["se"].get<Orthant<2>>(), Orthant<2>::se());
  CHECK_EQ(j["nw"].get<Orthant<2>>(), Orthant<2>::nw());
  CHECK_EQ(j["ne"].get<Orthant<2>>(), Orthant<2>::ne());
}

TEST_CASE("Test from_json for Orthant<3>")
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
  CHECK_EQ(j["null"].get<Orthant<3>>(), Orthant<3>::null());
  CHECK_EQ(j["bsw"].get<Orthant<3>>(), Orthant<3>::bsw());
  CHECK_EQ(j["bse"].get<Orthant<3>>(), Orthant<3>::bse());
  CHECK_EQ(j["bnw"].get<Orthant<3>>(), Orthant<3>::bnw());
  CHECK_EQ(j["bne"].get<Orthant<3>>(), Orthant<3>::bne());
  CHECK_EQ(j["tsw"].get<Orthant<3>>(), Orthant<3>::tsw());
  CHECK_EQ(j["tse"].get<Orthant<3>>(), Orthant<3>::tse());
  CHECK_EQ(j["tnw"].get<Orthant<3>>(), Orthant<3>::tnw());
  CHECK_EQ(j["tne"].get<Orthant<3>>(), Orthant<3>::tne());
}

TEST_CASE("Test to_json for Orthant<0>")
{
  nlohmann::json j;
  j["null"] = Orthant<0>::null();
  CHECK_EQ(j["null"], nullptr);
}

TEST_CASE("Test to_json for Orthant<1>")
{
  nlohmann::json j;
  j["null"] = Orthant<1>::null();
  j["lower"] = Orthant<1>::lower();
  j["upper"] = Orthant<1>::upper();
  CHECK_EQ(j["null"], nullptr);
  CHECK_EQ(j["lower"], "LOWER");
  CHECK_EQ(j["upper"], "UPPER");
}

TEST_CASE("Test to_json for Orthant<2>")
{
  nlohmann::json j;
  j["null"] = Orthant<2>::null();
  j["sw"] = Orthant<2>::sw();
  j["se"] = Orthant<2>::se();
  j["nw"] = Orthant<2>::nw();
  j["ne"] = Orthant<2>::ne();
  CHECK_EQ(j["null"], nullptr);
  CHECK_EQ(j["sw"], "SW");
  CHECK_EQ(j["se"], "SE");
  CHECK_EQ(j["nw"], "NW");
  CHECK_EQ(j["ne"], "NE");
}

TEST_CASE("Test to_json for Orthant<3>")
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
  CHECK_EQ(j["null"], nullptr);
  CHECK_EQ(j["bsw"], "BSW");
  CHECK_EQ(j["bse"], "BSE");
  CHECK_EQ(j["bnw"], "BNW");
  CHECK_EQ(j["bne"], "BNE");
  CHECK_EQ(j["tsw"], "TSW");
  CHECK_EQ(j["tse"], "TSE");
  CHECK_EQ(j["tnw"], "TNW");
  CHECK_EQ(j["tne"], "TNE");
}

TEST_CASE("Test from_json for Corner<2>")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["sw"] = "SW";
  j["se"] = "SE";
  j["nw"] = "NW";
  j["ne"] = "NE";
  CHECK_EQ(j["null"].get<Corner<2>>(), Corner<2>::null());
  CHECK_EQ(j["sw"].get<Corner<2>>(), Corner<2>::sw());
  CHECK_EQ(j["se"].get<Corner<2>>(), Corner<2>::se());
  CHECK_EQ(j["nw"].get<Corner<2>>(), Corner<2>::nw());
  CHECK_EQ(j["ne"].get<Corner<2>>(), Corner<2>::ne());
}

TEST_CASE("Test from_json for Corner<3>")
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
  CHECK_EQ(j["null"].get<Corner<3>>(), Corner<3>::null());
  CHECK_EQ(j["bsw"].get<Corner<3>>(), Corner<3>::bsw());
  CHECK_EQ(j["bse"].get<Corner<3>>(), Corner<3>::bse());
  CHECK_EQ(j["bnw"].get<Corner<3>>(), Corner<3>::bnw());
  CHECK_EQ(j["bne"].get<Corner<3>>(), Corner<3>::bne());
  CHECK_EQ(j["tsw"].get<Corner<3>>(), Corner<3>::tsw());
  CHECK_EQ(j["tse"].get<Corner<3>>(), Corner<3>::tse());
  CHECK_EQ(j["tnw"].get<Corner<3>>(), Corner<3>::tnw());
  CHECK_EQ(j["tne"].get<Corner<3>>(), Corner<3>::tne());
}

TEST_CASE("Test to_json for Corner<2>")
{
  nlohmann::json j;
  j["null"] = Corner<2>::null();
  j["sw"] = Corner<2>::sw();
  j["se"] = Corner<2>::se();
  j["nw"] = Corner<2>::nw();
  j["ne"] = Corner<2>::ne();
  CHECK_EQ(j["null"], nullptr);
  CHECK_EQ(j["sw"], "SW");
  CHECK_EQ(j["se"], "SE");
  CHECK_EQ(j["nw"], "NW");
  CHECK_EQ(j["ne"], "NE");
}

TEST_CASE("Test to_json for Corner<3>")
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
  CHECK_EQ(j["null"], nullptr);
  CHECK_EQ(j["bsw"], "BSW");
  CHECK_EQ(j["bse"], "BSE");
  CHECK_EQ(j["bnw"], "BNW");
  CHECK_EQ(j["bne"], "BNE");
  CHECK_EQ(j["tsw"], "TSW");
  CHECK_EQ(j["tse"], "TSE");
  CHECK_EQ(j["tnw"], "TNW");
  CHECK_EQ(j["tne"], "TNE");
}

TEST_CASE("PatchInfo from_json with children")
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
  CHECK_EQ(d.id, 9);
  CHECK_EQ(d.rank, 3);
  CHECK_EQ(d.refine_level, 329);
  CHECK_EQ(d.parent_id, 2);
  CHECK_EQ(d.parent_rank, 3);
  CHECK_EQ(d.orth_on_parent, Orthant<3>::tnw());
  CHECK_EQ(d.starts[0], 1);
  CHECK_EQ(d.starts[1], 2);
  CHECK_EQ(d.starts[2], 3);
  CHECK_EQ(d.spacings[0], 10);
  CHECK_EQ(d.spacings[1], 20);
  CHECK_EQ(d.spacings[2], 30);
  CHECK_EQ(d.ns[0], 1);
  CHECK_EQ(d.ns[1], 1);
  CHECK_EQ(d.ns[2], 1);
  CHECK_EQ(d.child_ids[0], 1);
  CHECK_EQ(d.child_ids[1], 2);
  CHECK_EQ(d.child_ids[2], 3);
  CHECK_EQ(d.child_ids[3], 4);
  CHECK_EQ(d.child_ids[4], 5);
  CHECK_EQ(d.child_ids[5], 6);
  CHECK_EQ(d.child_ids[6], 7);
  CHECK_EQ(d.child_ids[7], 8);
  CHECK_EQ(d.child_ranks[0], 0);
  CHECK_EQ(d.child_ranks[1], 1);
  CHECK_EQ(d.child_ranks[2], 2);
  CHECK_EQ(d.child_ranks[3], 3);
  CHECK_EQ(d.child_ranks[4], 4);
  CHECK_EQ(d.child_ranks[5], 5);
  CHECK_EQ(d.child_ranks[6], 6);
  CHECK_EQ(d.child_ranks[7], 7);
  CHECK_UNARY_FALSE(d.hasNbr(Side<3>::west()));
  CHECK_UNARY(d.hasNbr(Side<3>::east()));
  CHECK_EQ(d.getNbrType(Side<3>::east()), NbrType::Coarse);
  CHECK_UNARY(d.hasNbr(Side<3>::south()));
  CHECK_EQ(d.getNbrType(Side<3>::south()), NbrType::Fine);
  CHECK_UNARY(d.hasNbr(Side<3>::north()));
  CHECK_EQ(d.getNbrType(Side<3>::north()), NbrType::Normal);
  CHECK_UNARY_FALSE(d.hasNbr(Side<3>::bottom()));
  CHECK_UNARY_FALSE(d.hasNbr(Side<3>::top()));
}

TEST_CASE("PatchInfo to_json no children")
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

  CHECK_EQ(j["id"], d.id);
  CHECK_EQ(j["parent_id"], d.parent_id);
  CHECK_EQ(j["parent_rank"], d.parent_rank);
  CHECK_EQ(j["orth_on_parent"], "TNW");
  CHECK_EQ(j["rank"], d.rank);
  CHECK_EQ(j["child_ids"], nullptr);
  CHECK_EQ(j["child_ranks"], nullptr);

  REQUIRE_UNARY(j["starts"].is_array());
  REQUIRE_EQ(j["starts"].size(), 3);
  CHECK_EQ(j["starts"][0], d.starts[0]);
  CHECK_EQ(j["starts"][1], d.starts[1]);
  CHECK_EQ(j["starts"][2], d.starts[2]);

  REQUIRE_UNARY(j["lengths"].is_array());
  REQUIRE_EQ(j["lengths"].size(), 3);
  CHECK_EQ(j["lengths"][0], d.spacings[0] * d.ns[0]);
  CHECK_EQ(j["lengths"][1], d.spacings[1] * d.ns[1]);
  CHECK_EQ(j["lengths"][2], d.spacings[2] * d.ns[2]);

  REQUIRE_UNARY(j["nbrs"].is_array());
  REQUIRE_EQ(j["nbrs"].size(), 3);

  CHECK_EQ(j["nbrs"][0]["type"], "COARSE");
  CHECK_EQ(j["nbrs"][0]["side"], "EAST");

  CHECK_EQ(j["nbrs"][1]["type"], "FINE");
  CHECK_EQ(j["nbrs"][1]["side"], "SOUTH");

  CHECK_EQ(j["nbrs"][2]["type"], "NORMAL");
  CHECK_EQ(j["nbrs"][2]["side"], "NORTH");

  REQUIRE_UNARY(j["corner_nbrs"].is_array());
  REQUIRE_EQ(j["corner_nbrs"].size(), 3);

  CHECK_EQ(j["corner_nbrs"][0]["type"], "NORMAL");
  CHECK_EQ(j["corner_nbrs"][0]["corner"], "BSW");

  CHECK_EQ(j["corner_nbrs"][1]["type"], "FINE");
  CHECK_EQ(j["corner_nbrs"][1]["corner"], "BNW");

  CHECK_EQ(j["corner_nbrs"][2]["type"], "COARSE");
  CHECK_EQ(j["corner_nbrs"][2]["corner"], "TSE");

  REQUIRE_UNARY(j["edge_nbrs"].is_array());
  REQUIRE_EQ(j["edge_nbrs"].size(), 3);

  CHECK_EQ(j["edge_nbrs"][0]["type"], "COARSE");
  CHECK_EQ(j["edge_nbrs"][0]["edge"], "BN");

  CHECK_EQ(j["edge_nbrs"][1]["type"], "FINE");
  CHECK_EQ(j["edge_nbrs"][1]["edge"], "TW");

  CHECK_EQ(j["edge_nbrs"][2]["type"], "NORMAL");
  CHECK_EQ(j["edge_nbrs"][2]["edge"], "SW");
}

TEST_CASE("PatchInfo to_json no children no neighbors")
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

  CHECK_EQ(j["id"], d.id);
  CHECK_EQ(j["parent_id"], d.parent_id);
  CHECK_EQ(j["parent_rank"], d.parent_rank);
  CHECK_EQ(j["rank"], d.rank);
  CHECK_EQ(j["refine_level"], 329);
  CHECK_EQ(j["child_ids"], nullptr);
  CHECK_EQ(j["child_ranks"], nullptr);
  CHECK_EQ(j["orth_on_parent"], nullptr);

  REQUIRE_UNARY(j["starts"].is_array());
  REQUIRE_EQ(j["starts"].size(), 3);
  CHECK_EQ(j["starts"][0], d.starts[0]);
  CHECK_EQ(j["starts"][1], d.starts[1]);
  CHECK_EQ(j["starts"][2], d.starts[2]);

  REQUIRE_UNARY(j["lengths"].is_array());
  REQUIRE_EQ(j["lengths"].size(), 3);
  CHECK_EQ(j["lengths"][0], d.spacings[0] * d.ns[0]);
  CHECK_EQ(j["lengths"][1], d.spacings[1] * d.ns[1]);
  CHECK_EQ(j["lengths"][2], d.spacings[2] * d.ns[2]);

  REQUIRE_UNARY(j["nbrs"].is_array());
  REQUIRE_EQ(j["nbrs"].size(), 0);
}

TEST_CASE("PatchInfo to_json with children")
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

  CHECK_EQ(j["id"], d.id);
  CHECK_EQ(j["parent_id"], d.parent_id);
  CHECK_EQ(j["parent_rank"], d.parent_rank);
  CHECK_EQ(j["rank"], d.rank);
  CHECK_EQ(j["refine_level"], 329);

  REQUIRE_UNARY(j["child_ids"].is_array());
  REQUIRE_EQ(j["child_ids"].size(), 8);
  CHECK_EQ(j["child_ids"][0], d.child_ids[0]);
  CHECK_EQ(j["child_ids"][1], d.child_ids[1]);
  CHECK_EQ(j["child_ids"][2], d.child_ids[2]);
  CHECK_EQ(j["child_ids"][3], d.child_ids[3]);
  CHECK_EQ(j["child_ids"][4], d.child_ids[4]);
  CHECK_EQ(j["child_ids"][5], d.child_ids[5]);
  CHECK_EQ(j["child_ids"][6], d.child_ids[6]);
  CHECK_EQ(j["child_ids"][7], d.child_ids[7]);

  REQUIRE_UNARY(j["child_ranks"].is_array());
  REQUIRE_EQ(j["child_ranks"].size(), 8);
  CHECK_EQ(j["child_ranks"][0], d.child_ranks[0]);
  CHECK_EQ(j["child_ranks"][1], d.child_ranks[1]);
  CHECK_EQ(j["child_ranks"][2], d.child_ranks[2]);
  CHECK_EQ(j["child_ranks"][3], d.child_ranks[3]);
  CHECK_EQ(j["child_ranks"][4], d.child_ranks[4]);
  CHECK_EQ(j["child_ranks"][5], d.child_ranks[5]);
  CHECK_EQ(j["child_ranks"][6], d.child_ranks[6]);
  CHECK_EQ(j["child_ranks"][7], d.child_ranks[7]);

  REQUIRE_UNARY(j["starts"].is_array());
  REQUIRE_EQ(j["starts"].size(), 3);
  CHECK_EQ(j["starts"][0], d.starts[0]);
  CHECK_EQ(j["starts"][1], d.starts[1]);
  CHECK_EQ(j["starts"][2], d.starts[2]);

  REQUIRE_UNARY(j["lengths"].is_array());
  REQUIRE_EQ(j["lengths"].size(), 3);
  CHECK_EQ(j["lengths"][0], d.spacings[0] * d.ns[0]);
  CHECK_EQ(j["lengths"][1], d.spacings[1] * d.ns[1]);
  CHECK_EQ(j["lengths"][2], d.spacings[2] * d.ns[2]);

  REQUIRE_UNARY(j["nbrs"].is_array());
  REQUIRE_EQ(j["nbrs"].size(), 3);

  CHECK_EQ(j["nbrs"][0]["type"], "COARSE");
  CHECK_EQ(j["nbrs"][0]["side"], "EAST");

  CHECK_EQ(j["nbrs"][1]["type"], "FINE");
  CHECK_EQ(j["nbrs"][1]["side"], "SOUTH");

  CHECK_EQ(j["nbrs"][2]["type"], "NORMAL");
  CHECK_EQ(j["nbrs"][2]["side"], "NORTH");
}

TEST_CASE("PatchInfo from_json no children")
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
  CHECK_EQ(d.id, 9);
  CHECK_EQ(d.rank, 3);
  CHECK_EQ(d.refine_level, 329);
  CHECK_EQ(d.parent_id, 2);
  CHECK_EQ(d.parent_rank, 3);
  CHECK_EQ(d.orth_on_parent, Orthant<3>::null());
  CHECK_EQ(d.starts[0], 1);
  CHECK_EQ(d.starts[1], 2);
  CHECK_EQ(d.starts[2], 3);
  CHECK_EQ(d.spacings[0], 10);
  CHECK_EQ(d.spacings[1], 20);
  CHECK_EQ(d.spacings[2], 30);
  CHECK_EQ(d.ns[0], 1);
  CHECK_EQ(d.ns[1], 1);
  CHECK_EQ(d.ns[2], 1);
  CHECK_UNARY_FALSE(d.hasNbr(Side<3>::west()));
  CHECK_UNARY(d.hasNbr(Side<3>::east()));
  CHECK_EQ(d.getNbrType(Side<3>::east()), NbrType::Coarse);
  CHECK_UNARY(d.hasNbr(Side<3>::south()));
  CHECK_EQ(d.getNbrType(Side<3>::south()), NbrType::Fine);
  CHECK_UNARY(d.hasNbr(Side<3>::north()));
  CHECK_EQ(d.getNbrType(Side<3>::north()), NbrType::Normal);
  CHECK_UNARY_FALSE(d.hasNbr(Side<3>::bottom()));
  CHECK_UNARY_FALSE(d.hasNbr(Side<3>::top()));

  CHECK_UNARY(d.hasNbr(Corner<3>::bsw()));
  CHECK_EQ(d.getNbrType(Corner<3>::bsw()), NbrType::Normal);
  CHECK_UNARY_FALSE(d.hasNbr(Corner<3>::bse()));
  CHECK_UNARY(d.hasNbr(Corner<3>::bnw()));
  CHECK_EQ(d.getNbrType(Corner<3>::bnw()), NbrType::Fine);
  CHECK_UNARY_FALSE(d.hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(d.hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(d.hasNbr(Corner<3>::tse()));
  CHECK_EQ(d.getNbrType(Corner<3>::tse()), NbrType::Coarse);
  CHECK_UNARY_FALSE(d.hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(d.hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(d.hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(d.hasNbr(Edge::tn()));
  CHECK_UNARY(d.hasNbr(Edge::bn()));
  CHECK_EQ(d.getNbrType(Edge::bn()), NbrType::Coarse);
  CHECK_UNARY_FALSE(d.hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(d.hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(d.hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(d.hasNbr(Edge::be()));
  CHECK_UNARY(d.hasNbr(Edge::tw()));
  CHECK_EQ(d.getNbrType(Edge::tw()), NbrType::Fine);
  CHECK_UNARY(d.hasNbr(Edge::sw()));
  CHECK_EQ(d.getNbrType(Edge::sw()), NbrType::Normal);
  CHECK_UNARY_FALSE(d.hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(d.hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(d.hasNbr(Edge::nw()));
}

TEST_CASE("FineNbrInfo to_json")
{
  for (int ids_0 : { 1, 2 }) {
    for (int ids_1 : { 1, 2 }) {
      for (int ids_2 : { 1, 2 }) {
        for (int ids_3 : { 1, 2 }) {
          for (int ranks_0 : { 0, 1 }) {
            for (int ranks_1 : { 0, 1 }) {
              for (int ranks_2 : { 0, 1 }) {
                for (int ranks_3 : { 0, 1 }) {
                  FineNbrInfo<2> info;
                  info.ids[0] = ids_0;
                  info.ids[1] = ids_1;
                  info.ids[2] = ids_2;
                  info.ids[3] = ids_3;
                  info.ranks[0] = ranks_0;
                  info.ranks[1] = ranks_1;
                  info.ranks[2] = ranks_2;
                  info.ranks[3] = ranks_3;

                  nlohmann::json j = info;

                  CHECK_EQ(j["type"], "FINE");
                  REQUIRE_UNARY(j["ids"].is_array());
                  REQUIRE_EQ(j["ids"].size(), 4);
                  CHECK_EQ(j["ids"][0], info.ids[0]);
                  CHECK_EQ(j["ids"][1], info.ids[1]);
                  CHECK_EQ(j["ids"][2], info.ids[2]);
                  CHECK_EQ(j["ids"][3], info.ids[3]);
                  REQUIRE_UNARY(j["ranks"].is_array());
                  REQUIRE_EQ(j["ranks"].size(), 4);
                  CHECK_EQ(j["ranks"][0], info.ranks[0]);
                  CHECK_EQ(j["ranks"][1], info.ranks[1]);
                  CHECK_EQ(j["ranks"][2], info.ranks[2]);
                  CHECK_EQ(j["ranks"][3], info.ranks[3]);
                }
              }
            }
          }
        }
      }
    }
  }
}

TEST_CASE("FineNbrInfo from_json")
{
  for (int id1 : { 1, 2 }) {
    for (int id2 : { 1, 2 }) {
      for (int id3 : { 1, 2 }) {
        for (int id4 : { 1, 2 }) {
          for (int rank1 : { 0, 1 }) {
            for (int rank2 : { 0, 1 }) {
              for (int rank3 : { 0, 1 }) {
                for (int rank4 : { 0, 1 }) {

                  nlohmann::json j;
                  j["type"] = "NORMAL";
                  j["ids"] = { id1, id2, id3, id4 };
                  j["ranks"] = { rank1, rank2, rank3, rank4 };

                  FineNbrInfo<2> info = j.get<FineNbrInfo<2>>();
                  CHECK_EQ(info.ids[0], id1);
                  CHECK_EQ(info.ids[1], id2);
                  CHECK_EQ(info.ids[2], id3);
                  CHECK_EQ(info.ids[3], id4);
                  CHECK_EQ(info.ranks[0], rank1);
                  CHECK_EQ(info.ranks[1], rank2);
                  CHECK_EQ(info.ranks[2], rank3);
                  CHECK_EQ(info.ranks[3], rank4);
                }
              }
            }
          }
        }
      }
    }
  }
}

TEST_CASE("NbrType to_json")
{
  nlohmann::json j;
  j["normal"] = NbrType::Normal;
  j["coarse"] = NbrType::Coarse;
  j["fine"] = NbrType::Fine;
  CHECK_EQ(j["normal"], "NORMAL");
  CHECK_EQ(j["coarse"], "COARSE");
  CHECK_EQ(j["fine"], "FINE");
}

TEST_CASE("NbrType from_json")
{
  nlohmann::json j;
  j["normal"] = "NORMAL";
  j["coarse"] = "COARSE";
  j["fine"] = "FINE";
  CHECK_EQ(j["normal"].get<NbrType>(), NbrType::Normal);
  CHECK_EQ(j["coarse"].get<NbrType>(), NbrType::Coarse);
  CHECK_EQ(j["fine"].get<NbrType>(), NbrType::Fine);
}

TEST_CASE("NormalNbrInfo to_json")
{
  for (int id : { 1, 2, 3 }) {
    for (int rank : { 0, 1, 2 }) {
      NormalNbrInfo<2> info;
      info.id = id;
      info.rank = rank;

      nlohmann::json j = info;

      CHECK_EQ(j["type"], "NORMAL");
      REQUIRE_UNARY(j["ids"].is_array());
      CHECK_EQ(j["ids"].size(), 1);
      CHECK_EQ(j["ids"][0], info.id);
      REQUIRE_UNARY(j["ranks"].is_array());
      CHECK_EQ(j["ranks"].size(), 1);
      CHECK_EQ(j["ranks"][0], info.rank);
    }
  }
}

TEST_CASE("NormalNbrInfo from_json")
{
  for (int id : { 1, 2, 3 }) {
    for (int rank : { 0, 1, 2 }) {

      nlohmann::json j;
      j["type"] = "NORMAL";
      j["ids"] = { id };
      j["ranks"] = { rank };

      NormalNbrInfo<2> info = j.get<NormalNbrInfo<2>>();
      CHECK_EQ(info.id, id);
      CHECK_EQ(info.rank, rank);
    }
  }
}

TEST_CASE("CoarseNbrInfo to_json")
{
  for (int id : { 1, 2, 3 }) {
    for (int rank : { 0, 1, 2 }) {
      for (Orthant<2> orth_on_coarse : { Orthant<2>::sw(), Orthant<3>::se(), Orthant<2>::nw() }) {
        CoarseNbrInfo<2> info;
        info.id = id;
        info.rank = rank;
        info.orth_on_coarse = orth_on_coarse;

        ThunderEgg::tpl::nlohmann::json j = info;

        CHECK_EQ(j["type"], "COARSE");
        REQUIRE_UNARY(j["ids"].is_array());
        CHECK_EQ(j["ids"].size(), 1);
        CHECK_EQ(j["ids"][0], info.id);
        REQUIRE_UNARY(j["ranks"].is_array());
        CHECK_EQ(j["ranks"].size(), 1);
        CHECK_EQ(j["ranks"][0], info.rank);
        CHECK_EQ(j["orth_on_coarse"].get<Orthant<2>>(), info.orth_on_coarse);
      }
    }
  }
}

TEST_CASE("CoarseNbrInfo from_json")
{
  for (int id : { 1, 2, 3 }) {
    for (int rank : { 0, 1, 2 }) {
      for (Orthant<2> orth_on_coarse : { Orthant<2>::sw(), Orthant<3>::se(), Orthant<2>::nw() }) {

        ThunderEgg::tpl::nlohmann::json j;
        j["type"] = "COARSE";
        j["ids"] = { id };
        j["ranks"] = { rank };
        j["orth_on_coarse"] = orth_on_coarse;

        CoarseNbrInfo<2> info = j.get<CoarseNbrInfo<2>>();
        CHECK_EQ(info.id, id);
        CHECK_EQ(info.rank, rank);
        CHECK_EQ(info.orth_on_coarse, orth_on_coarse);
      }
    }
  }
}

TEST_CASE("NormalNbrInfo to_json")
{
  for (int id : { 1, 2, 3 }) {
    for (int rank : { 0, 1, 2 }) {
      NormalNbrInfo<2> info;
      info.id = id;
      info.rank = rank;

      nlohmann::json j = info;

      CHECK_EQ(j["type"], "NORMAL");
      REQUIRE_UNARY(j["ids"].is_array());
      CHECK_EQ(j["ids"].size(), 1);
      CHECK_EQ(j["ids"][0], info.id);
      REQUIRE_UNARY(j["ranks"].is_array());
      CHECK_EQ(j["ranks"].size(), 1);
      CHECK_EQ(j["ranks"][0], info.rank);
    }
  }
}

TEST_CASE("NormalNbrInfo from_json")
{
  for (int id : { 1, 2, 3 }) {
    for (int rank : { 0, 1, 2 }) {

      nlohmann::json j;
      j["type"] = "NORMAL";
      j["ids"] = { id };
      j["ranks"] = { rank };

      NormalNbrInfo<2> info = j.get<NormalNbrInfo<2>>();
      CHECK_EQ(info.id, id);
      CHECK_EQ(info.rank, rank);
    }
  }
}

TEST_CASE("PatchInfo<2> to_json no children")
{
  PatchInfo<2> d;
  d.id = 9;
  d.rank = 0;
  d.parent_id = 2;
  d.parent_rank = 3;
  d.orth_on_parent = Orthant<2>::nw();
  d.starts = { 1, 2 };
  d.spacings = { 0.1, 0.2 };
  d.ns = { 10, 20 };
  d.setNbrInfo(Side<2>::north(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Side<2>::east(), new CoarseNbrInfo<1>(2, Orthant<1>::lower()));
  d.setNbrInfo(Side<2>::south(), new FineNbrInfo<1>({ 3, 4 }));
  d.setNbrInfo(Corner<2>::sw(), new NormalNbrInfo<0>(1));
  d.setNbrInfo(Corner<2>::se(), new CoarseNbrInfo<0>(2, Orthant<0>(0)));
  d.setNbrInfo(Corner<2>::nw(), new FineNbrInfo<0>({ 1 }));

  nlohmann::json j = d;

  CHECK_EQ(j["id"], d.id);
  CHECK_EQ(j["parent_id"], d.parent_id);
  CHECK_EQ(j["parent_rank"], d.parent_rank);
  CHECK_EQ(j["orth_on_parent"], "TNW");
  CHECK_EQ(j["rank"], d.rank);
  CHECK_EQ(j["child_ids"], nullptr);
  CHECK_EQ(j["child_ranks"], nullptr);

  REQUIRE(j["starts"].is_array());
  REQUIRE(j["starts"].size() == 2);
  CHECK(j["starts"][0] == d.starts[0]);
  CHECK(j["starts"][1] == d.starts[1]);

  REQUIRE(j["lengths"].is_array());
  REQUIRE(j["lengths"].size() == 2);
  CHECK(j["lengths"][0] == d.spacings[0] * d.ns[0]);
  CHECK(j["lengths"][1] == d.spacings[1] * d.ns[1]);

  REQUIRE_UNARY(j["nbrs"].is_array());
  REQUIRE_EQ(j["nbrs"].size(), 3);

  CHECK_EQ(j["nbrs"][0]["type"], "COARSE");
  CHECK_EQ(j["nbrs"][0]["side"], "EAST");

  CHECK_EQ(j["nbrs"][1]["type"], "FINE");
  CHECK_EQ(j["nbrs"][1]["side"], "SOUTH");

  CHECK_EQ(j["nbrs"][2]["type"], "NORMAL");
  CHECK_EQ(j["nbrs"][2]["side"], "NORTH");

  REQUIRE_UNARY(j["corner_nbrs"].is_array());
  REQUIRE_EQ(j["corner_nbrs"].size(), 3);

  CHECK(j["corner_nbrs"][0]["type"] == "NORMAL");
  CHECK(j["corner_nbrs"][0]["corner"] == "SW");

  CHECK(j["corner_nbrs"][1]["type"] == "COARSE");
  CHECK(j["corner_nbrs"][1]["corner"] == "SE");

  CHECK(j["corner_nbrs"][2]["type"] == "FINE");
  CHECK(j["corner_nbrs"][2]["corner"] == "NW");
}

TEST_CASE("PatchInfo<2> to_json no children no neighbors")
{
  PatchInfo<2> d;
  d.id = 9;
  d.rank = 0;
  d.parent_id = 2;
  d.parent_rank = 3;
  d.refine_level = 329;
  d.starts = { 1, 2 };
  d.spacings = { 0.1, 0.2 };
  d.ns = { 10, 20 };

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
  REQUIRE(j["starts"].size() == 2);
  CHECK(j["starts"][0] == d.starts[0]);
  CHECK(j["starts"][1] == d.starts[1]);

  REQUIRE(j["lengths"].is_array());
  REQUIRE(j["lengths"].size() == 2);
  CHECK(j["lengths"][0] == d.spacings[0] * d.ns[0]);
  CHECK(j["lengths"][1] == d.spacings[1] * d.ns[1]);

  REQUIRE(j["nbrs"].is_array());
  REQUIRE(j["nbrs"].size() == 0);
}

TEST_CASE("PatchInfo<2> to_json with children")
{
  PatchInfo<2> d;
  d.id = 9;
  d.rank = 0;
  d.parent_id = 2;
  d.parent_rank = 3;
  d.refine_level = 329;
  d.starts = { 1, 2 };
  d.spacings = { 0.1, 0.2 };
  d.ns = { 10, 20 };
  d.child_ids = { 3, 4, 5, 6 };
  d.child_ranks = { 1, 2, 3, 4 };
  d.setNbrInfo(Side<2>::north(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Side<2>::east(), new CoarseNbrInfo<1>(2, Orthant<1>::lower()));
  d.setNbrInfo(Side<2>::south(), new FineNbrInfo<1>({ 3, 4 }));

  nlohmann::json j = d;

  CHECK(j["id"] == d.id);
  CHECK(j["parent_id"] == d.parent_id);
  CHECK(j["parent_rank"] == d.parent_rank);
  CHECK(j["rank"] == d.rank);
  CHECK(j["refine_level"] == 329);

  REQUIRE(j["child_ids"].is_array());
  REQUIRE(j["child_ids"].size() == 4);
  CHECK(j["child_ids"][0] == d.child_ids[0]);
  CHECK(j["child_ids"][1] == d.child_ids[1]);
  CHECK(j["child_ids"][2] == d.child_ids[2]);
  CHECK(j["child_ids"][3] == d.child_ids[3]);

  REQUIRE(j["child_ranks"].is_array());
  REQUIRE(j["child_ranks"].size() == 4);
  CHECK(j["child_ranks"][0] == d.child_ranks[0]);
  CHECK(j["child_ranks"][1] == d.child_ranks[1]);
  CHECK(j["child_ranks"][2] == d.child_ranks[2]);
  CHECK(j["child_ranks"][3] == d.child_ranks[3]);

  REQUIRE(j["starts"].is_array());
  REQUIRE(j["starts"].size() == 2);
  CHECK(j["starts"][0] == d.starts[0]);
  CHECK(j["starts"][1] == d.starts[1]);

  REQUIRE(j["lengths"].is_array());
  REQUIRE(j["lengths"].size() == 2);
  CHECK(j["lengths"][0] == d.spacings[0] * d.ns[0]);
  CHECK(j["lengths"][1] == d.spacings[1] * d.ns[1]);

  REQUIRE(j["nbrs"].is_array());
  REQUIRE(j["nbrs"].size() == 3);

  CHECK(j["nbrs"][0]["type"] == "COARSE");
  CHECK(j["nbrs"][0]["side"] == "EAST");

  CHECK(j["nbrs"][1]["type"] == "FINE");
  CHECK(j["nbrs"][1]["side"] == "SOUTH");

  CHECK(j["nbrs"][2]["type"] == "NORMAL");
  CHECK(j["nbrs"][2]["side"] == "NORTH");
}

TEST_CASE("PatchInfo<2> from_json no children")
{
  nlohmann::json j;
  j["id"] = 9;
  j["rank"] = 3;
  j["refine_level"] = 329;
  j["parent_id"] = 2;
  j["parent_rank"] = 3;
  j["starts"] = { 1, 2 };
  j["lengths"] = { 10, 20 };
  j["nbrs"] = { NormalNbrInfo<1>(1),
                CoarseNbrInfo<1>(2, Orthant<1>::lower()),
                FineNbrInfo<1>({ 3, 4 }) };
  j["nbrs"][0]["side"] = "NORTH";
  j["nbrs"][1]["side"] = "EAST";
  j["nbrs"][2]["side"] = "SOUTH";
  j["corner_nbrs"] = { NormalNbrInfo<0>(1),
                       CoarseNbrInfo<0>(2, Orthant<0>(0)),
                       FineNbrInfo<0>({ 1 }) };
  j["corner_nbrs"][0]["corner"] = "SW";
  j["corner_nbrs"][1]["corner"] = "SE";
  j["corner_nbrs"][2]["corner"] = "NW";

  PatchInfo<2> d = j.get<PatchInfo<2>>();
  CHECK(d.id == 9);
  CHECK(d.rank == 3);
  CHECK(d.refine_level == 329);
  CHECK(d.parent_id == 2);
  CHECK(d.parent_rank == 3);
  CHECK(d.orth_on_parent == Orthant<2>::null());
  CHECK(d.starts[0] == 1);
  CHECK(d.starts[1] == 2);
  CHECK(d.starts[2] == 3);
  CHECK(d.spacings[0] == 10);
  CHECK(d.spacings[1] == 20);
  CHECK(d.spacings[2] == 30);
  CHECK(d.ns[0] == 1);
  CHECK(d.ns[1] == 1);
  CHECK(d.ns[2] == 1);
  CHECK_FALSE(d.hasNbr(Side<2>::west()));
  CHECK(d.hasNbr(Side<2>::east()));
  CHECK(d.getNbrType(Side<2>::east()) == NbrType::Coarse);
  CHECK(d.hasNbr(Side<2>::south()));
  CHECK(d.getNbrType(Side<2>::south()) == NbrType::Fine);
  CHECK(d.hasNbr(Side<2>::north()));
  CHECK(d.getNbrType(Side<2>::north()) == NbrType::Normal);

  CHECK(d.hasNbr(Corner<2>::sw()));
  CHECK(d.getNbrType(Corner<2>::sw()) == NbrType::Normal);
  CHECK(d.hasNbr(Corner<2>::se()));
  CHECK(d.getNbrType(Corner<2>::se()) == NbrType::Coarse);
  CHECK(d.hasNbr(Corner<2>::nw()));
  CHECK(d.getNbrType(Corner<2>::nw()) == NbrType::Fine);
  CHECK_FALSE(d.hasNbr(Corner<2>::ne()));
}

TEST_CASE("PatchInfo<2> from_json with children")
{
  nlohmann::json j;
  j["id"] = 9;
  j["rank"] = 3;
  j["refine_level"] = 329;
  j["parent_id"] = 2;
  j["parent_rank"] = 3;
  j["orth_on_parent"] = "NW";
  j["starts"] = { 1, 2 };
  j["lengths"] = { 10, 20 };
  j["child_ids"] = { 1, 2, 3, 4 };
  j["child_ranks"] = { 0, 1, 2, 3 };
  j["nbrs"] = { NormalNbrInfo<1>(1),
                CoarseNbrInfo<1>(2, Orthant<1>::lower()),
                FineNbrInfo<1>({ 3, 4 }) };
  j["nbrs"][0]["side"] = "NORTH";
  j["nbrs"][1]["side"] = "EAST";
  j["nbrs"][2]["side"] = "SOUTH";

  PatchInfo<2> d = j.get<PatchInfo<2>>();
  CHECK(d.id == 9);
  CHECK(d.rank == 3);
  CHECK(d.refine_level == 329);
  CHECK(d.parent_id == 2);
  CHECK(d.parent_rank == 3);
  CHECK(d.orth_on_parent == Orthant<2>::nw());
  CHECK(d.starts[0] == 1);
  CHECK(d.starts[1] == 2);
  CHECK(d.spacings[0] == 10);
  CHECK(d.spacings[1] == 20);
  CHECK(d.ns[0] == 1);
  CHECK(d.ns[1] == 1);
  CHECK(d.child_ids[0] == 1);
  CHECK(d.child_ids[1] == 2);
  CHECK(d.child_ids[2] == 3);
  CHECK(d.child_ids[3] == 4);
  CHECK(d.child_ranks[0] == 0);
  CHECK(d.child_ranks[1] == 1);
  CHECK(d.child_ranks[2] == 2);
  CHECK(d.child_ranks[3] == 3);
  CHECK_FALSE(d.hasNbr(Side<2>::west()));
  CHECK(d.hasNbr(Side<2>::east()));
  CHECK(d.getNbrType(Side<2>::east()) == NbrType::Coarse);
  CHECK(d.hasNbr(Side<2>::south()));
  CHECK(d.getNbrType(Side<2>::south()) == NbrType::Fine);
  CHECK(d.hasNbr(Side<2>::north()));
  CHECK(d.getNbrType(Side<2>::north()) == NbrType::Normal);
}

TEST_CASE("PatchInfo<2> copy constructor")
{
  PatchInfo<2> d;
  d.id = 9;
  d.local_index = 10;
  d.global_index = 10;
  d.rank = 0;
  d.parent_id = 2;
  d.parent_rank = 3;
  d.num_ghost_cells = 239;
  d.refine_level = 329;
  d.starts = { 1, 2 };
  d.spacings = { 0.1, 0.2 };
  d.ns = { 10, 20 };
  d.child_ids = { 3, 4, 5, 6 };
  d.child_ranks = { 1, 2, 3, 4 };
  d.setNbrInfo(Side<2>::north(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Side<2>::east(), new CoarseNbrInfo<1>(2, Orthant<1>::lower()));
  d.setNbrInfo(Side<2>::south(), new FineNbrInfo<1>({ 3, 4 }));
  d.setNbrInfo(Corner<2>::sw(), new NormalNbrInfo<0>(1));
  d.setNbrInfo(Corner<2>::se(), new CoarseNbrInfo<0>(2, Orthant<0>(0)));
  d.setNbrInfo(Corner<2>::nw(), new FineNbrInfo<0>({ 1 }));

  PatchInfo<2> d2(d);

  CHECK_EQ(d.id, d2.id);
  CHECK_EQ(d.local_index, d2.global_index);
  CHECK_EQ(d.rank, d2.rank);
  CHECK_EQ(d.parent_id, d2.parent_id);
  CHECK_EQ(d.parent_rank, d2.parent_rank);
  CHECK_EQ(d.num_ghost_cells, d2.num_ghost_cells);
  CHECK_EQ(d.refine_level, d2.refine_level);
  CHECK_EQ(d.starts, d2.starts);
  CHECK_EQ(d.spacings, d2.spacings);
  CHECK_EQ(d.ns, d2.ns);
  CHECK_EQ(d.child_ids, d2.child_ids);
  CHECK_EQ(d.child_ranks, d2.child_ranks);

  for (Side<2> s : Side<2>::getValues()) {
    REQUIRE_EQ(d.hasNbr(s), d2.hasNbr(s));
    if (d.hasNbr(s)) {
      switch (d.getNbrType(s)) {
        case NbrType::Normal:
          CHECK_EQ(d.getNormalNbrInfo(s).id, d2.getNormalNbrInfo(s).id);
          CHECK_NE(&d.getNormalNbrInfo(s), &d2.getNormalNbrInfo(s));
          break;
        case NbrType::Fine:
          CHECK_EQ(d.getFineNbrInfo(s).ids[0], d2.getFineNbrInfo(s).ids[0]);
          CHECK_NE(&d.getFineNbrInfo(s), &d2.getFineNbrInfo(s));
          break;
        case NbrType::Coarse:
          CHECK_EQ(d.getCoarseNbrInfo(s).id, d2.getCoarseNbrInfo(s).id);
          CHECK_NE(&d.getCoarseNbrInfo(s), &d2.getCoarseNbrInfo(s));
          break;
      }
    }
  }
  for (Corner<2> c : Corner<2>::getValues()) {
    REQUIRE_EQ(d.hasNbr(c), d2.hasNbr(c));
    if (d.hasNbr(c)) {
      switch (d.getNbrType(c)) {
        case NbrType::Normal:
          CHECK_EQ(d.getNormalNbrInfo(c).id, d2.getNormalNbrInfo(c).id);
          CHECK_NE(&d.getNormalNbrInfo(c), &d2.getNormalNbrInfo(c));
          break;
        case NbrType::Fine:
          CHECK_EQ(d.getFineNbrInfo(c).ids[0], d2.getFineNbrInfo(c).ids[0]);
          CHECK_NE(&d.getFineNbrInfo(c), &d2.getFineNbrInfo(c));
          break;
        case NbrType::Coarse:
          CHECK_EQ(d.getCoarseNbrInfo(c).id, d2.getCoarseNbrInfo(c).id);
          CHECK_NE(&d.getCoarseNbrInfo(c), &d2.getCoarseNbrInfo(c));
          break;
      }
    }
  }
}