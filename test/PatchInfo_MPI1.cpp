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

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::tpl;

TEST_CASE("PatchInfo getNbrIds NormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<2>* nbr_info = new NormalNbrInfo<2>(2);
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE(ids.size() == 1);
  CHECK(ids[0] == 2);
}
TEST_CASE("PatchInfo getNbrRanks NormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<2>* nbr_info = new NormalNbrInfo<2>(2);
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE(ranks.size() == 1);
  CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists NormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<2>* nbr_info = new NormalNbrInfo<2>(2);
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_index == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist NormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<2>* nbr_info = new NormalNbrInfo<2>(2);
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_index == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists NormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<2>* nbr_info = new NormalNbrInfo<2>(2);
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK(nbr_info->global_index == 30);
}
TEST_CASE("PatchInfo getNbrIds CoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<2>* nbr_info = new CoarseNbrInfo<2>(2, Orthant<2>::sw());
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE(ids.size() == 1);
  CHECK(ids[0] == 2);
}
TEST_CASE("PatchInfo getNbrRanks CoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<2>* nbr_info = new CoarseNbrInfo<2>(2, Orthant<2>::sw());
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE(ranks.size() == 1);
  CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists CoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<2>* nbr_info = new CoarseNbrInfo<2>(2, Orthant<2>::sw());
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_index == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist CoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<2>* nbr_info = new CoarseNbrInfo<2>(2, Orthant<2>::sw());
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_index == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists CoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<2>* nbr_info = new CoarseNbrInfo<2>(2, Orthant<2>::sw());
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK(nbr_info->global_index == 30);
}
TEST_CASE("PatchInfo getNbrIds FineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<2>* nbr_info = new FineNbrInfo<2>({ 1, 2, 3, 4 });
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE(ids.size() == 4);
  CHECK(ids[0] == 1);
  CHECK(ids[1] == 2);
  CHECK(ids[2] == 3);
  CHECK(ids[3] == 4);
}
TEST_CASE("PatchInfo getNbrRanks FineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<2>* nbr_info = new FineNbrInfo<2>({ 1, 2, 3, 4 });
  nbr_info->ranks = { 3, 4, 6, 7 };
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE(ranks.size() == 4);
  CHECK(ranks[0] == 3);
  CHECK(ranks[1] == 4);
  CHECK(ranks[2] == 6);
  CHECK(ranks[3] == 7);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists FineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<2>* nbr_info = new FineNbrInfo<2>({ 2, 3, 4, 5 });
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  id_to_local_index_map[3] = 31;
  id_to_local_index_map[4] = 32;
  id_to_local_index_map[5] = 33;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_indexes[0] == 30);
  CHECK(nbr_info->local_indexes[1] == 31);
  CHECK(nbr_info->local_indexes[2] == 32);
  CHECK(nbr_info->local_indexes[3] == 33);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist FineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<2>* nbr_info = new FineNbrInfo<2>({ 2, 3, 4, 5 });
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  id_to_local_index_map[3] = 31;
  id_to_local_index_map[4] = 32;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_indexes[0] == 30);
  CHECK(nbr_info->local_indexes[1] == 31);
  CHECK(nbr_info->local_indexes[2] == 32);
  CHECK(nbr_info->local_indexes[3] == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists FineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<2>* nbr_info = new FineNbrInfo<2>({ 2, 3, 4, 5 });
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  id_to_global_index_map[3] = 31;
  id_to_global_index_map[4] = 32;
  id_to_global_index_map[5] = 33;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK(nbr_info->global_indexes[0] == 30);
  CHECK(nbr_info->global_indexes[1] == 31);
  CHECK(nbr_info->global_indexes[2] == 32);
  CHECK(nbr_info->global_indexes[3] == 33);
}
TEST_CASE("PatchInfo getNbrIds EdgeNormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE(ids.size() == 1);
  CHECK(ids[0] == 2);
}
TEST_CASE("PatchInfo getNbrRanks EdgeNormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE(ranks.size() == 1);
  CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists EdgeNormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_index == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist EdgeNormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_index == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists EdgeNormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK(nbr_info->global_index == 30);
}
TEST_CASE("PatchInfo getNbrIds EdgeCoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE(ids.size() == 1);
  CHECK(ids[0] == 2);
}
TEST_CASE("PatchInfo getNbrRanks EdgeCoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE(ranks.size() == 1);
  CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists EdgeCoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_index == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist EdgeCoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_index == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists EdgeCoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK(nbr_info->global_index == 30);
}
TEST_CASE("PatchInfo getNbrIds EdgeFineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 1, 2 });
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE(ids.size() == 2);
  CHECK(ids[0] == 1);
  CHECK(ids[1] == 2);
}
TEST_CASE("PatchInfo getNbrRanks EdgeFineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 1, 2 });
  nbr_info->ranks = { 3, 4 };
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE(ranks.size() == 2);
  CHECK(ranks[0] == 3);
  CHECK(ranks[1] == 4);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists EdgeFineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 2, 3 });
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  id_to_local_index_map[3] = 31;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_indexes[0] == 30);
  CHECK(nbr_info->local_indexes[1] == 31);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist EdgeFineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 2, 3 });
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_indexes[0] == 30);
  CHECK(nbr_info->local_indexes[1] == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists EdgeFineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 2, 3 });
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  id_to_global_index_map[3] = 31;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK(nbr_info->global_indexes[0] == 30);
  CHECK(nbr_info->global_indexes[1] == 31);
}
TEST_CASE("PatchInfo getNbrIds CornerNormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE(ids.size() == 1);
  CHECK(ids[0] == 2);
}
TEST_CASE("PatchInfo getNbrRanks CornerNormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE(ranks.size() == 1);
  CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists CornerNormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_index == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist CornerNormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_index == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists CornerNormalNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK(nbr_info->global_index == 30);
}
TEST_CASE("PatchInfo getNbrIds CornerCoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE(ids.size() == 1);
  CHECK(ids[0] == 2);
}
TEST_CASE("PatchInfo getNbrRanks CornerCoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE(ranks.size() == 1);
  CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists CornerCoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_index == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist CornerCoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_index == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists CornerCoarseNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK(nbr_info->global_index == 30);
}
TEST_CASE("PatchInfo getNbrIds CornerFineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 1 });
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE(ids.size() == 1);
  CHECK(ids[0] == 1);
}
TEST_CASE("PatchInfo getNbrRanks CornerFineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 1 });
  nbr_info->ranks = { 3 };
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE(ranks.size() == 1);
  CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists CornerFineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 2 });
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_indexes[0] == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist CornerFineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 2 });
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK(nbr_info->local_indexes[0] == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists CornerFineNbrInfo", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 2 });
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK(nbr_info->global_indexes[0] == 30);
}
TEST_CASE("PatchInfo Serialization/Deserialization", "[PatchInfo]")
{
  PatchInfo<3>* d_ptr = new PatchInfo<3>;
  PatchInfo<3>& d = *d_ptr;
  d.id = 0;
  d.setNbrInfo(Side<3>::north(), new NormalNbrInfo<2>(1));
  d.setNbrInfo(Side<3>::east(), new CoarseNbrInfo<2>(2, Orthant<2>::nw()));
  d.setNbrInfo(Side<3>::south(), new FineNbrInfo<2>({ 3, 4, 5, 6 }));
  d.setNbrInfo(Corner<3>::bsw(), new NormalNbrInfo<0>(1));
  d.setNbrInfo(Corner<3>::tse(), new CoarseNbrInfo<0>(2, Orthant<0>::null()));
  d.setNbrInfo(Corner<3>::bnw(), new FineNbrInfo<0>({ 1 }));
  d.setNbrInfo(Edge::sw(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Edge::bn(), new CoarseNbrInfo<1>(2, Orthant<1>::lower()));
  d.setNbrInfo(Edge::tw(), new FineNbrInfo<1>({ 1, 2 }));

  // serialize and then deserialize
  char* buff = new char[d.serialize(nullptr)];
  d.serialize(buff);
  delete d_ptr;
  PatchInfo<3> out;
  out.deserialize(buff);
  delete[] buff;

  // check that deserialized version has the same information
  REQUIRE(out.id == 0);

  REQUIRE(!out.hasNbr(Side<3>::west()));

  REQUIRE(out.hasNbr(Side<3>::east()));
  REQUIRE(out.getNbrType(Side<3>::east()) == NbrType::Coarse);
  REQUIRE(out.getCoarseNbrInfo(Side<3>::east()).id == 2);
  REQUIRE(out.getCoarseNbrInfo(Side<3>::east()).orth_on_coarse == Orthant<2>::nw());

  REQUIRE(out.hasNbr(Side<3>::south()));
  REQUIRE(out.getNbrType(Side<3>::south()) == NbrType::Fine);
  REQUIRE(out.getFineNbrInfo(Side<3>::south()).ids[0] == 3);
  REQUIRE(out.getFineNbrInfo(Side<3>::south()).ids[1] == 4);
  REQUIRE(out.getFineNbrInfo(Side<3>::south()).ids[2] == 5);
  REQUIRE(out.getFineNbrInfo(Side<3>::south()).ids[3] == 6);

  REQUIRE(out.hasNbr(Side<3>::north()));
  REQUIRE(out.getNbrType(Side<3>::north()) == NbrType::Normal);
  REQUIRE(out.getNormalNbrInfo(Side<3>::north()).id == 1);

  REQUIRE(!out.hasNbr(Side<3>::bottom()));
  REQUIRE(!out.hasNbr(Side<3>::top()));

  // Corners

  REQUIRE(out.hasNbr(Corner<3>::bsw()));
  REQUIRE(out.getNbrType(Corner<3>::bsw()) == NbrType::Normal);
  REQUIRE(out.getNormalNbrInfo(Corner<3>::bsw()).id == 1);

  REQUIRE(!out.hasNbr(Corner<3>::bse()));

  REQUIRE(out.hasNbr(Corner<3>::bnw()));
  REQUIRE(out.getNbrType(Corner<3>::bnw()) == NbrType::Fine);
  REQUIRE(out.getFineNbrInfo(Corner<3>::bnw()).ids[0] == 1);

  REQUIRE(!out.hasNbr(Corner<3>::bne()));
  REQUIRE(!out.hasNbr(Corner<3>::tsw()));

  REQUIRE(out.hasNbr(Corner<3>::tse()));
  REQUIRE(out.getNbrType(Corner<3>::tse()) == NbrType::Coarse);
  REQUIRE(out.getCoarseNbrInfo(Corner<3>::tse()).id == 2);
  REQUIRE(out.getCoarseNbrInfo(Corner<3>::tse()).orth_on_coarse == Orthant<0>::null());

  REQUIRE(!out.hasNbr(Corner<3>::tnw()));
  REQUIRE(!out.hasNbr(Corner<3>::tne()));

  // Edges

  REQUIRE(!out.hasNbr(Edge::bs()));
  REQUIRE(!out.hasNbr(Edge::tn()));

  REQUIRE(out.hasNbr(Edge::bn()));
  REQUIRE(out.getNbrType(Edge::bn()) == NbrType::Coarse);
  REQUIRE(out.getCoarseNbrInfo(Edge::bn()).id == 2);
  REQUIRE(out.getCoarseNbrInfo(Edge::bn()).orth_on_coarse == Orthant<1>::lower());

  REQUIRE(!out.hasNbr(Edge::ts()));
  REQUIRE(!out.hasNbr(Edge::bw()));
  REQUIRE(!out.hasNbr(Edge::te()));
  REQUIRE(!out.hasNbr(Edge::be()));

  REQUIRE(out.hasNbr(Edge::tw()));
  REQUIRE(out.getNbrType(Edge::tw()) == NbrType::Fine);
  REQUIRE(out.getFineNbrInfo(Edge::tw()).ids[0] == 1);
  REQUIRE(out.getFineNbrInfo(Edge::tw()).ids[1] == 2);

  REQUIRE(out.hasNbr(Edge::sw()));
  REQUIRE(out.getNbrType(Edge::sw()) == NbrType::Normal);
  REQUIRE(out.getNormalNbrInfo(Edge::sw()).id == 1);

  REQUIRE(!out.hasNbr(Edge::ne()));
  REQUIRE(!out.hasNbr(Edge::se()));
  REQUIRE(!out.hasNbr(Edge::nw()));
}
TEST_CASE("PatchInfo Default Values", "[PatchInfo]")
{
  PatchInfo<3> pinfo;
  CHECK(pinfo.id == 0);
  CHECK(pinfo.local_index == 0);
  CHECK(pinfo.global_index == 0);
  CHECK(pinfo.refine_level == -1);
  CHECK(pinfo.parent_id == -1);
  CHECK(pinfo.parent_rank == -1);
  for (int child_id : pinfo.child_ids) {
    CHECK(child_id == -1);
  }
  for (int child_rank : pinfo.child_ids) {
    CHECK(child_rank == -1);
  }
  CHECK(pinfo.num_ghost_cells == 0);
  CHECK(pinfo.rank == -1);
  CHECK(pinfo.orth_on_parent == Orthant<3>::null());
  for (int n : pinfo.ns) {
    CHECK(n == 1);
  }
  for (double start : pinfo.starts) {
    CHECK(start == 0);
  }
  for (double spacing : pinfo.spacings) {
    CHECK(spacing == 1);
  }
  CHECK_FALSE(pinfo.hasNbr());
}

TEST_CASE("PatchInfo copy constructor", "[PatchInfo]")
{
  PatchInfo<3> d;
  d.id = 9;
  d.local_index = 10;
  d.global_index = 10;
  d.rank = 0;
  d.parent_id = 2;
  d.parent_rank = 3;
  d.num_ghost_cells = 239;
  d.refine_level = 329;
  d.starts = { 1, 2, 3 };
  d.spacings = { 0.1, 0.2, 0.3 };
  d.ns = { 10, 20, 30 };
  d.child_ids = { 3, 4, 5, 6, 7, 8, 9, 10 };
  d.child_ranks = { 1, 2, 3, 4, 5, 6, 7, 8 };
  d.setNbrInfo(Side<3>::north(), new NormalNbrInfo<2>(1));
  d.setNbrInfo(Side<3>::east(), new CoarseNbrInfo<2>(2, Orthant<2>::nw()));
  d.setNbrInfo(Side<3>::south(), new FineNbrInfo<2>({ 3, 4, 5, 6 }));
  d.setNbrInfo(Corner<3>::bsw(), new NormalNbrInfo<0>(1));
  d.setNbrInfo(Corner<3>::tse(), new CoarseNbrInfo<0>(2, Orthant<0>(0)));
  d.setNbrInfo(Corner<3>::bnw(), new FineNbrInfo<0>({ 1 }));
  d.setNbrInfo(Edge::sw(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Edge::bn(), new CoarseNbrInfo<1>(2, Orthant<1>::lower()));
  d.setNbrInfo(Edge::tw(), new FineNbrInfo<1>({ 1, 2 }));

  PatchInfo<3> d2(d);

  CHECK(d.id == d2.id);
  CHECK(d.local_index == d2.global_index);
  CHECK(d.rank == d2.rank);
  CHECK(d.parent_id == d2.parent_id);
  CHECK(d.parent_rank == d2.parent_rank);
  CHECK(d.num_ghost_cells == d2.num_ghost_cells);
  CHECK(d.refine_level == d2.refine_level);
  CHECK(d.starts == d2.starts);
  CHECK(d.spacings == d2.spacings);
  CHECK(d.ns == d2.ns);
  CHECK(d.child_ids == d2.child_ids);
  CHECK(d.child_ranks == d2.child_ranks);

  for (Side<3> s : Side<3>::getValues()) {
    REQUIRE(d.hasNbr(s) == d2.hasNbr(s));
    if (d.hasNbr(s)) {
      switch (d.getNbrType(s)) {
        case NbrType::Normal:
          CHECK(d.getNormalNbrInfo(s).id == d2.getNormalNbrInfo(s).id);
          CHECK(&d.getNormalNbrInfo(s) != &d2.getNormalNbrInfo(s));
          break;
        case NbrType::Fine:
          CHECK(d.getFineNbrInfo(s).ids[0] == d2.getFineNbrInfo(s).ids[0]);
          CHECK(&d.getFineNbrInfo(s) != &d2.getFineNbrInfo(s));
          break;
        case NbrType::Coarse:
          CHECK(d.getCoarseNbrInfo(s).id == d2.getCoarseNbrInfo(s).id);
          CHECK(&d.getCoarseNbrInfo(s) != &d2.getCoarseNbrInfo(s));
          break;
      }
    }
  }
  for (Corner<3> c : Corner<3>::getValues()) {
    REQUIRE(d.hasNbr(c) == d2.hasNbr(c));
    if (d.hasNbr(c)) {
      switch (d.getNbrType(c)) {
        case NbrType::Normal:
          CHECK(d.getNormalNbrInfo(c).id == d2.getNormalNbrInfo(c).id);
          CHECK(&d.getNormalNbrInfo(c) != &d2.getNormalNbrInfo(c));
          break;
        case NbrType::Fine:
          CHECK(d.getFineNbrInfo(c).ids[0] == d2.getFineNbrInfo(c).ids[0]);
          CHECK(&d.getFineNbrInfo(c) != &d2.getFineNbrInfo(c));
          break;
        case NbrType::Coarse:
          CHECK(d.getCoarseNbrInfo(c).id == d2.getCoarseNbrInfo(c).id);
          CHECK(&d.getCoarseNbrInfo(c) != &d2.getCoarseNbrInfo(c));
          break;
      }
    }
  }
  for (Edge c : Edge::getValues()) {
    REQUIRE(d.hasNbr(c) == d2.hasNbr(c));
    if (d.hasNbr(c)) {
      switch (d.getNbrType(c)) {
        case NbrType::Normal:
          CHECK(d.getNormalNbrInfo(c).id == d2.getNormalNbrInfo(c).id);
          CHECK(&d.getNormalNbrInfo(c) != &d2.getNormalNbrInfo(c));
          break;
        case NbrType::Fine:
          CHECK(d.getFineNbrInfo(c).ids[0] == d2.getFineNbrInfo(c).ids[0]);
          CHECK(&d.getFineNbrInfo(c) != &d2.getFineNbrInfo(c));
          break;
        case NbrType::Coarse:
          CHECK(d.getCoarseNbrInfo(c).id == d2.getCoarseNbrInfo(c).id);
          CHECK(&d.getCoarseNbrInfo(c) != &d2.getCoarseNbrInfo(c));
          break;
      }
    }
  }
}
TEST_CASE("PatchInfo copy assignment", "[PatchInfo]")
{
  PatchInfo<3> d;
  d.id = 9;
  d.local_index = 10;
  d.global_index = 10;
  d.rank = 0;
  d.parent_id = 2;
  d.parent_rank = 3;
  d.num_ghost_cells = 239;
  d.refine_level = 329;
  d.starts = { 1, 2, 3 };
  d.spacings = { 0.1, 0.2, 0.3 };
  d.ns = { 10, 20, 30 };
  d.child_ids = { 3, 4, 5, 6, 7, 8, 9, 10 };
  d.child_ranks = { 1, 2, 3, 4, 5, 6, 7, 8 };
  d.setNbrInfo(Side<3>::north(), new NormalNbrInfo<2>(1));
  d.setNbrInfo(Side<3>::east(), new CoarseNbrInfo<2>(2, Orthant<2>::nw()));
  d.setNbrInfo(Side<3>::south(), new FineNbrInfo<2>({ 3, 4, 5, 6 }));
  d.setNbrInfo(Corner<3>::bsw(), new NormalNbrInfo<0>(1));
  d.setNbrInfo(Corner<3>::tse(), new CoarseNbrInfo<0>(2, Orthant<0>(0)));
  d.setNbrInfo(Corner<3>::bnw(), new FineNbrInfo<0>({ 1 }));
  d.setNbrInfo(Edge::sw(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Edge::bn(), new CoarseNbrInfo<1>(2, Orthant<1>::lower()));
  d.setNbrInfo(Edge::tw(), new FineNbrInfo<1>({ 1, 2 }));

  PatchInfo<3> d2;
  d2 = d;

  CHECK(d.id == d2.id);
  CHECK(d.local_index == d2.global_index);
  CHECK(d.rank == d2.rank);
  CHECK(d.parent_id == d2.parent_id);
  CHECK(d.parent_rank == d2.parent_rank);
  CHECK(d.num_ghost_cells == d2.num_ghost_cells);
  CHECK(d.refine_level == d2.refine_level);
  CHECK(d.starts == d2.starts);
  CHECK(d.spacings == d2.spacings);
  CHECK(d.ns == d2.ns);
  CHECK(d.child_ids == d2.child_ids);
  CHECK(d.child_ranks == d2.child_ranks);

  for (Side<3> s : Side<3>::getValues()) {
    REQUIRE(d.hasNbr(s) == d2.hasNbr(s));
    if (d.hasNbr(s)) {
      switch (d.getNbrType(s)) {
        case NbrType::Normal:
          CHECK(d.getNormalNbrInfo(s).id == d2.getNormalNbrInfo(s).id);
          CHECK(&d.getNormalNbrInfo(s) != &d2.getNormalNbrInfo(s));
          break;
        case NbrType::Fine:
          CHECK(d.getFineNbrInfo(s).ids[0] == d2.getFineNbrInfo(s).ids[0]);
          CHECK(&d.getFineNbrInfo(s) != &d2.getFineNbrInfo(s));
          break;
        case NbrType::Coarse:
          CHECK(d.getCoarseNbrInfo(s).id == d2.getCoarseNbrInfo(s).id);
          CHECK(&d.getCoarseNbrInfo(s) != &d2.getCoarseNbrInfo(s));
          break;
      }
    }
  }
  for (Corner<3> c : Corner<3>::getValues()) {
    REQUIRE(d.hasNbr(c) == d2.hasNbr(c));
    if (d.hasNbr(c)) {
      switch (d.getNbrType(c)) {
        case NbrType::Normal:
          CHECK(d.getNormalNbrInfo(c).id == d2.getNormalNbrInfo(c).id);
          CHECK(&d.getNormalNbrInfo(c) != &d2.getNormalNbrInfo(c));
          break;
        case NbrType::Fine:
          CHECK(d.getFineNbrInfo(c).ids[0] == d2.getFineNbrInfo(c).ids[0]);
          CHECK(&d.getFineNbrInfo(c) != &d2.getFineNbrInfo(c));
          break;
        case NbrType::Coarse:
          CHECK(d.getCoarseNbrInfo(c).id == d2.getCoarseNbrInfo(c).id);
          CHECK(&d.getCoarseNbrInfo(c) != &d2.getCoarseNbrInfo(c));
          break;
      }
    }
  }
  for (Edge c : Edge::getValues()) {
    REQUIRE(d.hasNbr(c) == d2.hasNbr(c));
    if (d.hasNbr(c)) {
      switch (d.getNbrType(c)) {
        case NbrType::Normal:
          CHECK(d.getNormalNbrInfo(c).id == d2.getNormalNbrInfo(c).id);
          CHECK(&d.getNormalNbrInfo(c) != &d2.getNormalNbrInfo(c));
          break;
        case NbrType::Fine:
          CHECK(d.getFineNbrInfo(c).ids[0] == d2.getFineNbrInfo(c).ids[0]);
          CHECK(&d.getFineNbrInfo(c) != &d2.getFineNbrInfo(c));
          break;
        case NbrType::Coarse:
          CHECK(d.getCoarseNbrInfo(c).id == d2.getCoarseNbrInfo(c).id);
          CHECK(&d.getCoarseNbrInfo(c) != &d2.getCoarseNbrInfo(c));
          break;
      }
    }
  }
}
