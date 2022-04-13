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

TEST_CASE("PatchInfo<2> getNbrIds NormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Side<1>::west(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 2);
}
TEST_CASE("PatchInfo<2> getNbrRanks NormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes exists NormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, 30);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes does not exist NormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, -1);
}
TEST_CASE("PatchInfo<2> setNeighborGlobalIndexes exists NormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_index, 30);
}
TEST_CASE("PatchInfo<2> getNbrIds CoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::sw());
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 2);
}
TEST_CASE("PatchInfo<2> getNbrRanks CoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::sw());
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes exists CoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::sw());
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, 30);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes does not exist CoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::sw());
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, -1);
}
TEST_CASE("PatchInfo<2> setNeighborGlobalIndexes exists CoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::sw());
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_index, 30);
}
TEST_CASE("PatchInfo<2> getNbrIds FineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 1, 2, 3, 4 });
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 4);
  CHECK_EQ(ids[0], 1);
  CHECK_EQ(ids[1], 2);
  CHECK_EQ(ids[2], 3);
  CHECK_EQ(ids[3], 4);
}
TEST_CASE("PatchInfo<2> getNbrRanks FineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 1, 2, 3, 4 });
  nbr_info->ranks = { 3, 4, 6, 7 };
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 4);
  CHECK_EQ(ranks[0], 3);
  CHECK_EQ(ranks[1], 4);
  CHECK_EQ(ranks[2], 6);
  CHECK_EQ(ranks[3], 7);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes exists FineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 2, 3, 4, 5 });
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  id_to_local_index_map[3] = 31;
  id_to_local_index_map[4] = 32;
  id_to_local_index_map[5] = 33;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_indexes[0], 30);
  CHECK_EQ(nbr_info->local_indexes[1], 31);
  CHECK_EQ(nbr_info->local_indexes[2], 32);
  CHECK_EQ(nbr_info->local_indexes[3], 33);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes does not exist FineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 2, 3, 4, 5 });
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  id_to_local_index_map[3] = 31;
  id_to_local_index_map[4] = 32;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_indexes[0], 30);
  CHECK_EQ(nbr_info->local_indexes[1], 31);
  CHECK_EQ(nbr_info->local_indexes[2], 32);
  CHECK_EQ(nbr_info->local_indexes[3], -1);
}
TEST_CASE("PatchInfo<2> setNeighborGlobalIndexes exists FineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 2, 3, 4, 5 });
  pinfo.setNbrInfo(Side<2>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  id_to_global_index_map[3] = 31;
  id_to_global_index_map[4] = 32;
  id_to_global_index_map[5] = 33;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_indexes[0], 30);
  CHECK_EQ(nbr_info->global_indexes[1], 31);
  CHECK_EQ(nbr_info->global_indexes[2], 32);
  CHECK_EQ(nbr_info->global_indexes[3], 33);
}
TEST_CASE("PatchInfo<2> getNbrIds EdgeNormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 2);
}
TEST_CASE("PatchInfo<2> getNbrRanks EdgeNormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes exists EdgeNormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, 30);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes does not exist EdgeNormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, -1);
}
TEST_CASE("PatchInfo<2> setNeighborGlobalIndexes exists EdgeNormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_index, 30);
}
TEST_CASE("PatchInfo<2> getNbrIds EdgeCoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 2);
}
TEST_CASE("PatchInfo<2> getNbrRanks EdgeCoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes exists EdgeCoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, 30);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes does not exist EdgeCoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, -1);
}
TEST_CASE("PatchInfo<2> setNeighborGlobalIndexes exists EdgeCoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_index, 30);
}
TEST_CASE("PatchInfo<2> getNbrIds EdgeFineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 1, 2 });
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 2);
  CHECK_EQ(ids[0], 1);
  CHECK_EQ(ids[1], 2);
}
TEST_CASE("PatchInfo<2> getNbrRanks EdgeFineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 1, 2 });
  nbr_info->ranks = { 3, 4 };
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 2);
  CHECK_EQ(ranks[0], 3);
  CHECK_EQ(ranks[1], 4);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes exists EdgeFineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 2, 3 });
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  id_to_local_index_map[3] = 31;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_indexes[0], 30);
  CHECK_EQ(nbr_info->local_indexes[1], 31);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes does not exist EdgeFineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 2, 3 });
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_indexes[0], 30);
  CHECK_EQ(nbr_info->local_indexes[1], -1);
}
TEST_CASE("PatchInfo<2> setNeighborGlobalIndexes exists EdgeFineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 2, 3 });
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  id_to_global_index_map[3] = 31;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_indexes[0], 30);
  CHECK_EQ(nbr_info->global_indexes[1], 31);
}
TEST_CASE("PatchInfo<2> getNbrIds CornerNormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 2);
}
TEST_CASE("PatchInfo<2> getNbrRanks CornerNormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes exists CornerNormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, 30);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes does not exist CornerNormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, -1);
}
TEST_CASE("PatchInfo<2> setNeighborGlobalIndexes exists CornerNormalNbrInfo")
{
  PatchInfo<2> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_index, 30);
}
TEST_CASE("PatchInfo<2> getNbrIds CornerCoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 2);
}
TEST_CASE("PatchInfo<2> getNbrRanks CornerCoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes exists CornerCoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, 30);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes does not exist CornerCoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, -1);
}
TEST_CASE("PatchInfo<2> setNeighborGlobalIndexes exists CornerCoarseNbrInfo")
{
  PatchInfo<2> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_index, 30);
}
TEST_CASE("PatchInfo<2> getNbrIds CornerFineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 1 });
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 1);
}
TEST_CASE("PatchInfo<2> getNbrRanks CornerFineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 1 });
  nbr_info->ranks = { 3 };
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes exists CornerFineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 2 });
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_indexes[0], 30);
}
TEST_CASE("PatchInfo<2> setNeighborLocalIndexes does not exist CornerFineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 2 });
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_indexes[0], -1);
}
TEST_CASE("PatchInfo<2> setNeighborGlobalIndexes exists CornerFineNbrInfo")
{
  PatchInfo<2> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 2 });
  pinfo.setNbrInfo(Corner<2>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_indexes[0], 30);
}
TEST_CASE("PatchInfo<2> Serialization/Deserialization")
{
  PatchInfo<2>* d_ptr = new PatchInfo<2>;
  PatchInfo<2>& d = *d_ptr;
  d.id = 0;
  d.setNbrInfo(Side<2>::north(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Side<2>::east(), new CoarseNbrInfo<1>(2, Orthant<1>::nw()));
  d.setNbrInfo(Side<2>::south(), new FineNbrInfo<1>({ 3, 4, 5, 6 }));
  d.setNbrInfo(Corner<2>::bsw(), new NormalNbrInfo<0>(1));
  d.setNbrInfo(Corner<2>::tse(), new CoarseNbrInfo<0>(2, Orthant<0>::null()));
  d.setNbrInfo(Corner<2>::bnw(), new FineNbrInfo<0>({ 1 }));
  d.setNbrInfo(Edge::sw(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Edge::bn(), new CoarseNbrInfo<1>(2, Orthant<1>::lower()));
  d.setNbrInfo(Edge::tw(), new FineNbrInfo<1>({ 1, 2 }));

  // serialize and then deserialize
  char* buff = new char[d.serialize(nullptr)];
  d.serialize(buff);
  delete d_ptr;
  PatchInfo<2> out;
  out.deserialize(buff);
  delete[] buff;

  // check that deserialized version has the same information
  REQUIRE_EQ(out.id, 0);

  REQUIRE_UNARY(!out.hasNbr(Side<2>::west()));

  REQUIRE_UNARY(out.hasNbr(Side<2>::east()));
  REQUIRE_EQ(out.getNbrType(Side<2>::east()), NbrType::Coarse);
  REQUIRE_EQ(out.getCoarseNbrInfo(Side<2>::east()).id, 2);
  REQUIRE_EQ(out.getCoarseNbrInfo(Side<2>::east()).orth_on_coarse, Orthant<1>::nw());

  REQUIRE_UNARY(out.hasNbr(Side<2>::south()));
  REQUIRE_EQ(out.getNbrType(Side<2>::south()), NbrType::Fine);
  REQUIRE_EQ(out.getFineNbrInfo(Side<2>::south()).ids[0], 3);
  REQUIRE_EQ(out.getFineNbrInfo(Side<2>::south()).ids[1], 4);
  REQUIRE_EQ(out.getFineNbrInfo(Side<2>::south()).ids[2], 5);
  REQUIRE_EQ(out.getFineNbrInfo(Side<2>::south()).ids[3], 6);

  REQUIRE_UNARY(out.hasNbr(Side<2>::north()));
  REQUIRE_EQ(out.getNbrType(Side<2>::north()), NbrType::Normal);
  REQUIRE_EQ(out.getNormalNbrInfo(Side<2>::north()).id, 1);

  REQUIRE_UNARY(!out.hasNbr(Side<2>::bottom()));
  REQUIRE_UNARY(!out.hasNbr(Side<2>::top()));

  // Corners

  REQUIRE_UNARY(out.hasNbr(Corner<2>::bsw()));
  REQUIRE_EQ(out.getNbrType(Corner<2>::bsw()), NbrType::Normal);
  REQUIRE_EQ(out.getNormalNbrInfo(Corner<2>::bsw()).id, 1);

  REQUIRE_UNARY(!out.hasNbr(Corner<2>::bse()));

  REQUIRE_UNARY(out.hasNbr(Corner<2>::bnw()));
  REQUIRE_EQ(out.getNbrType(Corner<2>::bnw()), NbrType::Fine);
  REQUIRE_EQ(out.getFineNbrInfo(Corner<2>::bnw()).ids[0], 1);

  REQUIRE_UNARY(!out.hasNbr(Corner<2>::bne()));
  REQUIRE_UNARY(!out.hasNbr(Corner<2>::tsw()));

  REQUIRE_UNARY(out.hasNbr(Corner<2>::tse()));
  REQUIRE_EQ(out.getNbrType(Corner<2>::tse()), NbrType::Coarse);
  REQUIRE_EQ(out.getCoarseNbrInfo(Corner<2>::tse()).id, 2);
  REQUIRE_EQ(out.getCoarseNbrInfo(Corner<2>::tse()).orth_on_coarse, Orthant<0>::null());

  REQUIRE_UNARY(!out.hasNbr(Corner<2>::tnw()));
  REQUIRE_UNARY(!out.hasNbr(Corner<2>::tne()));

  // Edges

  REQUIRE_UNARY(!out.hasNbr(Edge::bs()));
  REQUIRE_UNARY(!out.hasNbr(Edge::tn()));

  REQUIRE_UNARY(out.hasNbr(Edge::bn()));
  REQUIRE_EQ(out.getNbrType(Edge::bn()), NbrType::Coarse);
  REQUIRE_EQ(out.getCoarseNbrInfo(Edge::bn()).id, 2);
  REQUIRE_EQ(out.getCoarseNbrInfo(Edge::bn()).orth_on_coarse, Orthant<1>::lower());

  REQUIRE_UNARY(!out.hasNbr(Edge::ts()));
  REQUIRE_UNARY(!out.hasNbr(Edge::bw()));
  REQUIRE_UNARY(!out.hasNbr(Edge::te()));
  REQUIRE_UNARY(!out.hasNbr(Edge::be()));

  REQUIRE_UNARY(out.hasNbr(Edge::tw()));
  REQUIRE_EQ(out.getNbrType(Edge::tw()), NbrType::Fine);
  REQUIRE_EQ(out.getFineNbrInfo(Edge::tw()).ids[0], 1);
  REQUIRE_EQ(out.getFineNbrInfo(Edge::tw()).ids[1], 2);

  REQUIRE_UNARY(out.hasNbr(Edge::sw()));
  REQUIRE_EQ(out.getNbrType(Edge::sw()), NbrType::Normal);
  REQUIRE_EQ(out.getNormalNbrInfo(Edge::sw()).id, 1);

  REQUIRE_UNARY(!out.hasNbr(Edge::ne()));
  REQUIRE_UNARY(!out.hasNbr(Edge::se()));
  REQUIRE_UNARY(!out.hasNbr(Edge::nw()));
}
TEST_CASE("PatchInfo<2> Default Values")
{
  PatchInfo<2> pinfo;
  CHECK_EQ(pinfo.id, 0);
  CHECK_EQ(pinfo.local_index, 0);
  CHECK_EQ(pinfo.global_index, 0);
  CHECK_EQ(pinfo.refine_level, -1);
  CHECK_EQ(pinfo.parent_id, -1);
  CHECK_EQ(pinfo.parent_rank, -1);
  for (int child_id : pinfo.child_ids) {
    CHECK_EQ(child_id, -1);
  }
  for (int child_rank : pinfo.child_ids) {
    CHECK_EQ(child_rank, -1);
  }
  CHECK_EQ(pinfo.num_ghost_cells, 0);
  CHECK_EQ(pinfo.rank, -1);
  CHECK_EQ(pinfo.orth_on_parent, Orthant<2>::null());
  for (int n : pinfo.ns) {
    CHECK_EQ(n, 1);
  }
  for (double start : pinfo.starts) {
    CHECK_EQ(start, 0);
  }
  for (double spacing : pinfo.spacings) {
    CHECK_EQ(spacing, 1);
  }
  CHECK_UNARY_FALSE(pinfo.hasNbr());
}
TEST_CASE("PatchInfo<2> to_json no children")
{
  PatchInfo<2> d;
  d.id = 9;
  d.rank = 0;
  d.parent_id = 2;
  d.parent_rank = 3;
  d.orth_on_parent = Orthant<2>::tnw();
  d.starts = { 1, 2, 3 };
  d.spacings = { 0.1, 0.2, 0.3 };
  d.ns = { 10, 20, 30 };
  d.setNbrInfo(Side<2>::north(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Side<2>::east(), new CoarseNbrInfo<1>(2, Orthant<1>::nw()));
  d.setNbrInfo(Side<2>::south(), new FineNbrInfo<1>({ 3, 4, 5, 6 }));
  d.setNbrInfo(Corner<2>::bsw(), new NormalNbrInfo<0>(1));
  d.setNbrInfo(Corner<2>::tse(), new CoarseNbrInfo<0>(2, Orthant<0>(0)));
  d.setNbrInfo(Corner<2>::bnw(), new FineNbrInfo<0>({ 1 }));
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
TEST_CASE("PatchInfo<2> to_json no children no neighbors")
{
  PatchInfo<2> d;
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
TEST_CASE("PatchInfo<2> to_json with children")
{
  PatchInfo<2> d;
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
  d.setNbrInfo(Side<2>::north(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Side<2>::east(), new CoarseNbrInfo<1>(2, Orthant<1>::nw()));
  d.setNbrInfo(Side<2>::south(), new FineNbrInfo<1>({ 3, 4, 5, 6 }));

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
TEST_CASE("PatchInfo<2> from_json no children")
{
  nlohmann::json j;
  j["id"] = 9;
  j["rank"] = 3;
  j["refine_level"] = 329;
  j["parent_id"] = 2;
  j["parent_rank"] = 3;
  j["starts"] = { 1, 2, 3 };
  j["lengths"] = { 10, 20, 30 };
  j["nbrs"] = { NormalNbrInfo<1>(1), CoarseNbrInfo<1>(2, Orthant<1>::nw()), FineNbrInfo<1>({ 3, 4, 5, 6 }) };
  j["nbrs"][0]["side"] = "NORTH";
  j["nbrs"][1]["side"] = "EAST";
  j["nbrs"][2]["side"] = "SOUTH";
  j["corner_nbrs"] = { NormalNbrInfo<0>(1), CoarseNbrInfo<0>(2, Orthant<0>(0)), FineNbrInfo<0>({ 1 }) };
  j["corner_nbrs"][0]["corner"] = "BSW";
  j["corner_nbrs"][1]["corner"] = "TSE";
  j["corner_nbrs"][2]["corner"] = "BNW";
  j["edge_nbrs"] = { NormalNbrInfo<1>(1), CoarseNbrInfo<1>(2, Orthant<1>::lower()), FineNbrInfo<1>({ 1, 2 }) };
  j["edge_nbrs"][0]["edge"] = "SW";
  j["edge_nbrs"][1]["edge"] = "BN";
  j["edge_nbrs"][2]["edge"] = "TW";

  PatchInfo<2> d = j.get<PatchInfo<2>>();
  CHECK_EQ(d.id, 9);
  CHECK_EQ(d.rank, 3);
  CHECK_EQ(d.refine_level, 329);
  CHECK_EQ(d.parent_id, 2);
  CHECK_EQ(d.parent_rank, 3);
  CHECK_EQ(d.orth_on_parent, Orthant<2>::null());
  CHECK_EQ(d.starts[0], 1);
  CHECK_EQ(d.starts[1], 2);
  CHECK_EQ(d.starts[2], 3);
  CHECK_EQ(d.spacings[0], 10);
  CHECK_EQ(d.spacings[1], 20);
  CHECK_EQ(d.spacings[2], 30);
  CHECK_EQ(d.ns[0], 1);
  CHECK_EQ(d.ns[1], 1);
  CHECK_EQ(d.ns[2], 1);
  CHECK_UNARY_FALSE(d.hasNbr(Side<2>::west()));
  CHECK_UNARY(d.hasNbr(Side<2>::east()));
  CHECK_EQ(d.getNbrType(Side<2>::east()), NbrType::Coarse);
  CHECK_UNARY(d.hasNbr(Side<2>::south()));
  CHECK_EQ(d.getNbrType(Side<2>::south()), NbrType::Fine);
  CHECK_UNARY(d.hasNbr(Side<2>::north()));
  CHECK_EQ(d.getNbrType(Side<2>::north()), NbrType::Normal);
  CHECK_UNARY_FALSE(d.hasNbr(Side<2>::bottom()));
  CHECK_UNARY_FALSE(d.hasNbr(Side<2>::top()));

  CHECK_UNARY(d.hasNbr(Corner<2>::bsw()));
  CHECK_EQ(d.getNbrType(Corner<2>::bsw()), NbrType::Normal);
  CHECK_UNARY_FALSE(d.hasNbr(Corner<2>::bse()));
  CHECK_UNARY(d.hasNbr(Corner<2>::bnw()));
  CHECK_EQ(d.getNbrType(Corner<2>::bnw()), NbrType::Fine);
  CHECK_UNARY_FALSE(d.hasNbr(Corner<2>::bne()));
  CHECK_UNARY_FALSE(d.hasNbr(Corner<2>::tsw()));
  CHECK_UNARY(d.hasNbr(Corner<2>::tse()));
  CHECK_EQ(d.getNbrType(Corner<2>::tse()), NbrType::Coarse);
  CHECK_UNARY_FALSE(d.hasNbr(Corner<2>::tnw()));
  CHECK_UNARY_FALSE(d.hasNbr(Corner<2>::tne()));

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
TEST_CASE("PatchInfo<2> from_json with children")
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
  j["nbrs"] = { NormalNbrInfo<1>(1), CoarseNbrInfo<1>(2, Orthant<1>::nw()), FineNbrInfo<1>({ 3, 4, 5, 6 }) };
  j["nbrs"][0]["side"] = "NORTH";
  j["nbrs"][1]["side"] = "EAST";
  j["nbrs"][2]["side"] = "SOUTH";

  PatchInfo<2> d = j.get<PatchInfo<2>>();
  CHECK_EQ(d.id, 9);
  CHECK_EQ(d.rank, 3);
  CHECK_EQ(d.refine_level, 329);
  CHECK_EQ(d.parent_id, 2);
  CHECK_EQ(d.parent_rank, 3);
  CHECK_EQ(d.orth_on_parent, Orthant<2>::tnw());
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
  CHECK_UNARY_FALSE(d.hasNbr(Side<2>::west()));
  CHECK_UNARY(d.hasNbr(Side<2>::east()));
  CHECK_EQ(d.getNbrType(Side<2>::east()), NbrType::Coarse);
  CHECK_UNARY(d.hasNbr(Side<2>::south()));
  CHECK_EQ(d.getNbrType(Side<2>::south()), NbrType::Fine);
  CHECK_UNARY(d.hasNbr(Side<2>::north()));
  CHECK_EQ(d.getNbrType(Side<2>::north()), NbrType::Normal);
  CHECK_UNARY_FALSE(d.hasNbr(Side<2>::bottom()));
  CHECK_UNARY_FALSE(d.hasNbr(Side<2>::top()));
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
  d.starts = { 1, 2, 3 };
  d.spacings = { 0.1, 0.2, 0.3 };
  d.ns = { 10, 20, 30 };
  d.child_ids = { 3, 4, 5, 6, 7, 8, 9, 10 };
  d.child_ranks = { 1, 2, 3, 4, 5, 6, 7, 8 };
  d.setNbrInfo(Side<2>::north(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Side<2>::east(), new CoarseNbrInfo<1>(2, Orthant<1>::nw()));
  d.setNbrInfo(Side<2>::south(), new FineNbrInfo<1>({ 3, 4, 5, 6 }));
  d.setNbrInfo(Corner<2>::bsw(), new NormalNbrInfo<0>(1));
  d.setNbrInfo(Corner<2>::tse(), new CoarseNbrInfo<0>(2, Orthant<0>(0)));
  d.setNbrInfo(Corner<2>::bnw(), new FineNbrInfo<0>({ 1 }));
  d.setNbrInfo(Edge::sw(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Edge::bn(), new CoarseNbrInfo<1>(2, Orthant<1>::lower()));
  d.setNbrInfo(Edge::tw(), new FineNbrInfo<1>({ 1, 2 }));

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
TEST_CASE("PatchInfo<2> copy assignment")
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
  d.starts = { 1, 2, 3 };
  d.spacings = { 0.1, 0.2, 0.3 };
  d.ns = { 10, 20, 30 };
  d.child_ids = { 3, 4, 5, 6, 7, 8, 9, 10 };
  d.child_ranks = { 1, 2, 3, 4, 5, 6, 7, 8 };
  d.setNbrInfo(Side<2>::north(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Side<2>::east(), new CoarseNbrInfo<1>(2, Orthant<1>::nw()));
  d.setNbrInfo(Side<2>::south(), new FineNbrInfo<1>({ 3, 4, 5, 6 }));
  d.setNbrInfo(Corner<2>::bsw(), new NormalNbrInfo<0>(1));
  d.setNbrInfo(Corner<2>::tse(), new CoarseNbrInfo<0>(2, Orthant<0>(0)));
  d.setNbrInfo(Corner<2>::bnw(), new FineNbrInfo<0>({ 1 }));
  d.setNbrInfo(Edge::sw(), new NormalNbrInfo<1>(1));
  d.setNbrInfo(Edge::bn(), new CoarseNbrInfo<1>(2, Orthant<1>::lower()));
  d.setNbrInfo(Edge::tw(), new FineNbrInfo<1>({ 1, 2 }));

  PatchInfo<2> d2;
  d2 = d;

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