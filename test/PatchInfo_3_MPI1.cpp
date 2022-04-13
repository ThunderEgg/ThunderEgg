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

TEST_CASE("PatchInfo<3> getNbrIds NormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<2>* nbr_info = new NormalNbrInfo<2>(2);
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 2);
}
TEST_CASE("PatchInfo<3> getNbrRanks NormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<2>* nbr_info = new NormalNbrInfo<2>(2);
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes exists NormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<2>* nbr_info = new NormalNbrInfo<2>(2);
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, 30);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes does not exist NormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<2>* nbr_info = new NormalNbrInfo<2>(2);
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, -1);
}
TEST_CASE("PatchInfo<3> setNeighborGlobalIndexes exists NormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<2>* nbr_info = new NormalNbrInfo<2>(2);
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_index, 30);
}
TEST_CASE("PatchInfo<3> getNbrIds CoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<2>* nbr_info = new CoarseNbrInfo<2>(2, Orthant<2>::sw());
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 2);
}
TEST_CASE("PatchInfo<3> getNbrRanks CoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<2>* nbr_info = new CoarseNbrInfo<2>(2, Orthant<2>::sw());
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes exists CoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<2>* nbr_info = new CoarseNbrInfo<2>(2, Orthant<2>::sw());
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, 30);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes does not exist CoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<2>* nbr_info = new CoarseNbrInfo<2>(2, Orthant<2>::sw());
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, -1);
}
TEST_CASE("PatchInfo<3> setNeighborGlobalIndexes exists CoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<2>* nbr_info = new CoarseNbrInfo<2>(2, Orthant<2>::sw());
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_index, 30);
}
TEST_CASE("PatchInfo<3> getNbrIds FineNbrInfo")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<2>* nbr_info = new FineNbrInfo<2>({ 1, 2, 3, 4 });
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 4);
  CHECK_EQ(ids[0], 1);
  CHECK_EQ(ids[1], 2);
  CHECK_EQ(ids[2], 3);
  CHECK_EQ(ids[3], 4);
}
TEST_CASE("PatchInfo<3> getNbrRanks FineNbrInfo")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<2>* nbr_info = new FineNbrInfo<2>({ 1, 2, 3, 4 });
  nbr_info->ranks = { 3, 4, 6, 7 };
  pinfo.setNbrInfo(Side<3>::west(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 4);
  CHECK_EQ(ranks[0], 3);
  CHECK_EQ(ranks[1], 4);
  CHECK_EQ(ranks[2], 6);
  CHECK_EQ(ranks[3], 7);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes exists FineNbrInfo")
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

  CHECK_EQ(nbr_info->local_indexes[0], 30);
  CHECK_EQ(nbr_info->local_indexes[1], 31);
  CHECK_EQ(nbr_info->local_indexes[2], 32);
  CHECK_EQ(nbr_info->local_indexes[3], 33);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes does not exist FineNbrInfo")
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

  CHECK_EQ(nbr_info->local_indexes[0], 30);
  CHECK_EQ(nbr_info->local_indexes[1], 31);
  CHECK_EQ(nbr_info->local_indexes[2], 32);
  CHECK_EQ(nbr_info->local_indexes[3], -1);
}
TEST_CASE("PatchInfo<3> setNeighborGlobalIndexes exists FineNbrInfo")
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

  CHECK_EQ(nbr_info->global_indexes[0], 30);
  CHECK_EQ(nbr_info->global_indexes[1], 31);
  CHECK_EQ(nbr_info->global_indexes[2], 32);
  CHECK_EQ(nbr_info->global_indexes[3], 33);
}
TEST_CASE("PatchInfo<3> getNbrIds EdgeNormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 2);
}
TEST_CASE("PatchInfo<3> getNbrRanks EdgeNormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes exists EdgeNormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, 30);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes does not exist EdgeNormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, -1);
}
TEST_CASE("PatchInfo<3> setNeighborGlobalIndexes exists EdgeNormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_index, 30);
}
TEST_CASE("PatchInfo<3> getNbrIds EdgeCoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 2);
}
TEST_CASE("PatchInfo<3> getNbrRanks EdgeCoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes exists EdgeCoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, 30);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes does not exist EdgeCoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, -1);
}
TEST_CASE("PatchInfo<3> setNeighborGlobalIndexes exists EdgeCoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_index, 30);
}
TEST_CASE("PatchInfo<3> getNbrIds EdgeFineNbrInfo")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 1, 2 });
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 2);
  CHECK_EQ(ids[0], 1);
  CHECK_EQ(ids[1], 2);
}
TEST_CASE("PatchInfo<3> getNbrRanks EdgeFineNbrInfo")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 1, 2 });
  nbr_info->ranks = { 3, 4 };
  pinfo.setNbrInfo(Edge::bs(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 2);
  CHECK_EQ(ranks[0], 3);
  CHECK_EQ(ranks[1], 4);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes exists EdgeFineNbrInfo")
{
  PatchInfo<3> pinfo;
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
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes does not exist EdgeFineNbrInfo")
{
  PatchInfo<3> pinfo;
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
TEST_CASE("PatchInfo<3> setNeighborGlobalIndexes exists EdgeFineNbrInfo")
{
  PatchInfo<3> pinfo;
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
TEST_CASE("PatchInfo<3> getNbrIds CornerNormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 2);
}
TEST_CASE("PatchInfo<3> getNbrRanks CornerNormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes exists CornerNormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, 30);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes does not exist CornerNormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, -1);
}
TEST_CASE("PatchInfo<3> setNeighborGlobalIndexes exists CornerNormalNbrInfo")
{
  PatchInfo<3> pinfo;
  NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_index, 30);
}
TEST_CASE("PatchInfo<3> getNbrIds CornerCoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 2);
}
TEST_CASE("PatchInfo<3> getNbrRanks CornerCoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  nbr_info->rank = 3;
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes exists CornerCoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, 30);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes does not exist CornerCoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_index, -1);
}
TEST_CASE("PatchInfo<3> setNeighborGlobalIndexes exists CornerCoarseNbrInfo")
{
  PatchInfo<3> pinfo;
  CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_index, 30);
}
TEST_CASE("PatchInfo<3> getNbrIds CornerFineNbrInfo")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 1 });
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ids = pinfo.getNbrIds();

  REQUIRE_EQ(ids.size(), 1);
  CHECK_EQ(ids[0], 1);
}
TEST_CASE("PatchInfo<3> getNbrRanks CornerFineNbrInfo")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 1 });
  nbr_info->ranks = { 3 };
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();

  REQUIRE_EQ(ranks.size(), 1);
  CHECK_EQ(ranks[0], 3);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes exists CornerFineNbrInfo")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 2 });
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[2] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_indexes[0], 30);
}
TEST_CASE("PatchInfo<3> setNeighborLocalIndexes does not exist CornerFineNbrInfo")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 2 });
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_local_index_map;
  id_to_local_index_map[59] = 30;
  pinfo.setNeighborLocalIndexes(id_to_local_index_map);

  CHECK_EQ(nbr_info->local_indexes[0], -1);
}
TEST_CASE("PatchInfo<3> setNeighborGlobalIndexes exists CornerFineNbrInfo")
{
  PatchInfo<3> pinfo;
  FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 2 });
  pinfo.setNbrInfo(Corner<3>::bsw(), nbr_info);
  deque<int> ranks = pinfo.getNbrRanks();
  map<int, int> id_to_global_index_map;
  id_to_global_index_map[2] = 30;
  pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

  CHECK_EQ(nbr_info->global_indexes[0], 30);
}
TEST_CASE("PatchInfo<3> Serialization/Deserialization")
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
  REQUIRE_EQ(out.id, 0);

  REQUIRE_UNARY(!out.hasNbr(Side<3>::west()));

  REQUIRE_UNARY(out.hasNbr(Side<3>::east()));
  REQUIRE_EQ(out.getNbrType(Side<3>::east()), NbrType::Coarse);
  REQUIRE_EQ(out.getCoarseNbrInfo(Side<3>::east()).id, 2);
  REQUIRE_EQ(out.getCoarseNbrInfo(Side<3>::east()).orth_on_coarse, Orthant<2>::nw());

  REQUIRE_UNARY(out.hasNbr(Side<3>::south()));
  REQUIRE_EQ(out.getNbrType(Side<3>::south()), NbrType::Fine);
  REQUIRE_EQ(out.getFineNbrInfo(Side<3>::south()).ids[0], 3);
  REQUIRE_EQ(out.getFineNbrInfo(Side<3>::south()).ids[1], 4);
  REQUIRE_EQ(out.getFineNbrInfo(Side<3>::south()).ids[2], 5);
  REQUIRE_EQ(out.getFineNbrInfo(Side<3>::south()).ids[3], 6);

  REQUIRE_UNARY(out.hasNbr(Side<3>::north()));
  REQUIRE_EQ(out.getNbrType(Side<3>::north()), NbrType::Normal);
  REQUIRE_EQ(out.getNormalNbrInfo(Side<3>::north()).id, 1);

  REQUIRE_UNARY(!out.hasNbr(Side<3>::bottom()));
  REQUIRE_UNARY(!out.hasNbr(Side<3>::top()));

  // Corners

  REQUIRE_UNARY(out.hasNbr(Corner<3>::bsw()));
  REQUIRE_EQ(out.getNbrType(Corner<3>::bsw()), NbrType::Normal);
  REQUIRE_EQ(out.getNormalNbrInfo(Corner<3>::bsw()).id, 1);

  REQUIRE_UNARY(!out.hasNbr(Corner<3>::bse()));

  REQUIRE_UNARY(out.hasNbr(Corner<3>::bnw()));
  REQUIRE_EQ(out.getNbrType(Corner<3>::bnw()), NbrType::Fine);
  REQUIRE_EQ(out.getFineNbrInfo(Corner<3>::bnw()).ids[0], 1);

  REQUIRE_UNARY(!out.hasNbr(Corner<3>::bne()));
  REQUIRE_UNARY(!out.hasNbr(Corner<3>::tsw()));

  REQUIRE_UNARY(out.hasNbr(Corner<3>::tse()));
  REQUIRE_EQ(out.getNbrType(Corner<3>::tse()), NbrType::Coarse);
  REQUIRE_EQ(out.getCoarseNbrInfo(Corner<3>::tse()).id, 2);
  REQUIRE_EQ(out.getCoarseNbrInfo(Corner<3>::tse()).orth_on_coarse, Orthant<0>::null());

  REQUIRE_UNARY(!out.hasNbr(Corner<3>::tnw()));
  REQUIRE_UNARY(!out.hasNbr(Corner<3>::tne()));

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
TEST_CASE("PatchInfo<3> Default Values")
{
  PatchInfo<3> pinfo;
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
  CHECK_EQ(pinfo.orth_on_parent, Orthant<3>::null());
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

TEST_CASE("PatchInfo<3> copy constructor")
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

  for (Side<3> s : Side<3>::getValues()) {
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
  for (Corner<3> c : Corner<3>::getValues()) {
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
  for (Edge c : Edge::getValues()) {
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
TEST_CASE("PatchInfo<3> copy assignment")
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

  for (Side<3> s : Side<3>::getValues()) {
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
  for (Corner<3> c : Corner<3>::getValues()) {
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
  for (Edge c : Edge::getValues()) {
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
