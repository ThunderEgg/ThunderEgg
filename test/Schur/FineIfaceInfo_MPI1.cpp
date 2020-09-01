/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
 *  top-level directory.
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

#include "catch.hpp"
#include <ThunderEgg/Schur/FineIfaceInfo.h>
#include <algorithm>
using namespace std;
using namespace ThunderEgg;

template <typename Container, typename Value> bool contains(Container &deque, Value a)
{
	return find(deque.begin(), deque.end(), a) != deque.end();
}
TEST_CASE("Schur::FineIfaceInfo constructor", "[Schur::FineIfaceInfo]")
{
	for (Side<2> s : Side<2>::getValues()) {
		int           id                  = 1;
		array<int, 2> nbr_ids             = {2, 3};
		auto          pinfo               = make_shared<PatchInfo<2>>();
		pinfo->rank                       = 0;
		pinfo->id                         = id;
		pinfo->nbr_info[s.getIndex()]     = make_shared<FineNbrInfo<2>>(nbr_ids);
		pinfo->getFineNbrInfo(s).ranks[0] = 1;
		pinfo->getFineNbrInfo(s).ranks[1] = 2;
		Schur::FineIfaceInfo<2> iface_info(pinfo, s);
		INFO("Side: " << s);
		CHECK(iface_info.rank == 0);
		CHECK(iface_info.fine_ranks[0] == 1);
		CHECK(iface_info.fine_ranks[1] == 2);
		// check that the id is encoded as expected
		CHECK(iface_info.id / (int) Side<2>::num_sides == id);
		CHECK(iface_info.id % Side<2>::num_sides == s.getIndex());
		// check that iface belongs to nbr
		CHECK(iface_info.fine_ids[0] / (int) Side<2>::num_sides == nbr_ids[0]);
		CHECK(iface_info.fine_ids[0] % Side<2>::num_sides == s.opposite().getIndex());
		CHECK(iface_info.fine_ids[1] / (int) Side<2>::num_sides == nbr_ids[1]);
		CHECK(iface_info.fine_ids[1] % Side<2>::num_sides == s.opposite().getIndex());
		// local and global index should be set to -1
		CHECK(iface_info.local_index == -1);
		CHECK(iface_info.fine_local_indexes[0] == -1);
		CHECK(iface_info.fine_local_indexes[1] == -1);
		CHECK(iface_info.global_index == -1);
		CHECK(iface_info.fine_global_indexes[0] == -1);
		CHECK(iface_info.fine_global_indexes[1] == -1);
		CHECK(iface_info.nbr_info == pinfo->nbr_info[s.getIndex()]);
	}
}
TEST_CASE("Schur::FineIfaceInfo setLocalIndexesFromId", "[Schur::FineIfaceInfo]")
{
	for (Side<2> s : Side<2>::getValues()) {
		int           id                  = 1;
		array<int, 2> nbr_ids             = {2, 3};
		int           local_index         = GENERATE(0, 1, 2);
		auto          pinfo               = make_shared<PatchInfo<2>>();
		pinfo->id                         = id;
		pinfo->nbr_info[s.getIndex()]     = make_shared<FineNbrInfo<2>>(nbr_ids);
		pinfo->getFineNbrInfo(s).ranks[0] = 1;
		pinfo->getFineNbrInfo(s).ranks[1] = 2;
		Schur::FineIfaceInfo<2> iface_info(pinfo, s);
		// local and global index should be set to -1

		map<int, int> id_local_index_map;
		id_local_index_map[iface_info.id]          = local_index;
		id_local_index_map[iface_info.fine_ids[0]] = local_index + 1;
		id_local_index_map[iface_info.fine_ids[1]] = local_index + 2;
		iface_info.setLocalIndexesFromId(id_local_index_map);
		CHECK(iface_info.local_index == local_index);
		CHECK(iface_info.fine_local_indexes[0] == local_index + 1);
		CHECK(iface_info.fine_local_indexes[1] == local_index + 2);
	}
}
TEST_CASE("Schur::FineIfaceInfo setGlobalIndexesFromLocalIndex", "[Schur::FineIfaceInfo]")
{
	for (Side<2> s : Side<2>::getValues()) {
		int           id                  = 1;
		array<int, 2> nbr_ids             = {2, 3};
		int           local_index         = 0;
		int           global_index        = GENERATE(0, 1, 2);
		auto          pinfo               = make_shared<PatchInfo<2>>();
		pinfo->id                         = id;
		pinfo->nbr_info[s.getIndex()]     = make_shared<FineNbrInfo<2>>(nbr_ids);
		pinfo->getFineNbrInfo(s).ranks[0] = 1;
		pinfo->getFineNbrInfo(s).ranks[1] = 2;
		Schur::FineIfaceInfo<2> iface_info(pinfo, s);
		// local and global index should be set to -1

		map<int, int> id_local_index_map;
		id_local_index_map[iface_info.id]          = local_index;
		id_local_index_map[iface_info.fine_ids[0]] = local_index + 1;
		id_local_index_map[iface_info.fine_ids[1]] = local_index + 2;
		iface_info.setLocalIndexesFromId(id_local_index_map);

		map<int, int> local_global_index_map;
		local_global_index_map[local_index]     = global_index;
		local_global_index_map[local_index + 1] = global_index + 1;
		local_global_index_map[local_index + 2] = global_index + 2;
		iface_info.setGlobalIndexesFromLocalIndex(local_global_index_map);
		CHECK(iface_info.local_index == local_index);
		CHECK(iface_info.fine_global_indexes[0] == global_index + 1);
		CHECK(iface_info.fine_global_indexes[1] == global_index + 2);
	}
}