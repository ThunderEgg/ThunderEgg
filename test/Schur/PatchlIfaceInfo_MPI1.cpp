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
#include <ThunderEgg/Schur/PatchIfaceInfo.h>
#include <algorithm>
using namespace std;
using namespace ThunderEgg;

TEST_CASE("Schur::PatchIfaceInfo default constructor", "[Schur::PatchIfaceInfo]")
{
	Schur::PatchIfaceInfo<2> piinfo;
	CHECK(piinfo.pinfo == nullptr);
	for (Side<2> s : Side<2>::getValues()) {
		CHECK(piinfo.iface_info[s.getIndex()] == nullptr);
		CHECK(piinfo.getIfaceInfo(s) == nullptr);
		CHECK(piinfo.getNormalIfaceInfo(s) == nullptr);
		CHECK(piinfo.getCoarseIfaceInfo(s) == nullptr);
		CHECK(piinfo.getFineIfaceInfo(s) == nullptr);
	}
}
TEST_CASE("Schur::PatchIfaceInfo setIfaceInfo with NormalIfaceInfo", "[Schur::PatchIfaceInfo]")
{
	for (Side<2> side_to_set : Side<2>::getValues()) {
		// setup
		int  id                                   = 1;
		int  nbr_id                               = 2;
		auto pinfo                                = make_shared<PatchInfo<2>>();
		pinfo->rank                               = 0;
		pinfo->id                                 = id;
		pinfo->nbr_info[side_to_set.getIndex()]   = make_shared<NormalNbrInfo<2>>(nbr_id);
		pinfo->getNormalNbrInfo(side_to_set).rank = 1;
		auto iface_info = make_shared<Schur::NormalIfaceInfo<2>>(pinfo, side_to_set);

		Schur::PatchIfaceInfo<2> piinfo;
		piinfo.setIfaceInfo(side_to_set, iface_info);
		CHECK(piinfo.pinfo == nullptr);
		for (Side<2> s : Side<2>::getValues()) {
			if (s == side_to_set) {
				CHECK(piinfo.iface_info[s.getIndex()] == iface_info);
				CHECK(piinfo.getIfaceInfo(s) == iface_info);
				CHECK(piinfo.getNormalIfaceInfo(s) == iface_info);
			} else {
				CHECK(piinfo.iface_info[s.getIndex()] == nullptr);
				CHECK(piinfo.getIfaceInfo(s) == nullptr);
				CHECK(piinfo.getNormalIfaceInfo(s) == nullptr);
			}
			CHECK(piinfo.getCoarseIfaceInfo(s) == nullptr);
			CHECK(piinfo.getFineIfaceInfo(s) == nullptr);
		}
	}
}
TEST_CASE("Schur::PatchIfaceInfo setIfaceInfo with FineIfaceInfo", "[Schur::PatchIfaceInfo]")
{
	for (Side<2> side_to_set : Side<2>::getValues()) {
		// setup
		int           id                            = 1;
		array<int, 2> nbr_ids                       = {2, 3};
		auto          pinfo                         = make_shared<PatchInfo<2>>();
		pinfo->id                                   = id;
		pinfo->nbr_info[side_to_set.getIndex()]     = make_shared<FineNbrInfo<2>>(nbr_ids);
		pinfo->getFineNbrInfo(side_to_set).ranks[0] = 1;
		pinfo->getFineNbrInfo(side_to_set).ranks[1] = 2;
		auto iface_info = make_shared<Schur::FineIfaceInfo<2>>(pinfo, side_to_set);

		Schur::PatchIfaceInfo<2> piinfo;
		piinfo.setIfaceInfo(side_to_set, iface_info);
		CHECK(piinfo.pinfo == nullptr);
		for (Side<2> s : Side<2>::getValues()) {
			if (s == side_to_set) {
				CHECK(piinfo.iface_info[s.getIndex()] == iface_info);
				CHECK(piinfo.getIfaceInfo(s) == iface_info);
				CHECK(piinfo.getFineIfaceInfo(s) == iface_info);
			} else {
				CHECK(piinfo.iface_info[s.getIndex()] == nullptr);
				CHECK(piinfo.getIfaceInfo(s) == nullptr);
				CHECK(piinfo.getFineIfaceInfo(s) == nullptr);
			}
			CHECK(piinfo.getNormalIfaceInfo(s) == nullptr);
			CHECK(piinfo.getCoarseIfaceInfo(s) == nullptr);
		}
	}
}
TEST_CASE("Schur::PatchIfaceInfo setIfaceInfo with CoarseIfaceInfo", "[Schur::PatchIfaceInfo]")
{
	for (Side<2> side_to_set : Side<2>::getValues()) {
		// setup
		int  id     = 1;
		int  nbr_id = 2;
		auto pinfo  = make_shared<PatchInfo<2>>();
		pinfo->id   = id;
		pinfo->nbr_info[side_to_set.getIndex()]
		= make_shared<CoarseNbrInfo<2>>(nbr_id, Orthant<1>::upper());
		pinfo->getCoarseNbrInfo(side_to_set).rank = 1;
		auto iface_info = make_shared<Schur::CoarseIfaceInfo<2>>(pinfo, side_to_set);

		Schur::PatchIfaceInfo<2> piinfo;
		piinfo.setIfaceInfo(side_to_set, iface_info);
		CHECK(piinfo.pinfo == nullptr);
		for (Side<2> s : Side<2>::getValues()) {
			if (s == side_to_set) {
				CHECK(piinfo.iface_info[s.getIndex()] == iface_info);
				CHECK(piinfo.getIfaceInfo(s) == iface_info);
				CHECK(piinfo.getCoarseIfaceInfo(s) == iface_info);
			} else {
				CHECK(piinfo.iface_info[s.getIndex()] == nullptr);
				CHECK(piinfo.getIfaceInfo(s) == nullptr);
				CHECK(piinfo.getCoarseIfaceInfo(s) == nullptr);
			}
			CHECK(piinfo.getNormalIfaceInfo(s) == nullptr);
			CHECK(piinfo.getFineIfaceInfo(s) == nullptr);
		}
	}
}
TEST_CASE("Schur::PatchIfaceInfo PatchInfo constructor", "[Schur::PatchIfaceInfo]")
{
	// setup
	int           id            = 1;
	int           nbr_id        = 2;
	array<int, 2> fine_nbr_ids  = {3, 4};
	int           coarse_nbr_id = 5;
	auto          pinfo         = make_shared<PatchInfo<2>>();
	pinfo->id                   = id;
	pinfo->nbr_info[Side<2>::west().getIndex()]
	= make_shared<CoarseNbrInfo<2>>(coarse_nbr_id, Orthant<1>::upper());
	pinfo->nbr_info[Side<2>::south().getIndex()] = make_shared<FineNbrInfo<2>>(fine_nbr_ids);
	pinfo->nbr_info[Side<2>::north().getIndex()] = make_shared<NormalNbrInfo<2>>(nbr_id);

	Schur::PatchIfaceInfo<2> piinfo(pinfo);

	CHECK(piinfo.pinfo == pinfo);

	CHECK(piinfo.iface_info[Side<2>::west().getIndex()] != nullptr);
	CHECK(piinfo.getIfaceInfo(Side<2>::west()) != nullptr);
	CHECK(piinfo.getNormalIfaceInfo(Side<2>::west()) == nullptr);
	CHECK(piinfo.getCoarseIfaceInfo(Side<2>::west()) != nullptr);
	CHECK(piinfo.getFineIfaceInfo(Side<2>::west()) == nullptr);

	CHECK(piinfo.iface_info[Side<2>::east().getIndex()] == nullptr);
	CHECK(piinfo.getIfaceInfo(Side<2>::east()) == nullptr);
	CHECK(piinfo.getNormalIfaceInfo(Side<2>::east()) == nullptr);
	CHECK(piinfo.getCoarseIfaceInfo(Side<2>::east()) == nullptr);
	CHECK(piinfo.getFineIfaceInfo(Side<2>::east()) == nullptr);

	CHECK(piinfo.iface_info[Side<2>::south().getIndex()] != nullptr);
	CHECK(piinfo.getIfaceInfo(Side<2>::south()) != nullptr);
	CHECK(piinfo.getNormalIfaceInfo(Side<2>::south()) == nullptr);
	CHECK(piinfo.getCoarseIfaceInfo(Side<2>::south()) == nullptr);
	CHECK(piinfo.getFineIfaceInfo(Side<2>::south()) != nullptr);

	CHECK(piinfo.iface_info[Side<2>::north().getIndex()] != nullptr);
	CHECK(piinfo.getIfaceInfo(Side<2>::north()) != nullptr);
	CHECK(piinfo.getNormalIfaceInfo(Side<2>::north()) != nullptr);
	CHECK(piinfo.getCoarseIfaceInfo(Side<2>::north()) == nullptr);
	CHECK(piinfo.getFineIfaceInfo(Side<2>::north()) == nullptr);
}
TEST_CASE("Schur::PatchIfaceInfo setLocalIndexesFromId", "[Schur::PatchIfaceInfo]")
{
	for (Side<2> side_to_set : Side<2>::getValues()) {
		// setup
		int  id                                   = 1;
		int  nbr_id                               = 2;
		auto pinfo                                = make_shared<PatchInfo<2>>();
		pinfo->rank                               = 0;
		pinfo->id                                 = id;
		pinfo->nbr_info[side_to_set.getIndex()]   = make_shared<NormalNbrInfo<2>>(nbr_id);
		pinfo->getNormalNbrInfo(side_to_set).rank = 1;
		auto iface_info = make_shared<Schur::NormalIfaceInfo<2>>(pinfo, side_to_set);

		Schur::PatchIfaceInfo<2> piinfo;
		piinfo.setIfaceInfo(side_to_set, iface_info);

		int           local_index = 1;
		map<int, int> id_local_index_map;
		id_local_index_map[piinfo.getIfaceInfo(side_to_set)->id] = local_index;

		piinfo.setLocalIndexesFromId(id_local_index_map);

		CHECK(piinfo.getIfaceInfo(side_to_set)->local_index == local_index);
	}
}
TEST_CASE("Schur::PatchIfaceInfo setGlobalIndexesFormLocalIndex", "[Schur::PatchIfaceInfo]")
{
	for (Side<2> side_to_set : Side<2>::getValues()) {
		// setup
		int  id                                   = 1;
		int  nbr_id                               = 2;
		auto pinfo                                = make_shared<PatchInfo<2>>();
		pinfo->rank                               = 0;
		pinfo->id                                 = id;
		pinfo->nbr_info[side_to_set.getIndex()]   = make_shared<NormalNbrInfo<2>>(nbr_id);
		pinfo->getNormalNbrInfo(side_to_set).rank = 1;
		auto iface_info = make_shared<Schur::NormalIfaceInfo<2>>(pinfo, side_to_set);

		Schur::PatchIfaceInfo<2> piinfo;
		piinfo.setIfaceInfo(side_to_set, iface_info);

		int           local_index = 1;
		map<int, int> id_local_index_map;
		id_local_index_map[piinfo.getIfaceInfo(side_to_set)->id] = local_index;
		piinfo.setLocalIndexesFromId(id_local_index_map);

		map<int, int> local_global_index_map;
		local_global_index_map[local_index] = 1;
		piinfo.setGlobalIndexesFromLocalIndex(local_global_index_map);

		CHECK(piinfo.getIfaceInfo(side_to_set)->global_index == 1);
	}
}
TEST_CASE("Schur::PatchIfaceInfo serialization", "[Schur::PatchIfaceInfo]")
{
	// setup
	int           id            = 1;
	int           nbr_id        = 2;
	array<int, 2> fine_nbr_ids  = {3, 4};
	int           coarse_nbr_id = 5;
	auto          pinfo         = make_shared<PatchInfo<2>>();
	pinfo->id                   = id;
	pinfo->nbr_info[Side<2>::west().getIndex()]
	= make_shared<CoarseNbrInfo<2>>(coarse_nbr_id, Orthant<1>::upper());
	pinfo->nbr_info[Side<2>::south().getIndex()] = make_shared<FineNbrInfo<2>>(fine_nbr_ids);
	pinfo->nbr_info[Side<2>::north().getIndex()] = make_shared<NormalNbrInfo<2>>(nbr_id);

	Schur::PatchIfaceInfo<2> piinfo(pinfo);

	char *buff = new char[piinfo.serialize(nullptr)];
	piinfo.serialize(buff);
	Schur::PatchIfaceInfo<2> out;
	out.deserialize(buff);

	CHECK(out.pinfo->id == pinfo->id);

	CHECK(out.iface_info[Side<2>::west().getIndex()] != nullptr);
	CHECK(out.getIfaceInfo(Side<2>::west()) != nullptr);
	CHECK(out.getNormalIfaceInfo(Side<2>::west()) == nullptr);
	CHECK(out.getCoarseIfaceInfo(Side<2>::west()) != nullptr);
	CHECK(out.getFineIfaceInfo(Side<2>::west()) == nullptr);

	CHECK(out.iface_info[Side<2>::east().getIndex()] == nullptr);
	CHECK(out.getIfaceInfo(Side<2>::east()) == nullptr);
	CHECK(out.getNormalIfaceInfo(Side<2>::east()) == nullptr);
	CHECK(out.getCoarseIfaceInfo(Side<2>::east()) == nullptr);
	CHECK(out.getFineIfaceInfo(Side<2>::east()) == nullptr);

	CHECK(out.iface_info[Side<2>::south().getIndex()] != nullptr);
	CHECK(out.getIfaceInfo(Side<2>::south()) != nullptr);
	CHECK(out.getNormalIfaceInfo(Side<2>::south()) == nullptr);
	CHECK(out.getCoarseIfaceInfo(Side<2>::south()) == nullptr);
	CHECK(out.getFineIfaceInfo(Side<2>::south()) != nullptr);

	CHECK(out.iface_info[Side<2>::north().getIndex()] != nullptr);
	CHECK(out.getIfaceInfo(Side<2>::north()) != nullptr);
	CHECK(out.getNormalIfaceInfo(Side<2>::north()) != nullptr);
	CHECK(out.getCoarseIfaceInfo(Side<2>::north()) == nullptr);
	CHECK(out.getFineIfaceInfo(Side<2>::north()) == nullptr);
}