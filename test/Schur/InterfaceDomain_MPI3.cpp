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

#include "../utils/DomainReader.h"
#include "catch.hpp"
#include <ThunderEgg/Schur/InterfaceDomain.h>
#include <limits>
using namespace std;
using namespace ThunderEgg;
TEST_CASE(
"Schur::InterfaceDomain<2> patch interface local indexes start from 0 2d_refined_complicated",
"[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	int min_local_index = numeric_limits<int>::max();

	for (auto patch : interface_domain.getPatchIfaceInfos()) {
		for (Side<2> s : Side<2>::getValues()) {
			if (patch->pinfo->hasNbr(s)) {
				min_local_index = min(min_local_index, patch->getIfaceInfo(s)->patch_local_index);
			}
		}
	}

	CHECK(min_local_index == 0);
}
TEST_CASE(
"Schur::InterfaceDomain<2> patch interface local indexes are contiguous 2d_refined_complicated",
"[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	set<int> local_indexes;

	for (auto patch : interface_domain.getPatchIfaceInfos()) {
		for (Side<2> s : Side<2>::getValues()) {
			if (patch->pinfo->hasNbr(s)) {
				local_indexes.insert(patch->getIfaceInfo(s)->patch_local_index);
			}
		}
	}

	int prev_local_index = -1;
	for (int local_index : local_indexes) {
		REQUIRE(local_index == prev_local_index + 1);
		prev_local_index = local_index;
	}
}
TEST_CASE(
"Schur::InterfaceDomain<2> each id has only one patch interface local index associated with it 2d_refined_complicated",
"[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	map<int, set<int>> id_to_local_indexes;

	for (auto piinfo : interface_domain.getPatchIfaceInfos()) {
		for (Side<2> s : Side<2>::getValues()) {
			if (piinfo->pinfo->hasNbr(s)) {
				id_to_local_indexes[piinfo->getIfaceInfo(s)->id].insert(
				piinfo->getIfaceInfo(s)->patch_local_index);
			}
		}
	}

	REQUIRE(id_to_local_indexes.size() > 0);
	for (auto pair : id_to_local_indexes) {
		INFO("ID " << pair.first);
		CHECK(pair.second.size() == 1);
	}
}
TEST_CASE(
"Schur::InterfaceDomain<2> each patch interface local index has only one id associated with it 2d_refined_complicated",
"[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	map<int, set<int>> local_index_to_ids;

	for (auto piinfo : interface_domain.getPatchIfaceInfos()) {
		for (Side<2> s : Side<2>::getValues()) {
			if (piinfo->pinfo->hasNbr(s)) {
				local_index_to_ids[piinfo->getIfaceInfo(s)->patch_local_index].insert(
				piinfo->getIfaceInfo(s)->id);
			}
		}
	}

	REQUIRE(local_index_to_ids.size() > 0);
	for (auto pair : local_index_to_ids) {
		INFO("Local Index " << pair.first);
		CHECK(pair.second.size() == 1);
	}
}
TEST_CASE("Schur::InterfaceDomain<2> row local indexes start from 0 2d_refined_complicated",
          "[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	int min_local_index = numeric_limits<int>::max();

	for (auto iface : interface_domain.getInterfaces()) {
		min_local_index = min(min_local_index, iface->local_index);

		for (auto patch : iface->patches) {
			for (Side<2> s : Side<2>::getValues()) {
				if (patch.piinfo->pinfo->hasNbr(s)) {
					min_local_index
					= min(min_local_index, patch.piinfo->getIfaceInfo(s)->row_local_index);
				}
			}
		}
	}

	CHECK(min_local_index == 0);
}
TEST_CASE("Schur::InterfaceDomain<2> row local indexes are contiguous 2d_refined_complicated",
          "[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	set<int> local_indexes;

	for (auto iface : interface_domain.getInterfaces()) {
		local_indexes.insert(iface->local_index);

		for (auto patch : iface->patches) {
			for (Side<2> s : Side<2>::getValues()) {
				if (patch.piinfo->pinfo->hasNbr(s)) {
					local_indexes.insert(patch.piinfo->getIfaceInfo(s)->row_local_index);
				}
			}
		}
	}

	int prev_local_index = -1;
	for (int local_index : local_indexes) {
		REQUIRE(local_index == prev_local_index + 1);
		prev_local_index = local_index;
	}
}
TEST_CASE(
"Schur::InterfaceDomain<2> each id has only one row local index associated with it 2d_refined_complicated",
"[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	map<int, set<int>> id_to_local_indexes;

	for (auto iface : interface_domain.getInterfaces()) {
		id_to_local_indexes[iface->id].insert(iface->local_index);

		for (auto patch : iface->patches) {
			for (Side<2> s : Side<2>::getValues()) {
				if (patch.piinfo->pinfo->hasNbr(s)) {
					auto iface_info = patch.piinfo->getIfaceInfo(s);
					id_to_local_indexes[iface_info->id].insert(iface_info->row_local_index);
				}
			}
		}
	}

	REQUIRE(id_to_local_indexes.size() > 0);
	for (auto pair : id_to_local_indexes) {
		INFO("ID " << pair.first);
		CHECK(pair.second.size() == 1);
	}
}
TEST_CASE(
"Schur::InterfaceDomain<2> each row local index has only one id associated with it 2d_refined_complicated",
"[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	map<int, set<int>> local_index_to_ids;

	for (auto iface : interface_domain.getInterfaces()) {
		local_index_to_ids[iface->local_index].insert(iface->id);

		for (auto patch : iface->patches) {
			for (Side<2> s : Side<2>::getValues()) {
				if (patch.piinfo->pinfo->hasNbr(s)) {
					auto iface_info = patch.piinfo->getIfaceInfo(s);
					local_index_to_ids[iface_info->row_local_index].insert(iface_info->id);
				}
			}
		}
	}

	REQUIRE(local_index_to_ids.size() > 0);
	for (auto pair : local_index_to_ids) {
		INFO("Local Index " << pair.first);
		CHECK(pair.second.size() == 1);
	}
}
TEST_CASE("Schur::InterfaceDomain<2> col local indexes start from 0 2d_refined_complicated",
          "[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	int min_local_index = numeric_limits<int>::max();

	for (auto iface : interface_domain.getInterfaces()) {
		min_local_index = min(min_local_index, iface->local_index);

		for (auto patch : iface->patches) {
			if (patch.type.isNormal() || patch.type.isFineToFine()
			    || patch.type.isCoarseToCoarse()) {
				for (Side<2> s : Side<2>::getValues()) {
					if (patch.piinfo->pinfo->hasNbr(s)) {
						if (patch.piinfo->pinfo->getNbrType(s) == NbrType::Normal) {
							auto iface_info = patch.piinfo->getNormalIfaceInfo(s);

							min_local_index = min(min_local_index, iface_info->col_local_index);

						} else if (patch.piinfo->pinfo->getNbrType(s) == NbrType::Fine) {
							auto iface_info = patch.piinfo->getFineIfaceInfo(s);

							min_local_index = min(min_local_index, iface_info->col_local_index);

							min_local_index
							= min(min_local_index, iface_info->fine_col_local_indexes[0]);

							min_local_index
							= min(min_local_index, iface_info->fine_col_local_indexes[1]);

						} else if (patch.piinfo->pinfo->getNbrType(s) == NbrType::Coarse) {
							auto iface_info = patch.piinfo->getCoarseIfaceInfo(s);

							min_local_index = min(min_local_index, iface_info->col_local_index);

							min_local_index
							= min(min_local_index, iface_info->coarse_col_local_index);
						}
					}
				}
			}
		}
	}
	CHECK(min_local_index == 0);
}
TEST_CASE("Schur::InterfaceDomain<2> col local indexes are contiguous 2d_refined_complicated",
          "[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	set<int> local_indexes;

	for (auto iface : interface_domain.getInterfaces()) {
		local_indexes.insert(iface->local_index);

		for (auto patch : iface->patches) {
			if (patch.type.isNormal() || patch.type.isFineToFine()
			    || patch.type.isCoarseToCoarse()) {
				for (Side<2> s : Side<2>::getValues()) {
					if (patch.piinfo->pinfo->hasNbr(s)) {
						if (patch.piinfo->pinfo->getNbrType(s) == NbrType::Normal) {
							auto iface_info = patch.piinfo->getNormalIfaceInfo(s);

							local_indexes.insert(iface_info->col_local_index);

						} else if (patch.piinfo->pinfo->getNbrType(s) == NbrType::Fine) {
							auto iface_info = patch.piinfo->getFineIfaceInfo(s);

							local_indexes.insert(iface_info->col_local_index);

							local_indexes.insert(iface_info->fine_col_local_indexes[0]);

							local_indexes.insert(iface_info->fine_col_local_indexes[1]);

						} else if (patch.piinfo->pinfo->getNbrType(s) == NbrType::Coarse) {
							auto iface_info = patch.piinfo->getCoarseIfaceInfo(s);

							local_indexes.insert(iface_info->col_local_index);

							local_indexes.insert(iface_info->coarse_col_local_index);
						}
					}
				}
			}
		}
	}
	int prev_local_index = -1;
	for (int local_index : local_indexes) {
		REQUIRE(local_index == prev_local_index + 1);
		prev_local_index = local_index;
	}
}
TEST_CASE(
"Schur::InterfaceDomain<2> each id has only one col local index associated with it 2d_refined_complicated",
"[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	map<int, set<int>> id_to_local_indexes;

	for (auto iface : interface_domain.getInterfaces()) {
		id_to_local_indexes[iface->id].insert(iface->local_index);

		for (auto patch : iface->patches) {
			if (patch.type.isNormal() || patch.type.isFineToFine()
			    || patch.type.isCoarseToCoarse()) {
				for (Side<2> s : Side<2>::getValues()) {
					if (patch.piinfo->pinfo->hasNbr(s)) {
						if (patch.piinfo->pinfo->getNbrType(s) == NbrType::Normal) {
							auto iface_info = patch.piinfo->getNormalIfaceInfo(s);

							id_to_local_indexes[iface_info->id].insert(iface_info->col_local_index);

						} else if (patch.piinfo->pinfo->getNbrType(s) == NbrType::Fine) {
							auto iface_info = patch.piinfo->getFineIfaceInfo(s);

							id_to_local_indexes[iface_info->id].insert(iface_info->col_local_index);

							id_to_local_indexes[iface_info->fine_ids[0]].insert(
							iface_info->fine_col_local_indexes[0]);

							id_to_local_indexes[iface_info->fine_ids[1]].insert(
							iface_info->fine_col_local_indexes[1]);

						} else if (patch.piinfo->pinfo->getNbrType(s) == NbrType::Coarse) {
							auto iface_info = patch.piinfo->getCoarseIfaceInfo(s);

							id_to_local_indexes[iface_info->id].insert(iface_info->col_local_index);

							id_to_local_indexes[iface_info->coarse_id].insert(
							iface_info->coarse_col_local_index);
						}
					}
				}
			}
		}
	}

	REQUIRE(id_to_local_indexes.size() > 0);
	for (auto pair : id_to_local_indexes) {
		INFO("ID " << pair.first);
		CHECK(pair.second.size() == 1);
	}
}
TEST_CASE(
"Schur::InterfaceDomain<2> each col local index has only one id associated with it 2d_refined_complicated",
"[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	map<int, set<int>> local_index_to_ids;

	for (auto iface : interface_domain.getInterfaces()) {
		local_index_to_ids[iface->local_index].insert(iface->id);

		for (auto patch : iface->patches) {
			if (patch.type.isNormal() || patch.type.isFineToFine()
			    || patch.type.isCoarseToCoarse()) {
				for (Side<2> s : Side<2>::getValues()) {
					if (patch.piinfo->pinfo->hasNbr(s)) {
						if (patch.piinfo->pinfo->getNbrType(s) == NbrType::Normal) {
							auto iface_info = patch.piinfo->getNormalIfaceInfo(s);
							local_index_to_ids[iface_info->col_local_index].insert(iface_info->id);
						} else if (patch.piinfo->pinfo->getNbrType(s) == NbrType::Fine) {
							auto iface_info = patch.piinfo->getFineIfaceInfo(s);
							local_index_to_ids[iface_info->col_local_index].insert(iface_info->id);
							local_index_to_ids[iface_info->fine_col_local_indexes[0]].insert(
							iface_info->fine_ids[0]);
							local_index_to_ids[iface_info->fine_col_local_indexes[1]].insert(
							iface_info->fine_ids[1]);
						} else if (patch.piinfo->pinfo->getNbrType(s) == NbrType::Coarse) {
							auto iface_info = patch.piinfo->getCoarseIfaceInfo(s);
							local_index_to_ids[iface_info->col_local_index].insert(iface_info->id);
							local_index_to_ids[iface_info->coarse_col_local_index].insert(
							iface_info->coarse_id);
						}
					}
				}
			}
		}
	}

	REQUIRE(local_index_to_ids.size() > 0);
	for (auto pair : local_index_to_ids) {
		INFO("Local Index " << pair.first);
		CHECK(pair.second.size() == 1);
	}
}
TEST_CASE("Schur::InterfaceDomain<2> global indexes start from 0 2d_refined_complicated",
          "[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	int min_global_index = numeric_limits<int>::max();

	for (auto iface : interface_domain.getInterfaces()) {
		min_global_index = min(min_global_index, iface->global_index);

		for (auto patch : iface->patches) {
			for (Side<2> s : Side<2>::getValues()) {
				if (patch.piinfo->pinfo->hasNbr(s)) {
					min_global_index
					= min(min_global_index, patch.piinfo->getIfaceInfo(s)->global_index);

					switch (patch.piinfo->pinfo->getNbrType(s)) {
						case NbrType::Coarse:
							min_global_index
							= min(min_global_index,
							      patch.piinfo->getCoarseIfaceInfo(s)->coarse_global_index);
							break;
						case NbrType::Fine:
							min_global_index
							= min(min_global_index,
							      patch.piinfo->getFineIfaceInfo(s)->fine_global_indexes[0]);
							min_global_index
							= min(min_global_index,
							      patch.piinfo->getFineIfaceInfo(s)->fine_global_indexes[1]);
							break;
						default:
							break;
					}
				}
			}
		}
	}

	for (auto patch : interface_domain.getPatchIfaceInfos()) {
		for (Side<2> s : Side<2>::getValues()) {
			if (patch->pinfo->hasNbr(s)) {
				min_global_index = min(min_global_index, patch->getIfaceInfo(s)->global_index);
			}
		}
	}

	int global_min_global_index;
	MPI_Allreduce(&min_global_index, &global_min_global_index, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

	CHECK(global_min_global_index == 0);
}
TEST_CASE("Schur::InterfaceDomain<2> global indexes are contiguous 2d_refined_complicated",
          "[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::set<int> global_indexes;

	for (auto iface : interface_domain.getInterfaces()) {
		global_indexes.insert(iface->global_index);
	}

	int num_ifaces = interface_domain.getNumLocalInterfaces();
	int prev_global_index;
	MPI_Scan(&num_ifaces, &prev_global_index, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	prev_global_index -= num_ifaces;
	prev_global_index -= 1;

	for (int global_index : global_indexes) {
		REQUIRE(global_index == prev_global_index + 1);
		prev_global_index = global_index;
	}
}
TEST_CASE(
"Schur::InterfaceDomain<2> each id has only one global index associated with it 2d_refined_complicated",
"[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	map<int, set<int>> id_to_global_indexes;

	vector<int> iface_ids_in(interface_domain.getNumGlobalInterfaces());
	for (auto iface : interface_domain.getInterfaces()) {
		if (iface->global_index >= 0 && iface->global_index < (int) iface_ids_in.size()) {
			iface_ids_in[iface->global_index] = iface->id;
		}
	}

	vector<int> iface_ids_out(interface_domain.getNumGlobalInterfaces());
	MPI_Allreduce(iface_ids_in.data(), iface_ids_out.data(), iface_ids_in.size(), MPI_INT, MPI_SUM,
	              MPI_COMM_WORLD);

	for (int global_index = 0; global_index < interface_domain.getNumGlobalInterfaces();
	     global_index++) {
		id_to_global_indexes[iface_ids_out[global_index]].insert(global_index);
	}

	for (auto iface : interface_domain.getInterfaces()) {
		id_to_global_indexes[iface->id].insert(iface->global_index);

		for (auto patch : iface->patches) {
			for (Side<2> s : Side<2>::getValues()) {
				if (patch.piinfo->pinfo->hasNbr(s)) {
					id_to_global_indexes[patch.piinfo->getIfaceInfo(s)->id].insert(
					patch.piinfo->getIfaceInfo(s)->global_index);

					switch (patch.piinfo->pinfo->getNbrType(s)) {
						case NbrType::Coarse:
							id_to_global_indexes[patch.piinfo->getCoarseIfaceInfo(s)->coarse_id]
							.insert(patch.piinfo->getCoarseIfaceInfo(s)->coarse_global_index);
							break;
						case NbrType::Fine:
							id_to_global_indexes[patch.piinfo->getFineIfaceInfo(s)->fine_ids[0]]
							.insert(patch.piinfo->getFineIfaceInfo(s)->fine_global_indexes[0]);
							id_to_global_indexes[patch.piinfo->getFineIfaceInfo(s)->fine_ids[1]]
							.insert(patch.piinfo->getFineIfaceInfo(s)->fine_global_indexes[1]);
							break;
						default:
							break;
					}
				}
			}
		}
	}

	for (auto piinfo : interface_domain.getPatchIfaceInfos()) {
		for (Side<2> s : Side<2>::getValues()) {
			if (piinfo->pinfo->hasNbr(s)) {
				id_to_global_indexes[piinfo->getIfaceInfo(s)->id].insert(
				piinfo->getIfaceInfo(s)->global_index);
			}
		}
	}

	REQUIRE(id_to_global_indexes.size() > 0);
	for (auto pair : id_to_global_indexes) {
		INFO("ID " << pair.first);
		string global_indexes;
		for (int global_index : pair.second) {
			global_indexes += to_string(global_index) + ", ";
		}
		INFO("GLOBAL_INDEXES: " << global_indexes);
		CHECK(pair.second.size() == 1);
	}
}
TEST_CASE(
"Schur::InterfaceDomain<2> each global index has only one id associated with it 2d_refined_complicated",
"[Schur::InterfaceDomain]")
{
	DomainReader<2> domain_reader("mesh_inputs/2d_refined_complicated_mpi3.json", {10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> interface_domain(domain);

	map<int, set<int>> global_index_to_ids;

	vector<int> iface_ids_in(interface_domain.getNumGlobalInterfaces());
	for (auto iface : interface_domain.getInterfaces()) {
		if (iface->global_index >= 0 && iface->global_index < (int) iface_ids_in.size()) {
			iface_ids_in[iface->global_index] = iface->id;
		}
	}

	vector<int> iface_ids_out(interface_domain.getNumGlobalInterfaces());
	MPI_Allreduce(iface_ids_in.data(), iface_ids_out.data(), iface_ids_in.size(), MPI_INT, MPI_MAX,
	              MPI_COMM_WORLD);

	for (int global_index = 0; global_index < interface_domain.getNumGlobalInterfaces();
	     global_index++) {
		global_index_to_ids[global_index].insert(iface_ids_out[global_index]);
	}
	for (auto iface : interface_domain.getInterfaces()) {
		global_index_to_ids[iface->global_index].insert(iface->id);

		for (auto patch : iface->patches) {
			for (Side<2> s : Side<2>::getValues()) {
				if (patch.piinfo->pinfo->hasNbr(s)) {
					global_index_to_ids[patch.piinfo->getIfaceInfo(s)->global_index].insert(
					patch.piinfo->getIfaceInfo(s)->id);

					switch (patch.piinfo->pinfo->getNbrType(s)) {
						case NbrType::Coarse:
							global_index_to_ids[patch.piinfo->getCoarseIfaceInfo(s)
							                    ->coarse_global_index]
							.insert(patch.piinfo->getCoarseIfaceInfo(s)->coarse_id);
							break;
						case NbrType::Fine:
							global_index_to_ids[patch.piinfo->getFineIfaceInfo(s)
							                    ->fine_global_indexes[0]]
							.insert(patch.piinfo->getFineIfaceInfo(s)->fine_ids[0]);
							global_index_to_ids[patch.piinfo->getFineIfaceInfo(s)
							                    ->fine_global_indexes[1]]
							.insert(patch.piinfo->getFineIfaceInfo(s)->fine_ids[1]);
							break;
						default:
							break;
					}
				}
			}
		}
	}

	for (auto piinfo : interface_domain.getPatchIfaceInfos()) {
		for (Side<2> s : Side<2>::getValues()) {
			if (piinfo->pinfo->hasNbr(s)) {
				global_index_to_ids[piinfo->getIfaceInfo(s)->global_index].insert(
				piinfo->getIfaceInfo(s)->id);
			}
		}
	}

	REQUIRE(global_index_to_ids.size() > 0);
	for (auto pair : global_index_to_ids) {
		INFO("Global Index " << pair.first);
		CHECK(pair.second.size() == 1);
	}
}