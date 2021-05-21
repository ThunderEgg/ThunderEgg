/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Gebneral Public Licenbse as published by
 *  the Free Software Foundation, either version 3 of the Licenbse, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be ubseful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Gebneral Public Licenbse for more details.
 *
 *  You should have received a copy of the GNU Gebneral Public Licenbse
 *  along with this program.  If not, bsee <https://www.gnu.org/licenbses/>.
 ***************************************************************************/

#include "P8estDomainGenerator_SHARED.h"
#include "mpi.h"
#include <algorithm>
#include <catch2/catch_approx.hpp>
using namespace std;
using namespace ThunderEgg;

#include <catch2/catch_test_macros.hpp>

PatchVector::PatchVector(const ThunderEgg::Domain<3> &domain, int max_level)
: max_level(max_level), pinfos(GetAllPatchesOnRank0(domain))
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		int n = 0b1 << max_level;
		array.resize(n);
		for (auto &sub : array) {
			sub.resize(n);
			for (auto &subsub : sub) {
				subsub.resize(n);
			}
		}
		for (const PatchInfo<3> &pinfo : pinfos) {
			int i          = pinfo.starts[0] * n;
			int j          = pinfo.starts[1] * n;
			int k          = pinfo.starts[2] * n;
			array[i][j][k] = &pinfo;
		}
	}
}
const ThunderEgg::PatchInfo<3> *PatchVector::operator[](const std::string &str) const
{
	int i = 0;
	int j = 0;
	int k = 0;

	int  curr_level = 0;
	auto iter       = str.cbegin();
	while (iter != str.cend()) {
		switch (*iter) {
			case 'w':
				i <<= 1;
				break;
			case 'e':
				i <<= 1;
				i |= 1;
				break;
			case 's':
				j <<= 1;
				break;
			case 'n':
				j <<= 1;
				j |= 1;
				break;
			case 'b':
				k <<= 1;
				break;
			case 't':
				k <<= 1;
				k |= 1;
				break;
			case '_':
				curr_level++;
				break;
			default:
				throw 0;
		}
		++iter;
	}
	curr_level++;
	while (curr_level < max_level) {
		i <<= 1;
		j <<= 1;
		k <<= 1;
		curr_level++;
	}
	const PatchInfo<3> *pinfo = array[i][j][k];
	REQUIRE(pinfo != nullptr);
	return pinfo;
}
std::vector<PatchInfo<3>> GetAllPatchesOnRank0(const Domain<3> &domain)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	std::vector<PatchInfo<3>> all_patches;

	if (rank != 0) {
		nlohmann::json patches        = domain.getPatchInfoVector();
		string         patches_string = patches.dump();
		MPI_Send(patches_string.data(), (int) patches_string.size() + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	} else {
		for (int i = 0; i < size - 1; i++) {
			MPI_Status status;
			MPI_Probe(1, 0, MPI_COMM_WORLD, &status);

			int buffer_size;
			MPI_Get_count(&status, MPI_CHAR, &buffer_size);

			char patches_string[buffer_size];
			MPI_Recv(patches_string, buffer_size, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status);

			nlohmann::json patches = nlohmann::json::parse(patches_string);
			if (patches != nullptr) {
				patches.get_to(all_patches);
				for (auto &patch : all_patches) {
					patch.ns = domain.getNs();
					for (int i = 0; i < 3; i++) {
						patch.spacings[i] /= patch.ns[i];
					}
				}
			}
		}
		for (auto patch : domain.getPatchInfoVector()) {
			all_patches.push_back(patch);
		}
	}
	return all_patches;
}

void Ident(int block_no, double unit_x, double unit_y, double unit_z, double &x, double &y, double &z)
{
	x = unit_x;
	y = unit_y;
	z = unit_z;
}

void CheckRootDomainNeighbors(const ThunderEgg::Domain<3> &domain)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<PatchInfo<3>> patches = GetAllPatchesOnRank0(domain);
	if (rank == 0) {
		const PatchInfo<3> &root_patch = patches[0];
		for (Side<3> s : Side<3>::getValues()) {
			CHECK_FALSE(root_patch.hasNbr(s));
		}
		for (Edge e : Edge::getValues()) {
			CHECK_FALSE(root_patch.hasNbr(e));
		}
		for (Corner<3> c : Corner<3>::getValues()) {
			CHECK_FALSE(root_patch.hasNbr(c));
		}
	}
}
namespace
{
std::string getString(int max_level, int i, int j, int k)
{
	std::string str;
	for (int n = 0; n < max_level; n++) {
		if (i & 0b1) {
			str = 'e' + str;
		} else {
			str = 'w' + str;
		}
		i >>= 1;
		if (j & 0b1) {
			str = 'n' + str;
		} else {
			str = 's' + str;
		}
		j >>= 1;
		if (k & 0b1) {
			str = 't' + str;
		} else {
			str = 'b' + str;
		}
		k >>= 1;
		str = '_' + str;
	}
	return str.substr(1, str.size());
}
Orthant<3> getChildOrthant(const string &str)
{
	int    val = 0;
	string sub = str.substr(str.size() - 3, str.size());
	if (sub[0] == 't') {
		val |= 0b1 << 2;
	}
	if (sub[1] == 'n') {
		val |= 0b1 << 1;
	}
	if (sub[2] == 'e') {
		val |= 0b1 << 0;
	}
	return Orthant<3>(val);
}
} // namespace
void CheckParentAndChildIdsAndRanks(const ThunderEgg::Domain<3> &coarser_domain, int coarser_max_level, const ThunderEgg::Domain<3> &finer_domain, int finer_max_level)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	PatchVector coarser_pvector(coarser_domain, coarser_max_level);
	PatchVector finer_pvector(finer_domain, finer_max_level);
	//
	if (rank == 0) {
		int n = 0b1 << finer_max_level;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					string child_str  = getString(finer_max_level, i, j, k);
					string parent_str = child_str.substr(0, std::max<int>((int) child_str.size() - 4, 0));

					const PatchInfo<3> *child_patch  = finer_pvector[child_str];
					const PatchInfo<3> *parent_patch = coarser_pvector[parent_str];

					Orthant<3> orth = getChildOrthant(child_str);
					CHECK(child_patch->parent_id == parent_patch->id);
					CHECK(parent_patch->child_ids[orth.getIndex()] == child_patch->id);

					CHECK(child_patch->parent_rank == parent_patch->rank);
					CHECK(parent_patch->child_ranks[orth.getIndex()] == child_patch->rank);

					CHECK(child_patch->orth_on_parent == orth);
				}
			}
		}
	}
}
void CheckParentAndChildIdsAndRanksRefined(const ThunderEgg::Domain<3> &coarser_domain, int coarser_max_level, const ThunderEgg::Domain<3> &finer_domain, int finer_max_level)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	PatchVector coarser_pvector(coarser_domain, coarser_max_level);
	PatchVector finer_pvector(finer_domain, finer_max_level);
	//
	if (rank == 0) {
		int n = 0b1 << finer_max_level;
		INFO("N: " << n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					if (i < 2 && j < 2 && k < 2) {
						string child_str  = getString(finer_max_level, i, j, k);
						string parent_str = child_str.substr(0, std::max<int>((int) child_str.size() - 4, 0));

						INFO("Child Patch: " << child_str);
						INFO("Praent Patch: " << parent_str);
						const PatchInfo<3> *child_patch  = finer_pvector[child_str];
						const PatchInfo<3> *parent_patch = coarser_pvector[parent_str];

						Orthant<3> orth = getChildOrthant(child_str);
						CHECK(child_patch->parent_id == parent_patch->id);
						CHECK(parent_patch->child_ids[orth.getIndex()] == child_patch->id);

						CHECK(child_patch->parent_rank == parent_patch->rank);
						CHECK(parent_patch->child_ranks[orth.getIndex()] == child_patch->rank);

						CHECK(child_patch->orth_on_parent == orth);
					} else if (i % 2 == 0 && j % 2 == 0 && k % 2 == 0) {
						string child_str  = getString(finer_max_level, i, j, k);
						string parent_str = child_str.substr(0, std::max<int>((int) child_str.size() - 4, 0));

						INFO("Patch: " << parent_str);
						const PatchInfo<3> *child_patch  = finer_pvector[parent_str];
						const PatchInfo<3> *parent_patch = coarser_pvector[parent_str];

						CHECK(child_patch->parent_id == parent_patch->id);
						CHECK(parent_patch->child_ids[0] == child_patch->id);
						for (int i = 1; i < 8; i++) {
							CHECK(parent_patch->child_ids[i] == -1);
						}

						CHECK(child_patch->parent_rank == parent_patch->rank);
						CHECK(parent_patch->child_ranks[0] == child_patch->rank);
						for (int i = 1; i < 8; i++) {
							CHECK(parent_patch->child_ranks[i] == -1);
						}

						CHECK(child_patch->orth_on_parent == Orthant<3>::null());
					}
				}
			}
		}
	}
}
void CheckParentIdsAndRanksNull(const ThunderEgg::Domain<3> &domain)
{
	for (const PatchInfo<3> &pinfo : domain.getPatchInfoVector()) {
		CHECK(pinfo.parent_id == -1);
		CHECK(pinfo.parent_rank == -1);
		CHECK(pinfo.orth_on_parent == Orthant<3>::null());
	}
}
void CheckChildIdsAndRanksNull(const ThunderEgg::Domain<3> &domain)
{
	for (const PatchInfo<3> &pinfo : domain.getPatchInfoVector()) {
		for (int i = 0; i < 8; i++) {
			CHECK(pinfo.child_ids[i] == -1);
			CHECK(pinfo.child_ranks[i] == -1);
		}
	}
}

void Check4x4x4DomainSideHasNeighbors(const PatchVector &domain);
void Check4x4x4DomainSideNeighborIds(const PatchVector &domain);
void Check4x4x4DomainSideNeighborRanks(const PatchVector &domain);
void Check4x4x4DomainEdgeHasNeighbors(const PatchVector &domain);
void Check4x4x4DomainEdgeNeighborIds(const PatchVector &domain);
void Check4x4x4DomainEdgeNeighborRanks(const PatchVector &domain);
void Check4x4x4DomainCornerHasNeighbors(const PatchVector &domain);
void Check4x4x4DomainCornerNeighborIds(const PatchVector &domain);
void Check4x4x4DomainCornerNeighborRanks(const PatchVector &domain);

void Check4x4x4DomainNeighbors(const ThunderEgg::Domain<3> &domain)
{
	PatchVector pv(domain, 2);
	int         rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		Check4x4x4DomainSideHasNeighbors(pv);
		Check4x4x4DomainSideNeighborIds(pv);
		Check4x4x4DomainSideNeighborRanks(pv);
		Check4x4x4DomainEdgeHasNeighbors(pv);
		Check4x4x4DomainEdgeNeighborIds(pv);
		Check4x4x4DomainEdgeNeighborRanks(pv);
		Check4x4x4DomainCornerHasNeighbors(pv);
		Check4x4x4DomainCornerNeighborIds(pv);
		Check4x4x4DomainCornerNeighborRanks(pv);
	}
}

void Check4x4x4RefinedBSWDomainSideHasNeighbors(const PatchVector &domain);
void Check4x4x4RefinedBSWDomainSideNeighborIds(const PatchVector &domain);
void Check4x4x4RefinedBSWDomainSideNeighborRanks(const PatchVector &domain);
void Check4x4x4RefinedBSWDomainSideNeighborOrths(const PatchVector &domain);
void Check4x4x4RefinedBSWDomainEdgeHasNeighbors(const PatchVector &domain);
void Check4x4x4RefinedBSWDomainEdgeNeighborIds(const PatchVector &domain);
void Check4x4x4RefinedBSWDomainEdgeNeighborRanks(const PatchVector &domain);
void Check4x4x4RefinedBSWDomainEdgeNeighborOrths(const PatchVector &domain);
void Check4x4x4RefinedBSWDomainCornerHasNeighbors(const PatchVector &domain);
void Check4x4x4RefinedBSWDomainCornerNeighborIds(const PatchVector &domain);
void Check4x4x4RefinedBSWDomainCornerNeighborRanks(const PatchVector &domain);
void Check4x4x4RefinedBSWDomainCornerNeighborOrths(const PatchVector &domain);

void Check4x4x4RefinedBSWDomainNeighbors(const ThunderEgg::Domain<3> &domain)
{
	PatchVector pv(domain, 3);
	int         rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		Check4x4x4RefinedBSWDomainSideHasNeighbors(pv);
		Check4x4x4RefinedBSWDomainSideNeighborIds(pv);
		Check4x4x4RefinedBSWDomainSideNeighborRanks(pv);
		Check4x4x4RefinedBSWDomainSideNeighborOrths(pv);
		Check4x4x4RefinedBSWDomainEdgeHasNeighbors(pv);
		Check4x4x4RefinedBSWDomainEdgeNeighborIds(pv);
		Check4x4x4RefinedBSWDomainEdgeNeighborRanks(pv);
		Check4x4x4RefinedBSWDomainEdgeNeighborOrths(pv);
		Check4x4x4RefinedBSWDomainCornerHasNeighbors(pv);
		Check4x4x4RefinedBSWDomainCornerNeighborIds(pv);
		Check4x4x4RefinedBSWDomainCornerNeighborRanks(pv);
		Check4x4x4RefinedBSWDomainCornerNeighborOrths(pv);
	}
}