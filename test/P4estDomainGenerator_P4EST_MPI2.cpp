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

#include <ThunderEgg/P4estDomainGenerator.h>

#include <p4est.h>
#include <p4est_mesh.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

namespace
{
std::vector<PatchInfo<2>> GetAllPatchesOnRank0(std::shared_ptr<const Domain<2>> domain)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	std::vector<PatchInfo<2>> all_patches;

	if (rank == 1) {
		nlohmann::json patches;
		for (auto patch : domain->getPatchInfoVector()) {
			patches.push_back(patch);
		}
		string patches_string = patches.dump();
		MPI_Send(patches_string.data(), (int) patches_string.size() + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	} else {
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
				patch.ns = domain->getNs();
				for (int i = 0; i < 2; i++) {
					patch.spacings[i] /= patch.ns[i];
				}
			}
		}
		for (auto patch : domain->getPatchInfoVector()) {
			all_patches.push_back(patch);
		}
	}
	return all_patches;
}
std::string GetAllPatchesJSONString(const std::vector<PatchInfo<2>> &patches)
{
	nlohmann::json patches_j = patches;
	return patches_j.dump(1);
}
} // namespace
TEST_CASE("P4estDomainGenerator 4x4 Uniform", "[p4estDomGen]")
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	p4est_connectivity_t *conn = p4est_connectivity_new_unitsquare();

	p4est_t *p4est = p4est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

	p4est_refine(
	p4est, false,
	[](p4est_t *p4est, p4est_topidx_t witch_tree, p4est_quadrant_t *quadrant) -> int { return 1; },
	nullptr);
	p4est_refine(
	p4est, false,
	[](p4est_t *p4est, p4est_topidx_t witch_tree, p4est_quadrant_t *quadrant) -> int { return 1; },
	nullptr);

	p4est_partition(p4est, true, nullptr);

	int    nx              = GENERATE(5, 10);
	int    ny              = GENERATE(5, 10);
	double scale_x         = GENERATE(0.5, 1.0);
	double scale_y         = GENERATE(0.5, 1.0);
	int    num_ghost_cells = GENERATE(0, 1, 2);

	INFO("nx: " << nx);
	INFO("ny: " << ny);
	INFO("scale_x: " << scale_x);
	INFO("scale_y: " << scale_y);
	INFO("num_ghost_cells: " << num_ghost_cells);

	P4estDomainGenerator::BlockMapFunc bmf
	= [&](int block_no, double unit_x, double unit_y, double &x, double &y) {
		  x = scale_x * unit_x;
		  y = scale_y * unit_y;
	  };

	P4estDomainGenerator dg(p4est, {nx, ny}, num_ghost_cells, bmf);

	auto domain_2 = dg.getFinestDomain();
	auto domain_1 = dg.getCoarserDomain();
	auto domain_0 = dg.getCoarserDomain();

	//SECTION("correct number of patches")
	{
		CHECK(domain_2->getNumGlobalPatches() == 16);
		CHECK(domain_1->getNumGlobalPatches() == 4);
		CHECK(domain_0->getNumGlobalPatches() == 1);
	}
	//SECTION("patches have correct spacings")
	{
		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
			CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.spacings[0] == Catch::Approx(scale_x * 1.0 / nx));
			CHECK(patch.spacings[1] == Catch::Approx(scale_y * 1.0 / ny));
		}
	}
	//SECTION("patches have correct ns")
	{
		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
		}
	}

	//SECTION("patches have refine_level set")
	{
		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.refine_level == 2);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.refine_level == 1);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.refine_level == 0);
		}
	}
	//SECTION("patches have ranks set")
	{
		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.rank == rank);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.rank == rank);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.rank == rank);
		}
	}
	//SECTION("patches have num_ghost_cells set")
	{
		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.num_ghost_cells == num_ghost_cells);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.num_ghost_cells == num_ghost_cells);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.num_ghost_cells == num_ghost_cells);
		}
	}

	std::vector<PatchInfo<2>> domain_2_patches = GetAllPatchesOnRank0(domain_2);
	std::vector<PatchInfo<2>> domain_1_patches = GetAllPatchesOnRank0(domain_1);
	std::vector<PatchInfo<2>> domain_0_patches = GetAllPatchesOnRank0(domain_0);

	INFO("nx: " << nx);
	INFO("ny: " << ny);
	INFO("scale_x: " << scale_x);
	INFO("scale_y: " << scale_y);
	INFO("num_ghost_cells: " << num_ghost_cells);

	const PatchInfo<2> *domain_2_sw_sw_patch = nullptr;
	const PatchInfo<2> *domain_2_sw_se_patch = nullptr;
	const PatchInfo<2> *domain_2_sw_nw_patch = nullptr;
	const PatchInfo<2> *domain_2_sw_ne_patch = nullptr;

	const PatchInfo<2> *domain_2_se_sw_patch = nullptr;
	const PatchInfo<2> *domain_2_se_se_patch = nullptr;
	const PatchInfo<2> *domain_2_se_nw_patch = nullptr;
	const PatchInfo<2> *domain_2_se_ne_patch = nullptr;

	const PatchInfo<2> *domain_2_nw_sw_patch = nullptr;
	const PatchInfo<2> *domain_2_nw_se_patch = nullptr;
	const PatchInfo<2> *domain_2_nw_nw_patch = nullptr;
	const PatchInfo<2> *domain_2_nw_ne_patch = nullptr;

	const PatchInfo<2> *domain_2_ne_sw_patch = nullptr;
	const PatchInfo<2> *domain_2_ne_se_patch = nullptr;
	const PatchInfo<2> *domain_2_ne_nw_patch = nullptr;
	const PatchInfo<2> *domain_2_ne_ne_patch = nullptr;

	const PatchInfo<2> *domain_1_sw_patch = nullptr;
	const PatchInfo<2> *domain_1_se_patch = nullptr;
	const PatchInfo<2> *domain_1_nw_patch = nullptr;
	const PatchInfo<2> *domain_1_ne_patch = nullptr;

	const PatchInfo<2> *domain_0_coarser_patch = nullptr;

	if (rank == 0) {
		for (PatchInfo<2> &patch : domain_2_patches) {
			double x = patch.starts[0];
			double y = patch.starts[1];
			if (x == Catch::Approx(0) && y == Catch::Approx(0)) {
				domain_2_sw_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0)) {
				domain_2_sw_se_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0)) {
				domain_2_se_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.75 * scale_x) && y == Catch::Approx(0)) {
				domain_2_se_se_patch = &patch;
			}

			if (x == Catch::Approx(0) && y == Catch::Approx(0.25 * scale_y)) {
				domain_2_sw_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.25 * scale_y)) {
				domain_2_sw_ne_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.25 * scale_y)) {
				domain_2_se_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.75 * scale_x) && y == Catch::Approx(0.25 * scale_y)) {
				domain_2_se_ne_patch = &patch;
			}

			if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y)) {
				domain_2_nw_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_2_nw_se_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_2_ne_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.75 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_2_ne_se_patch = &patch;
			}

			if (x == Catch::Approx(0) && y == Catch::Approx(0.75 * scale_y)) {
				domain_2_nw_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.75 * scale_y)) {
				domain_2_nw_ne_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.75 * scale_y)) {
				domain_2_ne_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.75 * scale_x) && y == Catch::Approx(0.75 * scale_y)) {
				domain_2_ne_ne_patch = &patch;
			}
		}

		for (PatchInfo<2> &patch : domain_1_patches) {
			double x = patch.starts[0];
			double y = patch.starts[1];
			if (x == Catch::Approx(0) && y == Catch::Approx(0)) {
				domain_1_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0)) {
				domain_1_se_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y)) {
				domain_1_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_1_ne_patch = &patch;
			}
		}

		domain_0_coarser_patch = &domain_0_patches[0];

		REQUIRE(domain_2_sw_sw_patch != nullptr);
		REQUIRE(domain_2_sw_se_patch != nullptr);
		REQUIRE(domain_2_sw_nw_patch != nullptr);
		REQUIRE(domain_2_sw_ne_patch != nullptr);

		REQUIRE(domain_2_se_sw_patch != nullptr);
		REQUIRE(domain_2_se_se_patch != nullptr);
		REQUIRE(domain_2_se_nw_patch != nullptr);
		REQUIRE(domain_2_se_ne_patch != nullptr);

		REQUIRE(domain_2_nw_sw_patch != nullptr);
		REQUIRE(domain_2_nw_se_patch != nullptr);
		REQUIRE(domain_2_nw_nw_patch != nullptr);
		REQUIRE(domain_2_nw_ne_patch != nullptr);

		REQUIRE(domain_2_ne_sw_patch != nullptr);
		REQUIRE(domain_2_ne_se_patch != nullptr);
		REQUIRE(domain_2_ne_nw_patch != nullptr);
		REQUIRE(domain_2_ne_ne_patch != nullptr);

		REQUIRE(domain_1_sw_patch != nullptr);
		REQUIRE(domain_1_se_patch != nullptr);
		REQUIRE(domain_1_nw_patch != nullptr);
		REQUIRE(domain_1_ne_patch != nullptr);

		CHECK(domain_0_coarser_patch->starts[0] == Catch::Approx(0.0));
		CHECK(domain_0_coarser_patch->starts[1] == Catch::Approx(0.0));
	}

	//SECTION("parent ids are set correctly")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_sw_patch->parent_id == domain_1_sw_patch->id);
			CHECK(domain_2_sw_se_patch->parent_id == domain_1_sw_patch->id);
			CHECK(domain_2_sw_nw_patch->parent_id == domain_1_sw_patch->id);
			CHECK(domain_2_sw_ne_patch->parent_id == domain_1_sw_patch->id);

			CHECK(domain_2_se_sw_patch->parent_id == domain_1_se_patch->id);
			CHECK(domain_2_se_se_patch->parent_id == domain_1_se_patch->id);
			CHECK(domain_2_se_nw_patch->parent_id == domain_1_se_patch->id);
			CHECK(domain_2_se_ne_patch->parent_id == domain_1_se_patch->id);

			CHECK(domain_2_nw_sw_patch->parent_id == domain_1_nw_patch->id);
			CHECK(domain_2_nw_se_patch->parent_id == domain_1_nw_patch->id);
			CHECK(domain_2_nw_nw_patch->parent_id == domain_1_nw_patch->id);
			CHECK(domain_2_nw_ne_patch->parent_id == domain_1_nw_patch->id);

			CHECK(domain_2_ne_sw_patch->parent_id == domain_1_ne_patch->id);
			CHECK(domain_2_ne_se_patch->parent_id == domain_1_ne_patch->id);
			CHECK(domain_2_ne_nw_patch->parent_id == domain_1_ne_patch->id);
			CHECK(domain_2_ne_ne_patch->parent_id == domain_1_ne_patch->id);

			CHECK(domain_1_sw_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_se_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_nw_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_ne_patch->parent_id == domain_0_coarser_patch->id);

			CHECK(domain_0_coarser_patch->parent_id == -1);
		}
	}
	//SECTION("child ids are set correctly")
	{
		if (rank == 0) {
			for (int i = 0; i < 4; i++) {
				INFO("i: " << i);
				CHECK(domain_2_sw_sw_patch->child_ids[i] == -1);
				CHECK(domain_2_sw_se_patch->child_ids[i] == -1);
				CHECK(domain_2_sw_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_sw_ne_patch->child_ids[i] == -1);

				CHECK(domain_2_se_sw_patch->child_ids[i] == -1);
				CHECK(domain_2_se_se_patch->child_ids[i] == -1);
				CHECK(domain_2_se_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_se_ne_patch->child_ids[i] == -1);

				CHECK(domain_2_nw_sw_patch->child_ids[i] == -1);
				CHECK(domain_2_nw_se_patch->child_ids[i] == -1);
				CHECK(domain_2_nw_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_nw_ne_patch->child_ids[i] == -1);

				CHECK(domain_2_ne_sw_patch->child_ids[i] == -1);
				CHECK(domain_2_ne_se_patch->child_ids[i] == -1);
				CHECK(domain_2_ne_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_ne_ne_patch->child_ids[i] == -1);
			}

			CHECK(domain_1_sw_patch->child_ids[0] == domain_2_sw_sw_patch->id);
			CHECK(domain_1_sw_patch->child_ids[1] == domain_2_sw_se_patch->id);
			CHECK(domain_1_sw_patch->child_ids[2] == domain_2_sw_nw_patch->id);
			CHECK(domain_1_sw_patch->child_ids[3] == domain_2_sw_ne_patch->id);

			CHECK(domain_1_se_patch->child_ids[0] == domain_2_se_sw_patch->id);
			CHECK(domain_1_se_patch->child_ids[1] == domain_2_se_se_patch->id);
			CHECK(domain_1_se_patch->child_ids[2] == domain_2_se_nw_patch->id);
			CHECK(domain_1_se_patch->child_ids[3] == domain_2_se_ne_patch->id);

			CHECK(domain_1_nw_patch->child_ids[0] == domain_2_nw_sw_patch->id);
			CHECK(domain_1_nw_patch->child_ids[1] == domain_2_nw_se_patch->id);
			CHECK(domain_1_nw_patch->child_ids[2] == domain_2_nw_nw_patch->id);
			CHECK(domain_1_nw_patch->child_ids[3] == domain_2_nw_ne_patch->id);

			CHECK(domain_1_ne_patch->child_ids[0] == domain_2_ne_sw_patch->id);
			CHECK(domain_1_ne_patch->child_ids[1] == domain_2_ne_se_patch->id);
			CHECK(domain_1_ne_patch->child_ids[2] == domain_2_ne_nw_patch->id);
			CHECK(domain_1_ne_patch->child_ids[3] == domain_2_ne_ne_patch->id);

			CHECK(domain_0_coarser_patch->child_ids[0] == domain_1_sw_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[1] == domain_1_se_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[2] == domain_1_nw_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[3] == domain_1_ne_patch->id);
		}
	}
	//SECTION("orth on parent is set correctly")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_2_sw_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_2_sw_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_2_sw_ne_patch->orth_on_parent == Orthant<2>::ne());

			CHECK(domain_2_se_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_2_se_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_2_se_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_2_se_ne_patch->orth_on_parent == Orthant<2>::ne());

			CHECK(domain_2_nw_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_2_nw_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_2_nw_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_2_nw_ne_patch->orth_on_parent == Orthant<2>::ne());

			CHECK(domain_2_ne_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_2_ne_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_2_ne_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_2_ne_ne_patch->orth_on_parent == Orthant<2>::ne());

			CHECK(domain_1_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_1_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_1_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_1_ne_patch->orth_on_parent == Orthant<2>::ne());

			CHECK(domain_0_coarser_patch->orth_on_parent == Orthant<2>::null());
		}
	}
	//SECTION("parent ranks are set correctly")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_sw_patch->parent_rank == domain_1_sw_patch->rank);
			CHECK(domain_2_sw_se_patch->parent_rank == domain_1_sw_patch->rank);
			CHECK(domain_2_sw_nw_patch->parent_rank == domain_1_sw_patch->rank);
			CHECK(domain_2_sw_ne_patch->parent_rank == domain_1_sw_patch->rank);

			CHECK(domain_2_se_sw_patch->parent_rank == domain_1_se_patch->rank);
			CHECK(domain_2_se_se_patch->parent_rank == domain_1_se_patch->rank);
			CHECK(domain_2_se_nw_patch->parent_rank == domain_1_se_patch->rank);
			CHECK(domain_2_se_ne_patch->parent_rank == domain_1_se_patch->rank);

			CHECK(domain_2_nw_sw_patch->parent_rank == domain_1_nw_patch->rank);
			CHECK(domain_2_nw_se_patch->parent_rank == domain_1_nw_patch->rank);
			CHECK(domain_2_nw_nw_patch->parent_rank == domain_1_nw_patch->rank);
			CHECK(domain_2_nw_ne_patch->parent_rank == domain_1_nw_patch->rank);

			CHECK(domain_2_ne_sw_patch->parent_rank == domain_1_ne_patch->rank);
			CHECK(domain_2_ne_se_patch->parent_rank == domain_1_ne_patch->rank);
			CHECK(domain_2_ne_nw_patch->parent_rank == domain_1_ne_patch->rank);
			CHECK(domain_2_ne_ne_patch->parent_rank == domain_1_ne_patch->rank);

			CHECK(domain_1_sw_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_se_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_nw_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_ne_patch->parent_rank == domain_0_coarser_patch->rank);

			CHECK(domain_0_coarser_patch->parent_rank == -1);
		}
	}
	//SECTION("child ranks are set correctly")
	{
		if (rank == 0) {
			for (int i = 0; i < 4; i++) {
				INFO("i: " << i);
				CHECK(domain_2_sw_sw_patch->child_ids[i] == -1);
				CHECK(domain_2_sw_se_patch->child_ids[i] == -1);
				CHECK(domain_2_sw_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_sw_ne_patch->child_ids[i] == -1);

				CHECK(domain_2_se_sw_patch->child_ids[i] == -1);
				CHECK(domain_2_se_se_patch->child_ids[i] == -1);
				CHECK(domain_2_se_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_se_ne_patch->child_ids[i] == -1);

				CHECK(domain_2_nw_sw_patch->child_ids[i] == -1);
				CHECK(domain_2_nw_se_patch->child_ids[i] == -1);
				CHECK(domain_2_nw_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_nw_ne_patch->child_ids[i] == -1);

				CHECK(domain_2_ne_sw_patch->child_ids[i] == -1);
				CHECK(domain_2_ne_se_patch->child_ids[i] == -1);
				CHECK(domain_2_ne_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_ne_ne_patch->child_ids[i] == -1);
			}

			CHECK(domain_1_sw_patch->child_ranks[0] == domain_2_sw_sw_patch->rank);
			CHECK(domain_1_sw_patch->child_ranks[1] == domain_2_sw_se_patch->rank);
			CHECK(domain_1_sw_patch->child_ranks[2] == domain_2_sw_nw_patch->rank);
			CHECK(domain_1_sw_patch->child_ranks[3] == domain_2_sw_ne_patch->rank);

			CHECK(domain_1_se_patch->child_ranks[0] == domain_2_se_sw_patch->rank);
			CHECK(domain_1_se_patch->child_ranks[1] == domain_2_se_se_patch->rank);
			CHECK(domain_1_se_patch->child_ranks[2] == domain_2_se_nw_patch->rank);
			CHECK(domain_1_se_patch->child_ranks[3] == domain_2_se_ne_patch->rank);

			CHECK(domain_1_nw_patch->child_ranks[0] == domain_2_nw_sw_patch->rank);
			CHECK(domain_1_nw_patch->child_ranks[1] == domain_2_nw_se_patch->rank);
			CHECK(domain_1_nw_patch->child_ranks[2] == domain_2_nw_nw_patch->rank);
			CHECK(domain_1_nw_patch->child_ranks[3] == domain_2_nw_ne_patch->rank);

			CHECK(domain_1_ne_patch->child_ranks[0] == domain_2_ne_sw_patch->rank);
			CHECK(domain_1_ne_patch->child_ranks[1] == domain_2_ne_se_patch->rank);
			CHECK(domain_1_ne_patch->child_ranks[2] == domain_2_ne_nw_patch->rank);
			CHECK(domain_1_ne_patch->child_ranks[3] == domain_2_ne_ne_patch->rank);

			CHECK(domain_0_coarser_patch->child_ranks[0] == domain_1_sw_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[1] == domain_1_se_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[2] == domain_1_nw_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[3] == domain_1_ne_patch->rank);
		}
	}
	//SECTION("nbr_info ids are correct")
	{
		if (rank == 0) {
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_sw_se_patch->id);
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Side<2>::south()));
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_sw_nw_patch->id);

			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_sw_sw_patch->id);
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_se_sw_patch->id);
			CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Side<2>::south()));
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_sw_ne_patch->id);

			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_sw_se_patch->id);
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_se_se_patch->id);
			CHECK_FALSE(domain_2_se_sw_patch->hasNbr(Side<2>::south()));
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_se_nw_patch->id);

			CHECK(domain_2_se_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_se_sw_patch->id);
			CHECK_FALSE(domain_2_se_se_patch->hasNbr(Side<2>::east()));
			CHECK_FALSE(domain_2_se_se_patch->hasNbr(Side<2>::south()));
			CHECK(domain_2_se_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_se_ne_patch->id);

			CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_sw_ne_patch->id);
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_sw_sw_patch->id);
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_nw_sw_patch->id);

			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_sw_nw_patch->id);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_se_nw_patch->id);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_sw_se_patch->id);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_nw_se_patch->id);

			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_sw_ne_patch->id);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_se_ne_patch->id);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_se_sw_patch->id);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_ne_sw_patch->id);

			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_se_nw_patch->id);
			CHECK_FALSE(domain_2_se_ne_patch->hasNbr(Side<2>::east()));
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_se_se_patch->id);
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_ne_se_patch->id);

			CHECK_FALSE(domain_2_nw_sw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_nw_se_patch->id);
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_sw_nw_patch->id);
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_nw_nw_patch->id);

			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_nw_sw_patch->id);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_ne_sw_patch->id);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_sw_ne_patch->id);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_nw_ne_patch->id);

			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_nw_se_patch->id);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_ne_se_patch->id);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_se_nw_patch->id);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_ne_nw_patch->id);

			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_ne_sw_patch->id);
			CHECK_FALSE(domain_2_ne_se_patch->hasNbr(Side<2>::east()));
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_se_ne_patch->id);
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_ne_ne_patch->id);

			CHECK_FALSE(domain_2_nw_nw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_2_nw_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_nw_ne_patch->id);
			CHECK(domain_2_nw_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_nw_sw_patch->id);
			CHECK_FALSE(domain_2_nw_nw_patch->hasNbr(Side<2>::north()));

			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_nw_patch->rank);
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_nw_patch->rank);
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_nw_se_patch->rank);
			CHECK_FALSE(domain_2_nw_ne_patch->hasNbr(Side<2>::north()));

			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_ne_patch->rank);
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_ne_patch->rank);
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_ne_sw_patch->rank);
			CHECK_FALSE(domain_2_ne_nw_patch->hasNbr(Side<2>::north()));

			CHECK(domain_2_ne_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_ne_nw_patch->rank);
			CHECK_FALSE(domain_2_ne_ne_patch->hasNbr(Side<2>::east()));
			CHECK(domain_2_ne_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_ne_se_patch->rank);
			CHECK_FALSE(domain_2_ne_ne_patch->hasNbr(Side<2>::north()));

			CHECK_FALSE(domain_1_sw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_se_patch->rank);
			CHECK_FALSE(domain_1_sw_patch->hasNbr(Side<2>::south()));
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_nw_patch->rank);

			CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_sw_patch->rank);
			CHECK_FALSE(domain_1_se_patch->hasNbr(Side<2>::east()));
			CHECK_FALSE(domain_1_se_patch->hasNbr(Side<2>::south()));
			CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_ne_patch->rank);

			CHECK_FALSE(domain_1_nw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_ne_patch->rank);
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_sw_patch->rank);
			CHECK_FALSE(domain_1_nw_patch->hasNbr(Side<2>::north()));

			CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_nw_patch->rank);
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Side<2>::east()));
			CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_se_patch->rank);
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Side<2>::north()));

			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::west()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::east()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::south()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::north()));
		}
	}
	//SECTION("nbr_info ranks are correct")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_sw_se_patch->rank);
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_sw_nw_patch->rank);

			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_sw_sw_patch->rank);
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_se_sw_patch->rank);
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_sw_ne_patch->rank);

			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_sw_se_patch->rank);
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_se_se_patch->rank);
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_se_nw_patch->rank);

			CHECK(domain_2_se_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_se_sw_patch->rank);
			CHECK(domain_2_se_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_se_ne_patch->rank);

			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_sw_sw_patch->rank);
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_nw_sw_patch->rank);

			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_sw_nw_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_se_nw_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_sw_se_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_nw_se_patch->rank);

			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_se_ne_patch->rank);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_se_sw_patch->rank);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_ne_sw_patch->rank);

			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_se_nw_patch->rank);
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_se_se_patch->rank);
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_ne_se_patch->rank);

			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_nw_se_patch->rank);
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_sw_nw_patch->rank);
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_nw_nw_patch->rank);

			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_sw_patch->rank);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_sw_patch->rank);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_nw_ne_patch->rank);

			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_se_patch->rank);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_se_patch->rank);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_se_nw_patch->rank);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_ne_nw_patch->rank);

			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_ne_sw_patch->rank);
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_se_ne_patch->rank);
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_ne_ne_patch->rank);

			CHECK(domain_2_nw_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_nw_ne_patch->rank);
			CHECK(domain_2_nw_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_nw_sw_patch->rank);

			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_nw_patch->rank);
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_nw_patch->rank);
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_nw_se_patch->rank);

			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_ne_patch->rank);
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_ne_patch->rank);
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_ne_sw_patch->rank);

			CHECK(domain_2_ne_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_ne_nw_patch->rank);
			CHECK(domain_2_ne_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_ne_se_patch->rank);

			CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_se_patch->rank);
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_nw_patch->rank);

			CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_sw_patch->rank);
			CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_ne_patch->rank);

			CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_ne_patch->rank);
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_sw_patch->rank);

			CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_nw_patch->rank);
			CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_se_patch->rank);
		}
	}
	//SECTION("corner_nbr_info ids are correct")
	{
		if (rank == 0) {
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_sw_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_sw_ne_patch->id);

			CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_sw_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_sw_nw_patch->id);
			CHECK(domain_2_sw_se_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_se_nw_patch->id);

			CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_sw_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_sw_se_patch->id);
			CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_sw_nw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_nw_se_patch->id);

			CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_sw_sw_patch->id);
			CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_se_sw_patch->id);
			CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_nw_sw_patch->id);
			CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_ne_sw_patch->id);

			CHECK_FALSE(domain_2_se_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_2_se_sw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_se_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_sw_ne_patch->id);
			CHECK(domain_2_se_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_se_ne_patch->id);

			CHECK_FALSE(domain_2_se_se_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_2_se_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_se_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_se_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_se_nw_patch->id);
			CHECK_FALSE(domain_2_se_se_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_2_se_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_sw_se_patch->id);
			CHECK(domain_2_se_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_se_se_patch->id);
			CHECK(domain_2_se_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_nw_se_patch->id);
			CHECK(domain_2_se_nw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_ne_se_patch->id);

			CHECK(domain_2_se_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_se_sw_patch->id);
			CHECK_FALSE(domain_2_se_ne_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_se_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_ne_sw_patch->id);
			CHECK_FALSE(domain_2_se_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_2_nw_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_nw_sw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_sw_ne_patch->id);
			CHECK_FALSE(domain_2_nw_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_nw_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_nw_ne_patch->id);

			CHECK(domain_2_nw_se_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_sw_nw_patch->id);
			CHECK(domain_2_nw_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_se_nw_patch->id);
			CHECK(domain_2_nw_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_nw_nw_patch->id);
			CHECK(domain_2_nw_se_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_ne_nw_patch->id);

			CHECK_FALSE(domain_2_nw_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_nw_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_nw_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_nw_se_patch->id);
			CHECK_FALSE(domain_2_nw_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_2_nw_nw_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_2_nw_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_nw_sw_patch->id);
			CHECK(domain_2_nw_ne_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_ne_sw_patch->id);
			CHECK_FALSE(domain_2_nw_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_2_nw_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_2_ne_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_sw_ne_patch->id);
			CHECK(domain_2_ne_sw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_se_ne_patch->id);
			CHECK(domain_2_ne_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_nw_ne_patch->id);
			CHECK(domain_2_ne_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_ne_ne_patch->id);

			CHECK(domain_2_ne_se_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_se_nw_patch->id);
			CHECK_FALSE(domain_2_ne_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_ne_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_ne_nw_patch->id);
			CHECK_FALSE(domain_2_ne_se_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_2_ne_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_nw_se_patch->id);
			CHECK(domain_2_ne_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_ne_se_patch->id);
			CHECK_FALSE(domain_2_ne_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_2_ne_nw_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_2_ne_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_ne_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_ne_sw_patch->id);
			CHECK_FALSE(domain_2_ne_ne_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_2_ne_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_2_ne_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_1_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_1_ne_patch->id);

			CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_1_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_1_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_1_nw_patch->id);
			CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_1_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_1_se_patch->id);
			CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_1_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_1_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_1_sw_patch->id);
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::ne()));
		}
	}
	//SECTION("corner_nbr_info ranks are correct")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_sw_ne_patch->rank);

			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_sw_nw_patch->rank);
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_se_nw_patch->rank);

			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_sw_se_patch->rank);
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_nw_se_patch->rank);

			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_sw_sw_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_se_sw_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_nw_sw_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_ne_sw_patch->rank);

			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_se_ne_patch->rank);

			CHECK(domain_2_se_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_se_nw_patch->rank);

			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_sw_se_patch->rank);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_se_se_patch->rank);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_nw_se_patch->rank);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_ne_se_patch->rank);

			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_se_sw_patch->rank);
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_ne_sw_patch->rank);

			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_nw_ne_patch->rank);

			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_sw_nw_patch->rank);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_se_nw_patch->rank);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_nw_nw_patch->rank);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_ne_nw_patch->rank);

			CHECK(domain_2_nw_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_nw_se_patch->rank);

			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_nw_sw_patch->rank);
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_ne_sw_patch->rank);

			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_se_ne_patch->rank);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_nw_ne_patch->rank);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_ne_ne_patch->rank);

			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_se_nw_patch->rank);
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_ne_nw_patch->rank);

			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_nw_se_patch->rank);
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_ne_se_patch->rank);

			CHECK(domain_2_ne_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_ne_sw_patch->rank);

			CHECK(domain_1_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_1_ne_patch->rank);

			CHECK(domain_1_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_1_nw_patch->rank);

			CHECK(domain_1_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_1_se_patch->rank);

			CHECK(domain_1_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_1_sw_patch->rank);
		}
	}
}
TEST_CASE("P4estDomainGenerator 2x2 Refined SW", "[p4estDomGen]")
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	p4est_connectivity_t *conn = p4est_connectivity_new_unitsquare();

	p4est_t *p4est = p4est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

	p4est_refine(
	p4est, false,
	[](p4est_t *p4est, p4est_topidx_t witch_tree, p4est_quadrant_t *quadrant) -> int { return 1; },
	nullptr);
	p4est_refine(
	p4est, false,
	[](p4est_t *p4est, p4est_topidx_t witch_tree, p4est_quadrant_t *quadrant) -> int {
		return (quadrant->x == 0 && quadrant->y == 0);
	},
	nullptr);

	p4est_partition(p4est, true, nullptr);

	int    nx              = GENERATE(5, 10);
	int    ny              = GENERATE(5, 10);
	double scale_x         = GENERATE(0.5, 1.0);
	double scale_y         = GENERATE(0.5, 1.0);
	int    num_ghost_cells = GENERATE(0, 1, 2);

	INFO("nx: " << nx);
	INFO("ny: " << ny);
	INFO("scale_x: " << scale_x);
	INFO("scale_y: " << scale_y);
	INFO("num_ghost_cells: " << num_ghost_cells);

	P4estDomainGenerator::BlockMapFunc bmf
	= [&](int block_no, double unit_x, double unit_y, double &x, double &y) {
		  x = scale_x * unit_x;
		  y = scale_y * unit_y;
	  };

	P4estDomainGenerator dg(p4est, {nx, ny}, num_ghost_cells, bmf);

	auto domain_2 = dg.getFinestDomain();
	auto domain_1 = dg.getCoarserDomain();
	auto domain_0 = dg.getCoarserDomain();

	//SECTION("patches have correct ns")
	{
		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
		}
	}
	//SECTION("patches have ranks set")
	{
		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.rank == rank);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.rank == rank);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.rank == rank);
		}
	}

	//SECTION("patches have num_ghost_cells set")
	{
		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.num_ghost_cells == num_ghost_cells);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.num_ghost_cells == num_ghost_cells);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.num_ghost_cells == num_ghost_cells);
		}
	}

	std::vector<PatchInfo<2>> domain_2_patches = GetAllPatchesOnRank0(domain_2);
	std::vector<PatchInfo<2>> domain_1_patches = GetAllPatchesOnRank0(domain_1);
	std::vector<PatchInfo<2>> domain_0_patches = GetAllPatchesOnRank0(domain_0);

	const PatchInfo<2> *domain_2_sw_sw_patch = nullptr;
	const PatchInfo<2> *domain_2_sw_se_patch = nullptr;
	const PatchInfo<2> *domain_2_sw_nw_patch = nullptr;
	const PatchInfo<2> *domain_2_sw_ne_patch = nullptr;
	const PatchInfo<2> *domain_2_se_patch    = nullptr;
	const PatchInfo<2> *domain_2_nw_patch    = nullptr;
	const PatchInfo<2> *domain_2_ne_patch    = nullptr;

	const PatchInfo<2> *domain_1_sw_patch = nullptr;
	const PatchInfo<2> *domain_1_se_patch = nullptr;
	const PatchInfo<2> *domain_1_nw_patch = nullptr;
	const PatchInfo<2> *domain_1_ne_patch = nullptr;

	const PatchInfo<2> *domain_0_patch = nullptr;

	if (rank == 0) {
		for (PatchInfo<2> &patch : domain_2_patches) {
			double x = patch.starts[0];
			double y = patch.starts[1];
			if (x == Catch::Approx(0) && y == Catch::Approx(0)) {
				domain_2_sw_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0)) {
				domain_2_sw_se_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.25 * scale_y)) {
				domain_2_sw_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.25 * scale_y)) {
				domain_2_sw_ne_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0)) {
				domain_2_se_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y)) {
				domain_2_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_2_ne_patch = &patch;
			}
		}

		for (PatchInfo<2> &patch : domain_1_patches) {
			double x = patch.starts[0];
			double y = patch.starts[1];
			if (x == Catch::Approx(0) && y == Catch::Approx(0)) {
				domain_1_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0)) {
				domain_1_se_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y)) {
				domain_1_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_1_ne_patch = &patch;
			}
		}

		domain_0_patch = &domain_0_patches[0];

		REQUIRE(domain_2_sw_sw_patch != nullptr);
		REQUIRE(domain_2_sw_se_patch != nullptr);
		REQUIRE(domain_2_sw_nw_patch != nullptr);
		REQUIRE(domain_2_sw_ne_patch != nullptr);
		REQUIRE(domain_2_se_patch != nullptr);
		REQUIRE(domain_2_nw_patch != nullptr);
		REQUIRE(domain_2_ne_patch != nullptr);

		REQUIRE(domain_1_sw_patch != nullptr);
		REQUIRE(domain_1_se_patch != nullptr);
		REQUIRE(domain_1_nw_patch != nullptr);
		REQUIRE(domain_1_ne_patch != nullptr);

		CHECK(domain_0_patch->starts[0] == Catch::Approx(0.0));
		CHECK(domain_0_patch->starts[1] == Catch::Approx(0.0));
	}

	//SECTION("patches have correct spacings")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_sw_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_2_sw_sw_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_2_sw_se_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_2_sw_se_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_2_sw_nw_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_2_sw_nw_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_2_sw_ne_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_2_sw_ne_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));

			CHECK(domain_2_se_patch->spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
			CHECK(domain_2_se_patch->spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
			CHECK(domain_2_nw_patch->spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
			CHECK(domain_2_nw_patch->spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
			CHECK(domain_2_ne_patch->spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
			CHECK(domain_2_ne_patch->spacings[1] == Catch::Approx(scale_y * 0.5 / ny));

			CHECK(domain_1_sw_patch->spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
			CHECK(domain_1_sw_patch->spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
			CHECK(domain_1_se_patch->spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
			CHECK(domain_1_se_patch->spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
			CHECK(domain_1_nw_patch->spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
			CHECK(domain_1_nw_patch->spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
			CHECK(domain_1_ne_patch->spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
			CHECK(domain_1_ne_patch->spacings[1] == Catch::Approx(scale_y * 0.5 / ny));

			CHECK(domain_0_patch->spacings[0] == Catch::Approx(scale_x * 1.0 / nx));
			CHECK(domain_0_patch->spacings[1] == Catch::Approx(scale_y * 1.0 / ny));
		}
	}

	//SECTION("patches have refine_level set")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_sw_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_2_sw_sw_patch->refine_level == 2);
			CHECK(domain_2_sw_se_patch->refine_level == 2);
			CHECK(domain_2_sw_nw_patch->refine_level == 2);
			CHECK(domain_2_sw_ne_patch->refine_level == 2);

			CHECK(domain_2_se_patch->refine_level == 1);
			CHECK(domain_2_nw_patch->refine_level == 1);
			CHECK(domain_2_ne_patch->refine_level == 1);

			CHECK(domain_1_sw_patch->refine_level == 1);
			CHECK(domain_1_se_patch->refine_level == 1);
			CHECK(domain_1_nw_patch->refine_level == 1);
			CHECK(domain_1_ne_patch->refine_level == 1);

			CHECK(domain_0_patch->refine_level == 0);
		}
	}

	//SECTION("parent ids are set correctly")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_sw_patch->parent_id == domain_1_sw_patch->id);
			CHECK(domain_2_sw_se_patch->parent_id == domain_1_sw_patch->id);
			CHECK(domain_2_sw_nw_patch->parent_id == domain_1_sw_patch->id);
			CHECK(domain_2_sw_ne_patch->parent_id == domain_1_sw_patch->id);

			CHECK(domain_2_se_patch->parent_id == domain_1_se_patch->id);
			CHECK(domain_2_nw_patch->parent_id == domain_1_nw_patch->id);
			CHECK(domain_2_ne_patch->parent_id == domain_1_ne_patch->id);

			CHECK(domain_1_sw_patch->parent_id == domain_0_patch->id);
			CHECK(domain_1_se_patch->parent_id == domain_0_patch->id);
			CHECK(domain_1_nw_patch->parent_id == domain_0_patch->id);
			CHECK(domain_1_ne_patch->parent_id == domain_0_patch->id);

			CHECK(domain_0_patch->parent_id == -1);
		}
	}
	//SECTION("child ids are set correctly")
	{
		if (rank == 0) {
			for (int i = 0; i < 4; i++) {
				CHECK(domain_2_sw_sw_patch->child_ids[i] == -1);
				CHECK(domain_2_sw_se_patch->child_ids[i] == -1);
				CHECK(domain_2_sw_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_sw_ne_patch->child_ids[i] == -1);
				CHECK(domain_2_se_patch->child_ids[i] == -1);
				CHECK(domain_2_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_ne_patch->child_ids[i] == -1);
			}

			CHECK(domain_1_sw_patch->child_ids[0] == domain_2_sw_sw_patch->id);
			CHECK(domain_1_sw_patch->child_ids[1] == domain_2_sw_se_patch->id);
			CHECK(domain_1_sw_patch->child_ids[2] == domain_2_sw_nw_patch->id);
			CHECK(domain_1_sw_patch->child_ids[3] == domain_2_sw_ne_patch->id);

			CHECK(domain_1_se_patch->child_ids[0] == domain_2_se_patch->id);
			CHECK(domain_1_nw_patch->child_ids[0] == domain_2_nw_patch->id);
			CHECK(domain_1_ne_patch->child_ids[0] == domain_2_ne_patch->id);

			for (int i = 1; i < 4; i++) {
				CHECK(domain_1_se_patch->child_ids[i] == -1);
				CHECK(domain_1_nw_patch->child_ids[i] == -1);
				CHECK(domain_1_ne_patch->child_ids[i] == -1);
			}

			CHECK(domain_0_patch->child_ids[0] == domain_1_sw_patch->id);
			CHECK(domain_0_patch->child_ids[1] == domain_1_se_patch->id);
			CHECK(domain_0_patch->child_ids[2] == domain_1_nw_patch->id);
			CHECK(domain_0_patch->child_ids[3] == domain_1_ne_patch->id);
		}
	}
	//SECTION("orth on parent is set correctly")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_2_sw_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_2_sw_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_2_sw_ne_patch->orth_on_parent == Orthant<2>::ne());
			CHECK(domain_2_se_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_2_nw_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_2_ne_patch->orth_on_parent == Orthant<2>::null());

			CHECK(domain_1_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_1_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_1_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_1_ne_patch->orth_on_parent == Orthant<2>::ne());

			CHECK(domain_0_patch->orth_on_parent == Orthant<2>::null());
		}
	}
	//SECTION("parent ranks are set correctly")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_sw_patch->parent_rank == domain_1_sw_patch->rank);
			CHECK(domain_2_sw_se_patch->parent_rank == domain_1_sw_patch->rank);
			CHECK(domain_2_sw_nw_patch->parent_rank == domain_1_sw_patch->rank);
			CHECK(domain_2_sw_ne_patch->parent_rank == domain_1_sw_patch->rank);

			CHECK(domain_2_se_patch->parent_rank == domain_1_se_patch->rank);
			CHECK(domain_2_nw_patch->parent_rank == domain_1_nw_patch->rank);
			CHECK(domain_2_ne_patch->parent_rank == domain_1_ne_patch->rank);

			CHECK(domain_1_sw_patch->parent_rank == domain_0_patch->rank);
			CHECK(domain_1_se_patch->parent_rank == domain_0_patch->rank);
			CHECK(domain_1_nw_patch->parent_rank == domain_0_patch->rank);
			CHECK(domain_1_ne_patch->parent_rank == domain_0_patch->rank);

			CHECK(domain_0_patch->parent_rank == -1);
		}
	}
	//SECTION("child ranks are set correctly")
	{
		if (rank == 0) {
			for (int i = 0; i < 4; i++) {
				CHECK(domain_2_sw_sw_patch->child_ranks[i] == -1);
				CHECK(domain_2_sw_se_patch->child_ranks[i] == -1);
				CHECK(domain_2_sw_nw_patch->child_ranks[i] == -1);
				CHECK(domain_2_sw_ne_patch->child_ranks[i] == -1);
				CHECK(domain_2_se_patch->child_ranks[i] == -1);
				CHECK(domain_2_nw_patch->child_ranks[i] == -1);
				CHECK(domain_2_ne_patch->child_ranks[i] == -1);
			}

			CHECK(domain_1_sw_patch->child_ranks[0] == domain_2_sw_sw_patch->rank);
			CHECK(domain_1_sw_patch->child_ranks[1] == domain_2_sw_se_patch->rank);
			CHECK(domain_1_sw_patch->child_ranks[2] == domain_2_sw_nw_patch->rank);
			CHECK(domain_1_sw_patch->child_ranks[3] == domain_2_sw_ne_patch->rank);

			CHECK(domain_1_se_patch->child_ranks[0] == domain_2_se_patch->rank);
			CHECK(domain_1_nw_patch->child_ranks[0] == domain_2_nw_patch->rank);
			CHECK(domain_1_ne_patch->child_ranks[0] == domain_2_ne_patch->rank);

			for (int i = 1; i < 4; i++) {
				CHECK(domain_1_se_patch->child_ranks[i] == -1);
				CHECK(domain_1_nw_patch->child_ranks[i] == -1);
				CHECK(domain_1_ne_patch->child_ranks[i] == -1);
			}

			CHECK(domain_0_patch->child_ranks[0] == domain_1_sw_patch->rank);
			CHECK(domain_0_patch->child_ranks[1] == domain_1_se_patch->rank);
			CHECK(domain_0_patch->child_ranks[2] == domain_1_nw_patch->rank);
			CHECK(domain_0_patch->child_ranks[3] == domain_1_ne_patch->rank);
		}
	}
	//SECTION("nbr_info ids are correct")
	{
		if (rank == 0) {
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_sw_se_patch->id);
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Side<2>::south()));
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_sw_nw_patch->id);

			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_sw_sw_patch->id);
			CHECK(domain_2_sw_se_patch->getCoarseNbrInfo(Side<2>::east()).id == domain_2_se_patch->id);
			CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Side<2>::south()));
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_sw_ne_patch->id);

			CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_sw_ne_patch->id);
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_sw_sw_patch->id);
			CHECK(domain_2_sw_nw_patch->getCoarseNbrInfo(Side<2>::north()).id == domain_2_nw_patch->id);

			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_sw_nw_patch->id);
			CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Side<2>::east()).id == domain_2_se_patch->id);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_sw_se_patch->id);
			CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Side<2>::north()).id == domain_2_nw_patch->id);

			CHECK(domain_2_se_patch->getFineNbrInfo(Side<2>::west()).ids[0] == domain_2_sw_se_patch->id);
			CHECK(domain_2_se_patch->getFineNbrInfo(Side<2>::west()).ids[1] == domain_2_sw_ne_patch->id);
			CHECK_FALSE(domain_2_se_patch->hasNbr(Side<2>::east()));
			CHECK_FALSE(domain_2_se_patch->hasNbr(Side<2>::south()));
			CHECK(domain_2_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_ne_patch->id);

			CHECK_FALSE(domain_2_nw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_2_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_ne_patch->id);
			CHECK(domain_2_nw_patch->getFineNbrInfo(Side<2>::south()).ids[0] == domain_2_sw_nw_patch->id);
			CHECK(domain_2_nw_patch->getFineNbrInfo(Side<2>::south()).ids[1] == domain_2_sw_ne_patch->id);
			CHECK_FALSE(domain_2_nw_patch->hasNbr(Side<2>::north()));

			CHECK(domain_2_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_nw_patch->id);
			CHECK_FALSE(domain_2_ne_patch->hasNbr(Side<2>::east()));
			CHECK(domain_2_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_se_patch->id);
			CHECK_FALSE(domain_2_ne_patch->hasNbr(Side<2>::north()));

			CHECK_FALSE(domain_1_sw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_1_se_patch->id);
			CHECK_FALSE(domain_1_sw_patch->hasNbr(Side<2>::south()));
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_1_nw_patch->id);

			CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_1_sw_patch->id);
			CHECK_FALSE(domain_1_se_patch->hasNbr(Side<2>::east()));
			CHECK_FALSE(domain_1_se_patch->hasNbr(Side<2>::south()));
			CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_1_ne_patch->id);

			CHECK_FALSE(domain_1_nw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_1_ne_patch->id);
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_1_sw_patch->id);
			CHECK_FALSE(domain_1_nw_patch->hasNbr(Side<2>::north()));

			CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_1_nw_patch->id);
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Side<2>::east()));
			CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_1_se_patch->id);
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Side<2>::north()));

			CHECK_FALSE(domain_0_patch->hasNbr(Side<2>::west()));
			CHECK_FALSE(domain_0_patch->hasNbr(Side<2>::east()));
			CHECK_FALSE(domain_0_patch->hasNbr(Side<2>::south()));
			CHECK_FALSE(domain_0_patch->hasNbr(Side<2>::north()));
		}
	}
	//SECTION("nbr_info ranks are correct")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_sw_se_patch->rank);
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_sw_nw_patch->rank);

			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_sw_sw_patch->rank);
			CHECK(domain_2_sw_se_patch->getCoarseNbrInfo(Side<2>::east()).rank == domain_2_se_patch->rank);
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_sw_ne_patch->rank);

			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_sw_sw_patch->rank);
			CHECK(domain_2_sw_nw_patch->getCoarseNbrInfo(Side<2>::north()).rank == domain_2_nw_patch->rank);

			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_sw_nw_patch->rank);
			CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Side<2>::east()).rank == domain_2_se_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_sw_se_patch->rank);
			CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Side<2>::north()).rank == domain_2_nw_patch->rank);

			CHECK(domain_2_se_patch->getFineNbrInfo(Side<2>::west()).ranks[0] == domain_2_sw_se_patch->rank);
			CHECK(domain_2_se_patch->getFineNbrInfo(Side<2>::west()).ranks[1] == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_ne_patch->rank);

			CHECK(domain_2_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_patch->rank);
			CHECK(domain_2_nw_patch->getFineNbrInfo(Side<2>::south()).ranks[0] == domain_2_sw_nw_patch->rank);
			CHECK(domain_2_nw_patch->getFineNbrInfo(Side<2>::south()).ranks[1] == domain_2_sw_ne_patch->rank);

			CHECK(domain_2_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_patch->rank);
			CHECK(domain_2_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_se_patch->rank);

			CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_se_patch->rank);
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_nw_patch->rank);

			CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_sw_patch->rank);
			CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_ne_patch->rank);

			CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_ne_patch->rank);
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_sw_patch->rank);

			CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_nw_patch->rank);
			CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_se_patch->rank);
		}
	}
	//SECTION("nbr_info orth_on_coarse are correct")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_se_patch->getCoarseNbrInfo(Side<2>::east()).orth_on_coarse == Orthant<1>::lower());

			CHECK(domain_2_sw_nw_patch->getCoarseNbrInfo(Side<2>::north()).orth_on_coarse == Orthant<1>::lower());

			CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Side<2>::east()).orth_on_coarse == Orthant<1>::upper());
			CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Side<2>::north()).orth_on_coarse == Orthant<1>::upper());
		}
	}
	//SECTION("corner_nbr_info ids are correct")
	{
		if (rank == 0) {
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_sw_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_sw_ne_patch->id);

			CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_sw_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_sw_nw_patch->id);
			CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_sw_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_sw_se_patch->id);
			CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_sw_sw_patch->id);
			CHECK_FALSE(domain_2_sw_ne_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_2_sw_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Corner<2>::ne()).id == domain_2_ne_patch->id);

			CHECK_FALSE(domain_2_se_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_2_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_nw_patch->id);
			CHECK_FALSE(domain_2_se_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_2_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_se_patch->id);
			CHECK_FALSE(domain_2_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_2_nw_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_2_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_ne_patch->getFineNbrInfo(Corner<2>::sw()).ids[0] == domain_2_sw_ne_patch->id);
			CHECK_FALSE(domain_2_ne_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_2_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_2_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_1_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_1_ne_patch->id);

			CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_1_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_1_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_1_nw_patch->id);
			CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_1_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_1_se_patch->id);
			CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_1_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_1_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_1_sw_patch->id);
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_0_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_0_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_0_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_0_patch->hasNbr(Corner<2>::ne()));
		}
	}
	//SECTION("corner_nbr_info ranks are correct")
	{
		if (rank == 0) {
			CHECK(domain_2_sw_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_sw_ne_patch->rank);

			CHECK(domain_2_sw_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_sw_nw_patch->rank);

			CHECK(domain_2_sw_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_sw_se_patch->rank);

			CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_sw_sw_patch->rank);
			CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Corner<2>::ne()).rank == domain_2_ne_patch->rank);

			CHECK(domain_2_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_nw_patch->rank);

			CHECK(domain_2_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_se_patch->rank);

			CHECK(domain_2_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_ne_patch->getFineNbrInfo(Corner<2>::sw()).ranks[0] == domain_2_sw_ne_patch->rank);

			CHECK(domain_1_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_1_ne_patch->rank);

			CHECK(domain_1_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_1_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_1_nw_patch->rank);

			CHECK(domain_1_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_1_se_patch->rank);

			CHECK(domain_1_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_1_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_1_sw_patch->rank);
		}
	}
}
TEST_CASE("P4estDomainGenerator 4x4 Refined SW", "[p4estDomGen]")
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	p4est_connectivity_t *conn = p4est_connectivity_new_unitsquare();

	p4est_t *p4est = p4est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

	p4est_refine(
	p4est, false,
	[](p4est_t *p4est, p4est_topidx_t witch_tree, p4est_quadrant_t *quadrant) -> int { return 1; },
	nullptr);
	p4est_refine(
	p4est, false,
	[](p4est_t *p4est, p4est_topidx_t witch_tree, p4est_quadrant_t *quadrant) -> int { return 1; },
	nullptr);
	p4est_refine(
	p4est, false,
	[](p4est_t *p4est, p4est_topidx_t witch_tree, p4est_quadrant_t *quadrant) -> int {
		return (quadrant->x == 0 && quadrant->y == 0);
	},
	nullptr);

	p4est_partition(p4est, true, nullptr);

	int    nx              = GENERATE(5, 10);
	int    ny              = GENERATE(5, 10);
	double scale_x         = GENERATE(0.5, 1.0);
	double scale_y         = GENERATE(0.5, 1.0);
	int    num_ghost_cells = GENERATE(0, 1, 2);

	INFO("nx: " << nx);
	INFO("ny: " << ny);
	INFO("scale_x: " << scale_x);
	INFO("scale_y: " << scale_y);
	INFO("num_ghost_cells: " << num_ghost_cells);

	P4estDomainGenerator::BlockMapFunc bmf
	= [&](int block_no, double unit_x, double unit_y, double &x, double &y) {
		  x = scale_x * unit_x;
		  y = scale_y * unit_y;
	  };

	P4estDomainGenerator dg(p4est, {nx, ny}, num_ghost_cells, bmf);

	auto domain_3 = dg.getFinestDomain();
	auto domain_2 = dg.getCoarserDomain();
	auto domain_1 = dg.getCoarserDomain();
	auto domain_0 = dg.getCoarserDomain();

	//SECTION("correct number of patches")
	{
		CHECK(domain_3->getNumGlobalPatches() == 19);
		CHECK(domain_2->getNumGlobalPatches() == 16);
		CHECK(domain_1->getNumGlobalPatches() == 4);
		CHECK(domain_0->getNumGlobalPatches() == 1);
	}

	//SECTION("patches have correct ns")
	{
		for (auto patch : domain_3->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
		}

		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
		}
	}

	//SECTION("patches have ranks set")
	{
		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.rank == rank);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.rank == rank);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.rank == rank);
		}
	}
	//SECTION("patches have num_ghost_cells set")
	{
		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.num_ghost_cells == num_ghost_cells);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.num_ghost_cells == num_ghost_cells);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.num_ghost_cells == num_ghost_cells);
		}
	}

	std::vector<PatchInfo<2>> domain_3_patches = GetAllPatchesOnRank0(domain_3);
	std::vector<PatchInfo<2>> domain_2_patches = GetAllPatchesOnRank0(domain_2);
	std::vector<PatchInfo<2>> domain_1_patches = GetAllPatchesOnRank0(domain_1);
	std::vector<PatchInfo<2>> domain_0_patches = GetAllPatchesOnRank0(domain_0);

	INFO("nx: " << nx);
	INFO("ny: " << ny);
	INFO("scale_x: " << scale_x);
	INFO("scale_y: " << scale_y);
	INFO("num_ghost_cells: " << num_ghost_cells);

	const PatchInfo<2> *domain_3_sw_sw_sw_patch = nullptr;
	const PatchInfo<2> *domain_3_sw_sw_se_patch = nullptr;
	const PatchInfo<2> *domain_3_sw_sw_nw_patch = nullptr;
	const PatchInfo<2> *domain_3_sw_sw_ne_patch = nullptr;
	const PatchInfo<2> *domain_3_sw_se_patch    = nullptr;
	const PatchInfo<2> *domain_3_sw_nw_patch    = nullptr;
	const PatchInfo<2> *domain_3_sw_ne_patch    = nullptr;

	const PatchInfo<2> *domain_3_se_sw_patch = nullptr;
	const PatchInfo<2> *domain_3_se_se_patch = nullptr;
	const PatchInfo<2> *domain_3_se_nw_patch = nullptr;
	const PatchInfo<2> *domain_3_se_ne_patch = nullptr;

	const PatchInfo<2> *domain_3_nw_sw_patch = nullptr;
	const PatchInfo<2> *domain_3_nw_se_patch = nullptr;
	const PatchInfo<2> *domain_3_nw_nw_patch = nullptr;
	const PatchInfo<2> *domain_3_nw_ne_patch = nullptr;

	const PatchInfo<2> *domain_3_ne_sw_patch = nullptr;
	const PatchInfo<2> *domain_3_ne_se_patch = nullptr;
	const PatchInfo<2> *domain_3_ne_nw_patch = nullptr;
	const PatchInfo<2> *domain_3_ne_ne_patch = nullptr;

	const PatchInfo<2> *domain_2_sw_sw_patch = nullptr;
	const PatchInfo<2> *domain_2_sw_se_patch = nullptr;
	const PatchInfo<2> *domain_2_sw_nw_patch = nullptr;
	const PatchInfo<2> *domain_2_sw_ne_patch = nullptr;

	const PatchInfo<2> *domain_2_se_sw_patch = nullptr;
	const PatchInfo<2> *domain_2_se_se_patch = nullptr;
	const PatchInfo<2> *domain_2_se_nw_patch = nullptr;
	const PatchInfo<2> *domain_2_se_ne_patch = nullptr;

	const PatchInfo<2> *domain_2_nw_sw_patch = nullptr;
	const PatchInfo<2> *domain_2_nw_se_patch = nullptr;
	const PatchInfo<2> *domain_2_nw_nw_patch = nullptr;
	const PatchInfo<2> *domain_2_nw_ne_patch = nullptr;

	const PatchInfo<2> *domain_2_ne_sw_patch = nullptr;
	const PatchInfo<2> *domain_2_ne_se_patch = nullptr;
	const PatchInfo<2> *domain_2_ne_nw_patch = nullptr;
	const PatchInfo<2> *domain_2_ne_ne_patch = nullptr;

	const PatchInfo<2> *domain_1_sw_patch = nullptr;
	const PatchInfo<2> *domain_1_se_patch = nullptr;
	const PatchInfo<2> *domain_1_nw_patch = nullptr;
	const PatchInfo<2> *domain_1_ne_patch = nullptr;

	const PatchInfo<2> *domain_0_coarser_patch = nullptr;

	if (rank == 0) {
		for (PatchInfo<2> &patch : domain_3_patches) {
			double x = patch.starts[0];
			double y = patch.starts[1];
			if (x == Catch::Approx(0) && y == Catch::Approx(0)) {
				domain_3_sw_sw_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.125 * scale_x) && y == Catch::Approx(0)) {
				domain_3_sw_sw_se_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.125 * scale_y)) {
				domain_3_sw_sw_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.125 * scale_x) && y == Catch::Approx(0.125 * scale_y)) {
				domain_3_sw_sw_ne_patch = &patch;
			}

			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0)) {
				domain_3_sw_se_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0)) {
				domain_3_se_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.75 * scale_x) && y == Catch::Approx(0)) {
				domain_3_se_se_patch = &patch;
			}

			if (x == Catch::Approx(0) && y == Catch::Approx(0.25 * scale_y)) {
				domain_3_sw_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.25 * scale_y)) {
				domain_3_sw_ne_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.25 * scale_y)) {
				domain_3_se_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.75 * scale_x) && y == Catch::Approx(0.25 * scale_y)) {
				domain_3_se_ne_patch = &patch;
			}

			if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y)) {
				domain_3_nw_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_3_nw_se_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_3_ne_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.75 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_3_ne_se_patch = &patch;
			}

			if (x == Catch::Approx(0) && y == Catch::Approx(0.75 * scale_y)) {
				domain_3_nw_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.75 * scale_y)) {
				domain_3_nw_ne_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.75 * scale_y)) {
				domain_3_ne_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.75 * scale_x) && y == Catch::Approx(0.75 * scale_y)) {
				domain_3_ne_ne_patch = &patch;
			}
		}

		for (PatchInfo<2> &patch : domain_2_patches) {
			double x = patch.starts[0];
			double y = patch.starts[1];
			if (x == Catch::Approx(0) && y == Catch::Approx(0)) {
				domain_2_sw_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0)) {
				domain_2_sw_se_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0)) {
				domain_2_se_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.75 * scale_x) && y == Catch::Approx(0)) {
				domain_2_se_se_patch = &patch;
			}

			if (x == Catch::Approx(0) && y == Catch::Approx(0.25 * scale_y)) {
				domain_2_sw_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.25 * scale_y)) {
				domain_2_sw_ne_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.25 * scale_y)) {
				domain_2_se_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.75 * scale_x) && y == Catch::Approx(0.25 * scale_y)) {
				domain_2_se_ne_patch = &patch;
			}

			if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y)) {
				domain_2_nw_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_2_nw_se_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_2_ne_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.75 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_2_ne_se_patch = &patch;
			}

			if (x == Catch::Approx(0) && y == Catch::Approx(0.75 * scale_y)) {
				domain_2_nw_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.75 * scale_y)) {
				domain_2_nw_ne_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.75 * scale_y)) {
				domain_2_ne_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.75 * scale_x) && y == Catch::Approx(0.75 * scale_y)) {
				domain_2_ne_ne_patch = &patch;
			}
		}

		for (PatchInfo<2> &patch : domain_1_patches) {
			double x = patch.starts[0];
			double y = patch.starts[1];
			if (x == Catch::Approx(0) && y == Catch::Approx(0)) {
				domain_1_sw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0)) {
				domain_1_se_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y)) {
				domain_1_nw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
				domain_1_ne_patch = &patch;
			}
		}

		domain_0_coarser_patch = &domain_0_patches[0];

		REQUIRE(domain_3_sw_sw_sw_patch != nullptr);
		REQUIRE(domain_3_sw_sw_se_patch != nullptr);
		REQUIRE(domain_3_sw_sw_nw_patch != nullptr);
		REQUIRE(domain_3_sw_sw_ne_patch != nullptr);
		REQUIRE(domain_3_sw_se_patch != nullptr);
		REQUIRE(domain_3_sw_nw_patch != nullptr);
		REQUIRE(domain_3_sw_ne_patch != nullptr);

		REQUIRE(domain_3_se_sw_patch != nullptr);
		REQUIRE(domain_3_se_se_patch != nullptr);
		REQUIRE(domain_3_se_nw_patch != nullptr);
		REQUIRE(domain_3_se_ne_patch != nullptr);

		REQUIRE(domain_3_nw_sw_patch != nullptr);
		REQUIRE(domain_3_nw_se_patch != nullptr);
		REQUIRE(domain_3_nw_nw_patch != nullptr);
		REQUIRE(domain_3_nw_ne_patch != nullptr);

		REQUIRE(domain_3_ne_sw_patch != nullptr);
		REQUIRE(domain_3_ne_se_patch != nullptr);
		REQUIRE(domain_3_ne_nw_patch != nullptr);
		REQUIRE(domain_3_ne_ne_patch != nullptr);

		REQUIRE(domain_2_sw_sw_patch != nullptr);
		REQUIRE(domain_2_sw_se_patch != nullptr);
		REQUIRE(domain_2_sw_nw_patch != nullptr);
		REQUIRE(domain_2_sw_ne_patch != nullptr);

		REQUIRE(domain_2_se_sw_patch != nullptr);
		REQUIRE(domain_2_se_se_patch != nullptr);
		REQUIRE(domain_2_se_nw_patch != nullptr);
		REQUIRE(domain_2_se_ne_patch != nullptr);

		REQUIRE(domain_2_nw_sw_patch != nullptr);
		REQUIRE(domain_2_nw_se_patch != nullptr);
		REQUIRE(domain_2_nw_nw_patch != nullptr);
		REQUIRE(domain_2_nw_ne_patch != nullptr);

		REQUIRE(domain_2_ne_sw_patch != nullptr);
		REQUIRE(domain_2_ne_se_patch != nullptr);
		REQUIRE(domain_2_ne_nw_patch != nullptr);
		REQUIRE(domain_2_ne_ne_patch != nullptr);

		REQUIRE(domain_1_sw_patch != nullptr);
		REQUIRE(domain_1_se_patch != nullptr);
		REQUIRE(domain_1_nw_patch != nullptr);
		REQUIRE(domain_1_ne_patch != nullptr);

		CHECK(domain_0_coarser_patch->starts[0] == Catch::Approx(0.0));
		CHECK(domain_0_coarser_patch->starts[1] == Catch::Approx(0.0));
	}

	//SECTION("patches have correct spacings")
	{
		if (rank == 0) {
			CHECK(domain_3_sw_sw_sw_patch->spacings[0] == Catch::Approx(scale_x * 0.125 / nx));
			CHECK(domain_3_sw_sw_sw_patch->spacings[1] == Catch::Approx(scale_y * 0.125 / ny));
			CHECK(domain_3_sw_sw_se_patch->spacings[0] == Catch::Approx(scale_x * 0.125 / nx));
			CHECK(domain_3_sw_sw_se_patch->spacings[1] == Catch::Approx(scale_y * 0.125 / ny));
			CHECK(domain_3_sw_sw_nw_patch->spacings[0] == Catch::Approx(scale_x * 0.125 / nx));
			CHECK(domain_3_sw_sw_nw_patch->spacings[1] == Catch::Approx(scale_y * 0.125 / ny));
			CHECK(domain_3_sw_sw_ne_patch->spacings[0] == Catch::Approx(scale_x * 0.125 / nx));
			CHECK(domain_3_sw_sw_ne_patch->spacings[1] == Catch::Approx(scale_y * 0.125 / ny));

			CHECK(domain_3_sw_se_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_sw_se_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_3_sw_nw_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_sw_nw_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_3_sw_ne_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_sw_ne_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));

			CHECK(domain_3_se_sw_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_se_sw_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_3_se_se_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_se_se_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_3_se_nw_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_se_nw_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_3_se_ne_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_se_ne_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));

			CHECK(domain_3_nw_sw_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_nw_sw_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_3_nw_se_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_nw_se_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_3_nw_nw_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_nw_nw_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_3_nw_ne_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_nw_ne_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));

			CHECK(domain_3_ne_sw_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_ne_sw_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_3_ne_se_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_ne_se_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_3_ne_nw_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_ne_nw_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
			CHECK(domain_3_ne_ne_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(domain_3_ne_ne_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
		}
		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
			CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
			CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.spacings[0] == Catch::Approx(scale_x * 1.0 / nx));
			CHECK(patch.spacings[1] == Catch::Approx(scale_y * 1.0 / ny));
		}
	}

	//SECTION("patches have refine_level set")
	{
		if (rank == 0) {
			CHECK(domain_3_sw_sw_sw_patch->refine_level == 3);
			CHECK(domain_3_sw_sw_se_patch->refine_level == 3);
			CHECK(domain_3_sw_sw_nw_patch->refine_level == 3);
			CHECK(domain_3_sw_sw_ne_patch->refine_level == 3);

			CHECK(domain_3_sw_se_patch->refine_level == 2);
			CHECK(domain_3_sw_nw_patch->refine_level == 2);
			CHECK(domain_3_sw_ne_patch->refine_level == 2);

			CHECK(domain_3_se_sw_patch->refine_level == 2);
			CHECK(domain_3_se_se_patch->refine_level == 2);
			CHECK(domain_3_se_nw_patch->refine_level == 2);
			CHECK(domain_3_se_ne_patch->refine_level == 2);

			CHECK(domain_3_nw_sw_patch->refine_level == 2);
			CHECK(domain_3_nw_se_patch->refine_level == 2);
			CHECK(domain_3_nw_nw_patch->refine_level == 2);
			CHECK(domain_3_nw_ne_patch->refine_level == 2);

			CHECK(domain_3_ne_sw_patch->refine_level == 2);
			CHECK(domain_3_ne_se_patch->refine_level == 2);
			CHECK(domain_3_ne_nw_patch->refine_level == 2);
			CHECK(domain_3_ne_ne_patch->refine_level == 2);
		}

		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.refine_level == 2);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.refine_level == 1);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.refine_level == 0);
		}
	}

	//SECTION("parent ids are set correctly")
	{
		if (rank == 0) {
			CHECK(domain_3_sw_sw_sw_patch->parent_id == domain_2_sw_sw_patch->id);
			CHECK(domain_3_sw_sw_se_patch->parent_id == domain_2_sw_sw_patch->id);
			CHECK(domain_3_sw_sw_nw_patch->parent_id == domain_2_sw_sw_patch->id);
			CHECK(domain_3_sw_sw_ne_patch->parent_id == domain_2_sw_sw_patch->id);

			CHECK(domain_3_sw_se_patch->parent_id == domain_2_sw_se_patch->id);
			CHECK(domain_3_sw_nw_patch->parent_id == domain_2_sw_nw_patch->id);
			CHECK(domain_3_sw_ne_patch->parent_id == domain_2_sw_ne_patch->id);

			CHECK(domain_3_se_sw_patch->parent_id == domain_2_se_sw_patch->id);
			CHECK(domain_3_se_se_patch->parent_id == domain_2_se_se_patch->id);
			CHECK(domain_3_se_nw_patch->parent_id == domain_2_se_nw_patch->id);
			CHECK(domain_3_se_ne_patch->parent_id == domain_2_se_ne_patch->id);

			CHECK(domain_3_nw_sw_patch->parent_id == domain_2_nw_sw_patch->id);
			CHECK(domain_3_nw_se_patch->parent_id == domain_2_nw_se_patch->id);
			CHECK(domain_3_nw_nw_patch->parent_id == domain_2_nw_nw_patch->id);
			CHECK(domain_3_nw_ne_patch->parent_id == domain_2_nw_ne_patch->id);

			CHECK(domain_3_ne_sw_patch->parent_id == domain_2_ne_sw_patch->id);
			CHECK(domain_3_ne_se_patch->parent_id == domain_2_ne_se_patch->id);
			CHECK(domain_3_ne_nw_patch->parent_id == domain_2_ne_nw_patch->id);
			CHECK(domain_3_ne_ne_patch->parent_id == domain_2_ne_ne_patch->id);

			CHECK(domain_2_sw_sw_patch->parent_id == domain_1_sw_patch->id);
			CHECK(domain_2_sw_se_patch->parent_id == domain_1_sw_patch->id);
			CHECK(domain_2_sw_nw_patch->parent_id == domain_1_sw_patch->id);
			CHECK(domain_2_sw_ne_patch->parent_id == domain_1_sw_patch->id);

			CHECK(domain_2_se_sw_patch->parent_id == domain_1_se_patch->id);
			CHECK(domain_2_se_se_patch->parent_id == domain_1_se_patch->id);
			CHECK(domain_2_se_nw_patch->parent_id == domain_1_se_patch->id);
			CHECK(domain_2_se_ne_patch->parent_id == domain_1_se_patch->id);

			CHECK(domain_2_nw_sw_patch->parent_id == domain_1_nw_patch->id);
			CHECK(domain_2_nw_se_patch->parent_id == domain_1_nw_patch->id);
			CHECK(domain_2_nw_nw_patch->parent_id == domain_1_nw_patch->id);
			CHECK(domain_2_nw_ne_patch->parent_id == domain_1_nw_patch->id);

			CHECK(domain_2_ne_sw_patch->parent_id == domain_1_ne_patch->id);
			CHECK(domain_2_ne_se_patch->parent_id == domain_1_ne_patch->id);
			CHECK(domain_2_ne_nw_patch->parent_id == domain_1_ne_patch->id);
			CHECK(domain_2_ne_ne_patch->parent_id == domain_1_ne_patch->id);

			CHECK(domain_1_sw_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_se_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_nw_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_ne_patch->parent_id == domain_0_coarser_patch->id);

			CHECK(domain_0_coarser_patch->parent_id == -1);
		}
	}
	//SECTION("child ids are set correctly")
	{
		if (rank == 0) {
			for (int i = 0; i < 4; i++) {
				INFO("i: " << i);
				CHECK(domain_3_sw_sw_sw_patch->child_ids[i] == -1);
				CHECK(domain_3_sw_sw_se_patch->child_ids[i] == -1);
				CHECK(domain_3_sw_sw_nw_patch->child_ids[i] == -1);
				CHECK(domain_3_sw_sw_ne_patch->child_ids[i] == -1);

				CHECK(domain_3_sw_se_patch->child_ids[i] == -1);
				CHECK(domain_3_sw_nw_patch->child_ids[i] == -1);
				CHECK(domain_3_sw_ne_patch->child_ids[i] == -1);

				CHECK(domain_3_se_sw_patch->child_ids[i] == -1);
				CHECK(domain_3_se_se_patch->child_ids[i] == -1);
				CHECK(domain_3_se_nw_patch->child_ids[i] == -1);
				CHECK(domain_3_se_ne_patch->child_ids[i] == -1);

				CHECK(domain_3_nw_sw_patch->child_ids[i] == -1);
				CHECK(domain_3_nw_se_patch->child_ids[i] == -1);
				CHECK(domain_3_nw_nw_patch->child_ids[i] == -1);
				CHECK(domain_3_nw_ne_patch->child_ids[i] == -1);

				CHECK(domain_3_ne_sw_patch->child_ids[i] == -1);
				CHECK(domain_3_ne_se_patch->child_ids[i] == -1);
				CHECK(domain_3_ne_nw_patch->child_ids[i] == -1);
				CHECK(domain_3_ne_ne_patch->child_ids[i] == -1);
			}

			CHECK(domain_2_sw_sw_patch->child_ids[0] == domain_3_sw_sw_sw_patch->id);
			CHECK(domain_2_sw_sw_patch->child_ids[1] == domain_3_sw_sw_se_patch->id);
			CHECK(domain_2_sw_sw_patch->child_ids[2] == domain_3_sw_sw_nw_patch->id);
			CHECK(domain_2_sw_sw_patch->child_ids[3] == domain_3_sw_sw_ne_patch->id);

			CHECK(domain_2_sw_se_patch->child_ids[0] == domain_3_sw_se_patch->id);
			CHECK(domain_2_sw_nw_patch->child_ids[0] == domain_3_sw_nw_patch->id);
			CHECK(domain_2_sw_ne_patch->child_ids[0] == domain_3_sw_ne_patch->id);

			CHECK(domain_2_se_sw_patch->child_ids[0] == domain_3_se_sw_patch->id);
			CHECK(domain_2_se_se_patch->child_ids[0] == domain_3_se_se_patch->id);
			CHECK(domain_2_se_nw_patch->child_ids[0] == domain_3_se_nw_patch->id);
			CHECK(domain_2_se_ne_patch->child_ids[0] == domain_3_se_ne_patch->id);

			CHECK(domain_2_nw_sw_patch->child_ids[0] == domain_3_nw_sw_patch->id);
			CHECK(domain_2_nw_se_patch->child_ids[0] == domain_3_nw_se_patch->id);
			CHECK(domain_2_nw_nw_patch->child_ids[0] == domain_3_nw_nw_patch->id);
			CHECK(domain_2_nw_ne_patch->child_ids[0] == domain_3_nw_ne_patch->id);

			CHECK(domain_2_ne_sw_patch->child_ids[0] == domain_3_ne_sw_patch->id);
			CHECK(domain_2_ne_se_patch->child_ids[0] == domain_3_ne_se_patch->id);
			CHECK(domain_2_ne_nw_patch->child_ids[0] == domain_3_ne_nw_patch->id);
			CHECK(domain_2_ne_ne_patch->child_ids[0] == domain_3_ne_ne_patch->id);

			for (int i = 1; i < 4; i++) {
				CHECK(domain_2_sw_se_patch->child_ids[i] == -1);
				CHECK(domain_2_sw_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_sw_ne_patch->child_ids[i] == -1);

				CHECK(domain_2_se_sw_patch->child_ids[i] == -1);
				CHECK(domain_2_se_se_patch->child_ids[i] == -1);
				CHECK(domain_2_se_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_se_ne_patch->child_ids[i] == -1);

				CHECK(domain_2_nw_sw_patch->child_ids[i] == -1);
				CHECK(domain_2_nw_se_patch->child_ids[i] == -1);
				CHECK(domain_2_nw_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_nw_ne_patch->child_ids[i] == -1);

				CHECK(domain_2_ne_sw_patch->child_ids[i] == -1);
				CHECK(domain_2_ne_se_patch->child_ids[i] == -1);
				CHECK(domain_2_ne_nw_patch->child_ids[i] == -1);
				CHECK(domain_2_ne_ne_patch->child_ids[i] == -1);
			}

			CHECK(domain_1_sw_patch->child_ids[0] == domain_2_sw_sw_patch->id);
			CHECK(domain_1_sw_patch->child_ids[1] == domain_2_sw_se_patch->id);
			CHECK(domain_1_sw_patch->child_ids[2] == domain_2_sw_nw_patch->id);
			CHECK(domain_1_sw_patch->child_ids[3] == domain_2_sw_ne_patch->id);

			CHECK(domain_1_se_patch->child_ids[0] == domain_2_se_sw_patch->id);
			CHECK(domain_1_se_patch->child_ids[1] == domain_2_se_se_patch->id);
			CHECK(domain_1_se_patch->child_ids[2] == domain_2_se_nw_patch->id);
			CHECK(domain_1_se_patch->child_ids[3] == domain_2_se_ne_patch->id);

			CHECK(domain_1_nw_patch->child_ids[0] == domain_2_nw_sw_patch->id);
			CHECK(domain_1_nw_patch->child_ids[1] == domain_2_nw_se_patch->id);
			CHECK(domain_1_nw_patch->child_ids[2] == domain_2_nw_nw_patch->id);
			CHECK(domain_1_nw_patch->child_ids[3] == domain_2_nw_ne_patch->id);

			CHECK(domain_1_ne_patch->child_ids[0] == domain_2_ne_sw_patch->id);
			CHECK(domain_1_ne_patch->child_ids[1] == domain_2_ne_se_patch->id);
			CHECK(domain_1_ne_patch->child_ids[2] == domain_2_ne_nw_patch->id);
			CHECK(domain_1_ne_patch->child_ids[3] == domain_2_ne_ne_patch->id);

			CHECK(domain_0_coarser_patch->child_ids[0] == domain_1_sw_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[1] == domain_1_se_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[2] == domain_1_nw_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[3] == domain_1_ne_patch->id);
		}
	}
	//SECTION("orth on parent is set correctly")
	{
		if (rank == 0) {
			CHECK(domain_3_sw_sw_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_3_sw_sw_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_3_sw_sw_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_3_sw_sw_ne_patch->orth_on_parent == Orthant<2>::ne());

			CHECK(domain_3_sw_se_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_3_sw_nw_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_3_sw_ne_patch->orth_on_parent == Orthant<2>::null());

			CHECK(domain_3_se_sw_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_3_se_se_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_3_se_nw_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_3_se_ne_patch->orth_on_parent == Orthant<2>::null());

			CHECK(domain_3_nw_sw_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_3_nw_se_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_3_nw_nw_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_3_nw_ne_patch->orth_on_parent == Orthant<2>::null());

			CHECK(domain_3_ne_sw_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_3_ne_se_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_3_ne_nw_patch->orth_on_parent == Orthant<2>::null());
			CHECK(domain_3_ne_ne_patch->orth_on_parent == Orthant<2>::null());

			CHECK(domain_2_sw_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_2_sw_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_2_sw_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_2_sw_ne_patch->orth_on_parent == Orthant<2>::ne());

			CHECK(domain_2_se_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_2_se_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_2_se_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_2_se_ne_patch->orth_on_parent == Orthant<2>::ne());

			CHECK(domain_2_nw_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_2_nw_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_2_nw_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_2_nw_ne_patch->orth_on_parent == Orthant<2>::ne());

			CHECK(domain_2_ne_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_2_ne_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_2_ne_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_2_ne_ne_patch->orth_on_parent == Orthant<2>::ne());

			CHECK(domain_1_sw_patch->orth_on_parent == Orthant<2>::sw());
			CHECK(domain_1_se_patch->orth_on_parent == Orthant<2>::se());
			CHECK(domain_1_nw_patch->orth_on_parent == Orthant<2>::nw());
			CHECK(domain_1_ne_patch->orth_on_parent == Orthant<2>::ne());

			CHECK(domain_0_coarser_patch->orth_on_parent == Orthant<2>::null());
		}
	}
	//SECTION("parent ranks are set correctly")
	{
		if (rank == 0) {
			CHECK(domain_3_sw_sw_sw_patch->parent_rank == domain_2_sw_sw_patch->rank);
			CHECK(domain_3_sw_sw_se_patch->parent_rank == domain_2_sw_sw_patch->rank);
			CHECK(domain_3_sw_sw_nw_patch->parent_rank == domain_2_sw_sw_patch->rank);
			CHECK(domain_3_sw_sw_ne_patch->parent_rank == domain_2_sw_sw_patch->rank);

			CHECK(domain_3_sw_se_patch->parent_rank == domain_2_sw_se_patch->rank);
			CHECK(domain_3_sw_nw_patch->parent_rank == domain_2_sw_nw_patch->rank);
			CHECK(domain_3_sw_ne_patch->parent_rank == domain_2_sw_ne_patch->rank);

			CHECK(domain_3_se_sw_patch->parent_rank == domain_2_se_sw_patch->rank);
			CHECK(domain_3_se_se_patch->parent_rank == domain_2_se_se_patch->rank);
			CHECK(domain_3_se_nw_patch->parent_rank == domain_2_se_nw_patch->rank);
			CHECK(domain_3_se_ne_patch->parent_rank == domain_2_se_ne_patch->rank);

			CHECK(domain_3_nw_sw_patch->parent_rank == domain_2_nw_sw_patch->rank);
			CHECK(domain_3_nw_se_patch->parent_rank == domain_2_nw_se_patch->rank);
			CHECK(domain_3_nw_nw_patch->parent_rank == domain_2_nw_nw_patch->rank);
			CHECK(domain_3_nw_ne_patch->parent_rank == domain_2_nw_ne_patch->rank);

			CHECK(domain_3_ne_sw_patch->parent_rank == domain_2_ne_sw_patch->rank);
			CHECK(domain_3_ne_se_patch->parent_rank == domain_2_ne_se_patch->rank);
			CHECK(domain_3_ne_nw_patch->parent_rank == domain_2_ne_nw_patch->rank);
			CHECK(domain_3_ne_ne_patch->parent_rank == domain_2_ne_ne_patch->rank);

			CHECK(domain_2_sw_sw_patch->parent_rank == domain_1_sw_patch->rank);
			CHECK(domain_2_sw_se_patch->parent_rank == domain_1_sw_patch->rank);
			CHECK(domain_2_sw_nw_patch->parent_rank == domain_1_sw_patch->rank);
			CHECK(domain_2_sw_ne_patch->parent_rank == domain_1_sw_patch->rank);

			CHECK(domain_2_se_sw_patch->parent_rank == domain_1_se_patch->rank);
			CHECK(domain_2_se_se_patch->parent_rank == domain_1_se_patch->rank);
			CHECK(domain_2_se_nw_patch->parent_rank == domain_1_se_patch->rank);
			CHECK(domain_2_se_ne_patch->parent_rank == domain_1_se_patch->rank);

			CHECK(domain_2_nw_sw_patch->parent_rank == domain_1_nw_patch->rank);
			CHECK(domain_2_nw_se_patch->parent_rank == domain_1_nw_patch->rank);
			CHECK(domain_2_nw_nw_patch->parent_rank == domain_1_nw_patch->rank);
			CHECK(domain_2_nw_ne_patch->parent_rank == domain_1_nw_patch->rank);

			CHECK(domain_2_ne_sw_patch->parent_rank == domain_1_ne_patch->rank);
			CHECK(domain_2_ne_se_patch->parent_rank == domain_1_ne_patch->rank);
			CHECK(domain_2_ne_nw_patch->parent_rank == domain_1_ne_patch->rank);
			CHECK(domain_2_ne_ne_patch->parent_rank == domain_1_ne_patch->rank);

			CHECK(domain_1_sw_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_se_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_nw_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_ne_patch->parent_rank == domain_0_coarser_patch->rank);

			CHECK(domain_0_coarser_patch->parent_rank == -1);
		}
	}
	//SECTION("child ranks are set correctly")
	{
		if (rank == 0) {
			for (int i = 0; i < 4; i++) {
				INFO("i: " << i);
				CHECK(domain_3_sw_sw_sw_patch->child_ranks[i] == -1);
				CHECK(domain_3_sw_sw_se_patch->child_ranks[i] == -1);
				CHECK(domain_3_sw_sw_nw_patch->child_ranks[i] == -1);
				CHECK(domain_3_sw_sw_ne_patch->child_ranks[i] == -1);

				CHECK(domain_3_sw_se_patch->child_ranks[i] == -1);
				CHECK(domain_3_sw_nw_patch->child_ranks[i] == -1);
				CHECK(domain_3_sw_ne_patch->child_ranks[i] == -1);

				CHECK(domain_3_se_sw_patch->child_ranks[i] == -1);
				CHECK(domain_3_se_se_patch->child_ranks[i] == -1);
				CHECK(domain_3_se_nw_patch->child_ranks[i] == -1);
				CHECK(domain_3_se_ne_patch->child_ranks[i] == -1);

				CHECK(domain_3_nw_sw_patch->child_ranks[i] == -1);
				CHECK(domain_3_nw_se_patch->child_ranks[i] == -1);
				CHECK(domain_3_nw_nw_patch->child_ranks[i] == -1);
				CHECK(domain_3_nw_ne_patch->child_ranks[i] == -1);

				CHECK(domain_3_ne_sw_patch->child_ranks[i] == -1);
				CHECK(domain_3_ne_se_patch->child_ranks[i] == -1);
				CHECK(domain_3_ne_nw_patch->child_ranks[i] == -1);
				CHECK(domain_3_ne_ne_patch->child_ranks[i] == -1);
			}

			CHECK(domain_2_sw_sw_patch->child_ranks[0] == domain_3_sw_sw_sw_patch->rank);
			CHECK(domain_2_sw_sw_patch->child_ranks[1] == domain_3_sw_sw_se_patch->rank);
			CHECK(domain_2_sw_sw_patch->child_ranks[2] == domain_3_sw_sw_nw_patch->rank);
			CHECK(domain_2_sw_sw_patch->child_ranks[3] == domain_3_sw_sw_ne_patch->rank);

			CHECK(domain_2_sw_se_patch->child_ranks[0] == domain_3_sw_se_patch->rank);
			CHECK(domain_2_sw_nw_patch->child_ranks[0] == domain_3_sw_nw_patch->rank);
			CHECK(domain_2_sw_ne_patch->child_ranks[0] == domain_3_sw_ne_patch->rank);

			CHECK(domain_2_se_sw_patch->child_ranks[0] == domain_3_se_sw_patch->rank);
			CHECK(domain_2_se_se_patch->child_ranks[0] == domain_3_se_se_patch->rank);
			CHECK(domain_2_se_nw_patch->child_ranks[0] == domain_3_se_nw_patch->rank);
			CHECK(domain_2_se_ne_patch->child_ranks[0] == domain_3_se_ne_patch->rank);

			CHECK(domain_2_nw_sw_patch->child_ranks[0] == domain_3_nw_sw_patch->rank);
			CHECK(domain_2_nw_se_patch->child_ranks[0] == domain_3_nw_se_patch->rank);
			CHECK(domain_2_nw_nw_patch->child_ranks[0] == domain_3_nw_nw_patch->rank);
			CHECK(domain_2_nw_ne_patch->child_ranks[0] == domain_3_nw_ne_patch->rank);

			CHECK(domain_2_ne_sw_patch->child_ranks[0] == domain_3_ne_sw_patch->rank);
			CHECK(domain_2_ne_se_patch->child_ranks[0] == domain_3_ne_se_patch->rank);
			CHECK(domain_2_ne_nw_patch->child_ranks[0] == domain_3_ne_nw_patch->rank);
			CHECK(domain_2_ne_ne_patch->child_ranks[0] == domain_3_ne_ne_patch->rank);

			for (int i = 1; i < 4; i++) {
				CHECK(domain_2_sw_se_patch->child_ranks[i] == -1);
				CHECK(domain_2_sw_nw_patch->child_ranks[i] == -1);
				CHECK(domain_2_sw_ne_patch->child_ranks[i] == -1);

				CHECK(domain_2_se_sw_patch->child_ranks[i] == -1);
				CHECK(domain_2_se_se_patch->child_ranks[i] == -1);
				CHECK(domain_2_se_nw_patch->child_ranks[i] == -1);
				CHECK(domain_2_se_ne_patch->child_ranks[i] == -1);

				CHECK(domain_2_nw_sw_patch->child_ranks[i] == -1);
				CHECK(domain_2_nw_se_patch->child_ranks[i] == -1);
				CHECK(domain_2_nw_nw_patch->child_ranks[i] == -1);
				CHECK(domain_2_nw_ne_patch->child_ranks[i] == -1);

				CHECK(domain_2_ne_sw_patch->child_ranks[i] == -1);
				CHECK(domain_2_ne_se_patch->child_ranks[i] == -1);
				CHECK(domain_2_ne_nw_patch->child_ranks[i] == -1);
				CHECK(domain_2_ne_ne_patch->child_ranks[i] == -1);
			}

			CHECK(domain_1_sw_patch->child_ranks[0] == domain_2_sw_sw_patch->rank);
			CHECK(domain_1_sw_patch->child_ranks[1] == domain_2_sw_se_patch->rank);
			CHECK(domain_1_sw_patch->child_ranks[2] == domain_2_sw_nw_patch->rank);
			CHECK(domain_1_sw_patch->child_ranks[3] == domain_2_sw_ne_patch->rank);

			CHECK(domain_1_se_patch->child_ranks[0] == domain_2_se_sw_patch->rank);
			CHECK(domain_1_se_patch->child_ranks[1] == domain_2_se_se_patch->rank);
			CHECK(domain_1_se_patch->child_ranks[2] == domain_2_se_nw_patch->rank);
			CHECK(domain_1_se_patch->child_ranks[3] == domain_2_se_ne_patch->rank);

			CHECK(domain_1_nw_patch->child_ranks[0] == domain_2_nw_sw_patch->rank);
			CHECK(domain_1_nw_patch->child_ranks[1] == domain_2_nw_se_patch->rank);
			CHECK(domain_1_nw_patch->child_ranks[2] == domain_2_nw_nw_patch->rank);
			CHECK(domain_1_nw_patch->child_ranks[3] == domain_2_nw_ne_patch->rank);

			CHECK(domain_1_ne_patch->child_ranks[0] == domain_2_ne_sw_patch->rank);
			CHECK(domain_1_ne_patch->child_ranks[1] == domain_2_ne_se_patch->rank);
			CHECK(domain_1_ne_patch->child_ranks[2] == domain_2_ne_nw_patch->rank);
			CHECK(domain_1_ne_patch->child_ranks[3] == domain_2_ne_ne_patch->rank);

			CHECK(domain_0_coarser_patch->child_ranks[0] == domain_1_sw_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[1] == domain_1_se_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[2] == domain_1_nw_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[3] == domain_1_ne_patch->rank);
		}
	}
	//SECTION("nbr_info ids are correct")
	{
		if (rank == 0) {
			CHECK_FALSE(domain_3_sw_sw_sw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_3_sw_sw_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_3_sw_sw_se_patch->id);
			CHECK_FALSE(domain_3_sw_sw_sw_patch->hasNbr(Side<2>::south()));
			CHECK(domain_3_sw_sw_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_sw_sw_nw_patch->id);

			CHECK(domain_3_sw_sw_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_3_sw_sw_sw_patch->id);
			CHECK(domain_3_sw_sw_se_patch->getCoarseNbrInfo(Side<2>::east()).id == domain_3_sw_se_patch->id);
			CHECK_FALSE(domain_3_sw_sw_se_patch->hasNbr(Side<2>::south()));
			CHECK(domain_3_sw_sw_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_sw_sw_ne_patch->id);

			CHECK_FALSE(domain_3_sw_sw_nw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_3_sw_sw_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_3_sw_sw_ne_patch->id);
			CHECK(domain_3_sw_sw_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_3_sw_sw_sw_patch->id);
			CHECK(domain_3_sw_sw_nw_patch->getCoarseNbrInfo(Side<2>::north()).id == domain_3_sw_nw_patch->id);

			CHECK(domain_3_sw_sw_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_3_sw_sw_nw_patch->id);
			CHECK(domain_3_sw_sw_ne_patch->getCoarseNbrInfo(Side<2>::east()).id == domain_3_sw_se_patch->id);
			CHECK(domain_3_sw_sw_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_3_sw_sw_se_patch->id);
			CHECK(domain_3_sw_sw_ne_patch->getCoarseNbrInfo(Side<2>::north()).id == domain_3_sw_nw_patch->id);

			CHECK(domain_3_sw_se_patch->getFineNbrInfo(Side<2>::west()).ids[0] == domain_3_sw_sw_se_patch->id);
			CHECK(domain_3_sw_se_patch->getFineNbrInfo(Side<2>::west()).ids[1] == domain_3_sw_sw_ne_patch->id);
			CHECK(domain_3_sw_se_patch->getNormalNbrInfo(Side<2>::east()).id == domain_3_se_sw_patch->id);
			CHECK_FALSE(domain_3_sw_se_patch->hasNbr(Side<2>::south()));
			CHECK(domain_3_sw_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_sw_ne_patch->id);

			CHECK_FALSE(domain_3_sw_nw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_3_sw_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_3_sw_ne_patch->id);
			CHECK(domain_3_sw_nw_patch->getFineNbrInfo(Side<2>::south()).ids[0] == domain_3_sw_sw_nw_patch->id);
			CHECK(domain_3_sw_nw_patch->getFineNbrInfo(Side<2>::south()).ids[1] == domain_3_sw_sw_ne_patch->id);
			CHECK(domain_3_sw_nw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_nw_sw_patch->id);

			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_3_sw_nw_patch->id);
			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Side<2>::east()).id == domain_3_se_nw_patch->id);
			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_3_sw_se_patch->id);
			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_nw_se_patch->id);

			CHECK(domain_3_se_sw_patch->getNormalNbrInfo(Side<2>::west()).id == domain_3_sw_se_patch->id);
			CHECK(domain_3_se_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_3_se_se_patch->id);
			CHECK_FALSE(domain_3_se_sw_patch->hasNbr(Side<2>::south()));
			CHECK(domain_3_se_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_se_nw_patch->id);

			CHECK(domain_3_se_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_3_se_sw_patch->id);
			CHECK_FALSE(domain_3_se_se_patch->hasNbr(Side<2>::east()));
			CHECK_FALSE(domain_3_se_se_patch->hasNbr(Side<2>::south()));
			CHECK(domain_3_se_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_se_ne_patch->id);

			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Side<2>::west()).id == domain_3_sw_ne_patch->id);
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_3_se_ne_patch->id);
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_3_se_sw_patch->id);
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_ne_sw_patch->id);

			CHECK(domain_3_se_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_3_se_nw_patch->id);
			CHECK_FALSE(domain_3_se_ne_patch->hasNbr(Side<2>::east()));
			CHECK(domain_3_se_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_3_se_se_patch->id);
			CHECK(domain_3_se_ne_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_ne_se_patch->id);

			CHECK_FALSE(domain_3_nw_sw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_3_nw_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_3_nw_se_patch->id);
			CHECK(domain_3_nw_sw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_3_sw_nw_patch->id);
			CHECK(domain_3_nw_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_nw_nw_patch->id);

			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_3_nw_sw_patch->id);
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Side<2>::east()).id == domain_3_ne_sw_patch->id);
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Side<2>::south()).id == domain_3_sw_ne_patch->id);
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_nw_ne_patch->id);

			CHECK_FALSE(domain_3_nw_nw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_3_nw_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_3_nw_ne_patch->id);
			CHECK(domain_3_nw_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_3_nw_sw_patch->id);
			CHECK_FALSE(domain_3_nw_nw_patch->hasNbr(Side<2>::north()));

			CHECK(domain_3_nw_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_nw_nw_patch->rank);
			CHECK(domain_3_nw_ne_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_ne_nw_patch->rank);
			CHECK(domain_3_nw_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_nw_se_patch->rank);
			CHECK_FALSE(domain_3_nw_ne_patch->hasNbr(Side<2>::north()));

			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Side<2>::west()).id == domain_3_nw_se_patch->id);
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_3_ne_se_patch->id);
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_3_se_nw_patch->id);
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_ne_nw_patch->id);

			CHECK(domain_3_ne_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_3_ne_sw_patch->id);
			CHECK_FALSE(domain_3_ne_se_patch->hasNbr(Side<2>::east()));
			CHECK(domain_3_ne_se_patch->getNormalNbrInfo(Side<2>::south()).id == domain_3_se_ne_patch->id);
			CHECK(domain_3_ne_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_3_ne_ne_patch->id);

			CHECK(domain_3_ne_nw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_nw_ne_patch->rank);
			CHECK(domain_3_ne_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_ne_ne_patch->rank);
			CHECK(domain_3_ne_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_ne_sw_patch->rank);
			CHECK_FALSE(domain_3_ne_nw_patch->hasNbr(Side<2>::north()));

			CHECK(domain_3_ne_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_ne_nw_patch->rank);
			CHECK_FALSE(domain_3_ne_ne_patch->hasNbr(Side<2>::east()));
			CHECK(domain_3_ne_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_ne_se_patch->rank);
			CHECK_FALSE(domain_3_ne_ne_patch->hasNbr(Side<2>::north()));

			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_sw_se_patch->id);
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Side<2>::south()));
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_sw_nw_patch->id);

			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_sw_sw_patch->id);
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_se_sw_patch->id);
			CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Side<2>::south()));
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_sw_ne_patch->id);

			CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_sw_ne_patch->id);
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_sw_sw_patch->id);
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_nw_sw_patch->id);

			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_sw_nw_patch->id);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_se_nw_patch->id);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_sw_se_patch->id);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_nw_se_patch->id);

			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_sw_se_patch->id);
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_se_se_patch->id);
			CHECK_FALSE(domain_2_se_sw_patch->hasNbr(Side<2>::south()));
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_se_nw_patch->id);

			CHECK(domain_2_se_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_se_sw_patch->id);
			CHECK_FALSE(domain_2_se_se_patch->hasNbr(Side<2>::east()));
			CHECK_FALSE(domain_2_se_se_patch->hasNbr(Side<2>::south()));
			CHECK(domain_2_se_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_se_ne_patch->id);

			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_sw_ne_patch->id);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_se_ne_patch->id);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_se_sw_patch->id);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_ne_sw_patch->id);

			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_se_nw_patch->id);
			CHECK_FALSE(domain_2_se_ne_patch->hasNbr(Side<2>::east()));
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_se_se_patch->id);
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_ne_se_patch->id);

			CHECK_FALSE(domain_2_nw_sw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_nw_se_patch->id);
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_sw_nw_patch->id);
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_nw_nw_patch->id);

			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_nw_sw_patch->id);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_ne_sw_patch->id);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_sw_ne_patch->id);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_nw_ne_patch->id);

			CHECK_FALSE(domain_2_nw_nw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_2_nw_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_nw_ne_patch->id);
			CHECK(domain_2_nw_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_nw_sw_patch->id);
			CHECK_FALSE(domain_2_nw_nw_patch->hasNbr(Side<2>::north()));

			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_nw_patch->rank);
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_nw_patch->rank);
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_nw_se_patch->rank);
			CHECK_FALSE(domain_2_nw_ne_patch->hasNbr(Side<2>::north()));

			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_nw_se_patch->id);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_ne_se_patch->id);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_se_nw_patch->id);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_ne_nw_patch->id);

			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_ne_sw_patch->id);
			CHECK_FALSE(domain_2_ne_se_patch->hasNbr(Side<2>::east()));
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_se_ne_patch->id);
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_ne_ne_patch->id);

			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_ne_patch->rank);
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_ne_patch->rank);
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_ne_sw_patch->rank);
			CHECK_FALSE(domain_2_ne_nw_patch->hasNbr(Side<2>::north()));

			CHECK(domain_2_ne_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_ne_nw_patch->rank);
			CHECK_FALSE(domain_2_ne_ne_patch->hasNbr(Side<2>::east()));
			CHECK(domain_2_ne_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_ne_se_patch->rank);
			CHECK_FALSE(domain_2_ne_ne_patch->hasNbr(Side<2>::north()));

			CHECK_FALSE(domain_1_sw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_se_patch->rank);
			CHECK_FALSE(domain_1_sw_patch->hasNbr(Side<2>::south()));
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_nw_patch->rank);

			CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_sw_patch->rank);
			CHECK_FALSE(domain_1_se_patch->hasNbr(Side<2>::east()));
			CHECK_FALSE(domain_1_se_patch->hasNbr(Side<2>::south()));
			CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_ne_patch->rank);

			CHECK_FALSE(domain_1_nw_patch->hasNbr(Side<2>::west()));
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_ne_patch->rank);
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_sw_patch->rank);
			CHECK_FALSE(domain_1_nw_patch->hasNbr(Side<2>::north()));

			CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_nw_patch->rank);
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Side<2>::east()));
			CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_se_patch->rank);
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Side<2>::north()));

			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::west()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::east()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::south()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::north()));
		}
	}
	//SECTION("nbr_info ranks are correct")
	{
		if (rank == 0) {
			CHECK(domain_3_sw_sw_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_sw_sw_se_patch->rank);
			CHECK(domain_3_sw_sw_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_sw_sw_nw_patch->rank);

			CHECK(domain_3_sw_sw_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_sw_sw_sw_patch->rank);
			CHECK(domain_3_sw_sw_se_patch->getCoarseNbrInfo(Side<2>::east()).rank == domain_3_sw_se_patch->rank);
			CHECK(domain_3_sw_sw_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_sw_sw_ne_patch->rank);

			CHECK(domain_3_sw_sw_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_sw_sw_ne_patch->rank);
			CHECK(domain_3_sw_sw_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_sw_sw_sw_patch->rank);
			CHECK(domain_3_sw_sw_nw_patch->getCoarseNbrInfo(Side<2>::north()).rank == domain_3_sw_nw_patch->rank);

			CHECK(domain_3_sw_sw_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_sw_sw_nw_patch->rank);
			CHECK(domain_3_sw_sw_ne_patch->getCoarseNbrInfo(Side<2>::east()).rank == domain_3_sw_nw_patch->rank);
			CHECK(domain_3_sw_sw_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_sw_sw_se_patch->rank);
			CHECK(domain_3_sw_sw_ne_patch->getCoarseNbrInfo(Side<2>::north()).rank == domain_3_sw_nw_patch->rank);

			CHECK(domain_3_sw_se_patch->getFineNbrInfo(Side<2>::west()).ranks[0] == domain_3_sw_sw_se_patch->rank);
			CHECK(domain_3_sw_se_patch->getFineNbrInfo(Side<2>::west()).ranks[1] == domain_3_sw_sw_ne_patch->rank);
			CHECK(domain_3_sw_se_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_se_sw_patch->rank);
			CHECK(domain_3_sw_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_sw_ne_patch->rank);

			CHECK(domain_3_sw_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_sw_ne_patch->rank);
			CHECK(domain_3_sw_nw_patch->getFineNbrInfo(Side<2>::south()).ranks[0] == domain_3_sw_sw_nw_patch->rank);
			CHECK(domain_3_sw_nw_patch->getFineNbrInfo(Side<2>::south()).ranks[1] == domain_3_sw_sw_ne_patch->rank);
			CHECK(domain_3_sw_nw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_nw_sw_patch->rank);

			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_sw_nw_patch->rank);
			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_se_nw_patch->rank);
			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_sw_se_patch->rank);
			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_nw_se_patch->rank);

			CHECK(domain_3_se_sw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_sw_se_patch->rank);
			CHECK(domain_3_se_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_se_se_patch->rank);
			CHECK(domain_3_se_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_se_nw_patch->rank);

			CHECK(domain_3_se_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_se_sw_patch->rank);
			CHECK(domain_3_se_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_se_ne_patch->rank);

			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_sw_ne_patch->rank);
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_se_ne_patch->rank);
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_se_sw_patch->rank);
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_ne_sw_patch->rank);

			CHECK(domain_3_se_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_se_nw_patch->rank);
			CHECK(domain_3_se_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_se_se_patch->rank);
			CHECK(domain_3_se_ne_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_ne_se_patch->rank);

			CHECK(domain_3_nw_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_nw_se_patch->rank);
			CHECK(domain_3_nw_sw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_sw_nw_patch->rank);
			CHECK(domain_3_nw_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_nw_nw_patch->rank);

			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_nw_sw_patch->rank);
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_ne_sw_patch->rank);
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_sw_ne_patch->rank);
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_nw_ne_patch->rank);

			CHECK(domain_3_nw_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_nw_ne_patch->rank);
			CHECK(domain_3_nw_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_nw_sw_patch->rank);

			CHECK(domain_3_nw_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_nw_nw_patch->rank);
			CHECK(domain_3_nw_ne_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_ne_nw_patch->rank);
			CHECK(domain_3_nw_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_nw_se_patch->rank);

			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_nw_se_patch->rank);
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_ne_se_patch->rank);
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_se_nw_patch->rank);
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_ne_nw_patch->rank);

			CHECK(domain_3_ne_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_ne_sw_patch->rank);
			CHECK(domain_3_ne_se_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_se_ne_patch->rank);
			CHECK(domain_3_ne_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_3_ne_ne_patch->rank);

			CHECK(domain_3_ne_nw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_nw_ne_patch->rank);
			CHECK(domain_3_ne_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_3_ne_ne_patch->rank);
			CHECK(domain_3_ne_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_ne_sw_patch->rank);

			CHECK(domain_3_ne_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_3_ne_nw_patch->rank);
			CHECK(domain_3_ne_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_3_ne_se_patch->rank);

			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_sw_se_patch->rank);
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_sw_nw_patch->rank);

			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_sw_sw_patch->rank);
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_se_sw_patch->rank);
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_sw_ne_patch->rank);

			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_sw_se_patch->rank);
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_se_se_patch->rank);
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_se_nw_patch->rank);

			CHECK(domain_2_se_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_se_sw_patch->rank);
			CHECK(domain_2_se_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_se_ne_patch->rank);

			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_sw_sw_patch->rank);
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_nw_sw_patch->rank);

			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_sw_nw_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_se_nw_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_sw_se_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_nw_se_patch->rank);

			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_se_ne_patch->rank);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_se_sw_patch->rank);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_ne_sw_patch->rank);

			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_se_nw_patch->rank);
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_se_se_patch->rank);
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_ne_se_patch->rank);

			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_nw_se_patch->rank);
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_sw_nw_patch->rank);
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_nw_nw_patch->rank);

			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_sw_patch->rank);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_sw_patch->rank);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_nw_ne_patch->rank);

			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_se_patch->rank);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_se_patch->rank);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_se_nw_patch->rank);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_ne_nw_patch->rank);

			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_ne_sw_patch->rank);
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_se_ne_patch->rank);
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_ne_ne_patch->rank);

			CHECK(domain_2_nw_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_nw_ne_patch->rank);
			CHECK(domain_2_nw_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_nw_sw_patch->rank);

			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_nw_patch->rank);
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_nw_patch->rank);
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_nw_se_patch->rank);

			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_ne_patch->rank);
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_ne_patch->rank);
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_ne_sw_patch->rank);

			CHECK(domain_2_ne_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_ne_nw_patch->rank);
			CHECK(domain_2_ne_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_ne_se_patch->rank);

			CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_se_patch->rank);
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_nw_patch->rank);

			CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_sw_patch->rank);
			CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_ne_patch->rank);

			CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_ne_patch->rank);
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_sw_patch->rank);

			CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_nw_patch->rank);
			CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_se_patch->rank);
		}
	}
	//SECTION("corner_nbr_info ids are correct")
	{
		if (rank == 0) {
			CHECK_FALSE(domain_3_sw_sw_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_3_sw_sw_sw_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_3_sw_sw_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_sw_sw_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_3_sw_sw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_3_sw_sw_ne_patch->id);

			CHECK_FALSE(domain_3_sw_sw_se_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_3_sw_sw_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_sw_sw_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_sw_sw_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_3_sw_sw_nw_patch->id);
			CHECK_FALSE(domain_3_sw_sw_se_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_3_sw_sw_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_sw_sw_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_sw_sw_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_3_sw_sw_se_patch->id);
			CHECK_FALSE(domain_3_sw_sw_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_3_sw_sw_nw_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_3_sw_sw_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_sw_sw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_3_sw_sw_sw_patch->id);
			CHECK_FALSE(domain_3_sw_sw_ne_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_3_sw_sw_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_sw_sw_ne_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_3_sw_sw_ne_patch->getCoarseNbrInfo(Corner<2>::ne()).id == domain_3_sw_ne_patch->id);

			CHECK_FALSE(domain_3_sw_se_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_3_sw_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_sw_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_sw_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_3_sw_nw_patch->id);
			CHECK(domain_3_sw_se_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_3_sw_se_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_3_se_nw_patch->id);

			CHECK_FALSE(domain_3_sw_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_sw_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_sw_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_3_sw_se_patch->id);
			CHECK_FALSE(domain_3_sw_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_sw_nw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_3_sw_nw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_3_nw_se_patch->id);

			CHECK(domain_3_sw_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_sw_ne_patch->getFineNbrInfo(Corner<2>::sw()).ids[0] == domain_3_sw_sw_ne_patch->id);
			CHECK(domain_3_sw_ne_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_3_se_sw_patch->id);
			CHECK(domain_3_sw_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_3_nw_sw_patch->id);
			CHECK(domain_3_sw_ne_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_3_ne_sw_patch->id);

			CHECK_FALSE(domain_3_se_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_3_se_sw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_se_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_se_sw_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_3_sw_ne_patch->id);
			CHECK(domain_3_se_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_3_se_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_3_se_ne_patch->id);

			CHECK_FALSE(domain_3_se_se_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_3_se_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_se_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_se_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_3_se_nw_patch->id);
			CHECK_FALSE(domain_3_se_se_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_3_se_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_3_sw_se_patch->id);
			CHECK(domain_3_se_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_3_se_se_patch->id);
			CHECK(domain_3_se_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_3_nw_se_patch->id);
			CHECK(domain_3_se_nw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_3_ne_se_patch->id);

			CHECK(domain_3_se_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_se_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_3_se_sw_patch->id);
			CHECK_FALSE(domain_3_se_ne_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_se_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_se_ne_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_3_ne_sw_patch->id);
			CHECK_FALSE(domain_3_se_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_3_nw_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_nw_sw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_nw_sw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_3_sw_ne_patch->id);
			CHECK_FALSE(domain_3_nw_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_nw_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_3_nw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_3_nw_ne_patch->id);

			CHECK(domain_3_nw_se_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_3_sw_nw_patch->id);
			CHECK(domain_3_nw_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_3_se_nw_patch->id);
			CHECK(domain_3_nw_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_3_nw_nw_patch->id);
			CHECK(domain_3_nw_se_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_3_ne_nw_patch->id);

			CHECK_FALSE(domain_3_nw_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_nw_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_nw_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_3_nw_se_patch->id);
			CHECK_FALSE(domain_3_nw_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_3_nw_nw_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_3_nw_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_nw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_3_nw_sw_patch->id);
			CHECK(domain_3_nw_ne_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_nw_ne_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_3_ne_sw_patch->id);
			CHECK_FALSE(domain_3_nw_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_3_nw_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_3_ne_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_3_sw_ne_patch->id);
			CHECK(domain_3_ne_sw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_3_se_ne_patch->id);
			CHECK(domain_3_ne_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_3_nw_ne_patch->id);
			CHECK(domain_3_ne_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_3_ne_ne_patch->id);

			CHECK(domain_3_ne_se_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_ne_se_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_3_se_nw_patch->id);
			CHECK_FALSE(domain_3_ne_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_ne_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_3_ne_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_3_ne_nw_patch->id);
			CHECK_FALSE(domain_3_ne_se_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_3_ne_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_ne_nw_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_3_nw_se_patch->id);
			CHECK(domain_3_ne_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_3_ne_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_3_ne_se_patch->id);
			CHECK_FALSE(domain_3_ne_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_3_ne_nw_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_3_ne_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_3_ne_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_3_ne_sw_patch->id);
			CHECK_FALSE(domain_3_ne_ne_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_3_ne_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_3_ne_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_sw_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_sw_ne_patch->id);

			CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_sw_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_sw_nw_patch->id);
			CHECK(domain_2_sw_se_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_se_nw_patch->id);

			CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_sw_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_sw_se_patch->id);
			CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_sw_nw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_nw_se_patch->id);

			CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_sw_sw_patch->id);
			CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_se_sw_patch->id);
			CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_nw_sw_patch->id);
			CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_ne_sw_patch->id);

			CHECK_FALSE(domain_2_se_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_2_se_sw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_se_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_sw_ne_patch->id);
			CHECK(domain_2_se_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_se_ne_patch->id);

			CHECK_FALSE(domain_2_se_se_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_2_se_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_se_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_se_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_se_nw_patch->id);
			CHECK_FALSE(domain_2_se_se_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_2_se_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_sw_se_patch->id);
			CHECK(domain_2_se_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_se_se_patch->id);
			CHECK(domain_2_se_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_nw_se_patch->id);
			CHECK(domain_2_se_nw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_ne_se_patch->id);

			CHECK(domain_2_se_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_se_sw_patch->id);
			CHECK_FALSE(domain_2_se_ne_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_se_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_ne_sw_patch->id);
			CHECK_FALSE(domain_2_se_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_2_nw_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_nw_sw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_sw_ne_patch->id);
			CHECK_FALSE(domain_2_nw_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_nw_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_nw_ne_patch->id);

			CHECK(domain_2_nw_se_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_sw_nw_patch->id);
			CHECK(domain_2_nw_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_se_nw_patch->id);
			CHECK(domain_2_nw_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_nw_nw_patch->id);
			CHECK(domain_2_nw_se_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_ne_nw_patch->id);

			CHECK_FALSE(domain_2_nw_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_nw_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_nw_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_nw_se_patch->id);
			CHECK_FALSE(domain_2_nw_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_2_nw_nw_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_2_nw_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_nw_sw_patch->id);
			CHECK(domain_2_nw_ne_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_ne_sw_patch->id);
			CHECK_FALSE(domain_2_nw_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_2_nw_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_2_ne_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_sw_ne_patch->id);
			CHECK(domain_2_ne_sw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_se_ne_patch->id);
			CHECK(domain_2_ne_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_nw_ne_patch->id);
			CHECK(domain_2_ne_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_ne_ne_patch->id);

			CHECK(domain_2_ne_se_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_se_nw_patch->id);
			CHECK_FALSE(domain_2_ne_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_ne_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_ne_nw_patch->id);
			CHECK_FALSE(domain_2_ne_se_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_2_ne_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_nw_se_patch->id);
			CHECK(domain_2_ne_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_ne_se_patch->id);
			CHECK_FALSE(domain_2_ne_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_2_ne_nw_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_2_ne_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_2_ne_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_ne_sw_patch->id);
			CHECK_FALSE(domain_2_ne_ne_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_2_ne_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_2_ne_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_1_sw_patch->hasNbr(Corner<2>::ne()));
			CHECK(domain_1_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_1_ne_patch->id);

			CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_1_se_patch->hasNbr(Corner<2>::nw()));
			CHECK(domain_1_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_1_nw_patch->id);
			CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_1_nw_patch->hasNbr(Corner<2>::se()));
			CHECK(domain_1_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_1_se_patch->id);
			CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::ne()));

			CHECK(domain_1_ne_patch->hasNbr(Corner<2>::sw()));
			CHECK(domain_1_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_1_sw_patch->id);
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::ne()));

			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::sw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::se()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::nw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::ne()));
		}
	}
	//SECTION("corner_nbr_info ranks are correct")
	{
		if (rank == 0) {
			CHECK(domain_3_sw_sw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_3_sw_sw_ne_patch->rank);

			CHECK(domain_3_sw_sw_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_3_sw_sw_nw_patch->rank);

			CHECK(domain_3_sw_sw_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_3_sw_sw_se_patch->rank);

			CHECK(domain_3_sw_sw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_3_sw_sw_sw_patch->rank);
			CHECK(domain_3_sw_sw_ne_patch->getCoarseNbrInfo(Corner<2>::ne()).rank == domain_3_sw_ne_patch->rank);

			CHECK(domain_3_sw_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_3_sw_nw_patch->rank);
			CHECK(domain_3_sw_se_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_3_se_nw_patch->rank);

			CHECK(domain_3_sw_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_3_sw_se_patch->rank);
			CHECK(domain_3_sw_nw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_3_nw_se_patch->rank);

			CHECK(domain_3_sw_ne_patch->getFineNbrInfo(Corner<2>::sw()).ranks[0] == domain_3_sw_sw_ne_patch->rank);
			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_3_se_sw_patch->rank);
			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_3_nw_sw_patch->rank);
			CHECK(domain_3_sw_ne_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_3_ne_sw_patch->rank);

			CHECK(domain_3_se_sw_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_3_sw_ne_patch->rank);
			CHECK(domain_3_se_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_3_se_ne_patch->rank);

			CHECK(domain_3_se_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_3_se_nw_patch->rank);

			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_3_sw_se_patch->rank);
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_3_se_se_patch->rank);
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_3_nw_se_patch->rank);
			CHECK(domain_3_se_nw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_3_ne_se_patch->rank);

			CHECK(domain_3_se_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_3_se_sw_patch->rank);
			CHECK(domain_3_se_ne_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_3_ne_sw_patch->rank);

			CHECK(domain_3_nw_sw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_3_sw_ne_patch->rank);
			CHECK(domain_3_nw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_3_nw_ne_patch->rank);

			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_3_sw_nw_patch->rank);
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_3_se_nw_patch->rank);
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_3_nw_nw_patch->rank);
			CHECK(domain_3_nw_se_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_3_ne_nw_patch->rank);

			CHECK(domain_3_nw_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_3_nw_se_patch->rank);

			CHECK(domain_3_nw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_3_nw_sw_patch->rank);
			CHECK(domain_3_nw_ne_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_3_ne_sw_patch->rank);

			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_3_sw_ne_patch->rank);
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_3_se_ne_patch->rank);
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_3_nw_ne_patch->rank);
			CHECK(domain_3_ne_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_3_ne_ne_patch->rank);

			CHECK(domain_3_ne_se_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_3_se_nw_patch->rank);
			CHECK(domain_3_ne_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_3_ne_nw_patch->rank);

			CHECK(domain_3_ne_nw_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_3_nw_se_patch->rank);
			CHECK(domain_3_ne_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_3_ne_se_patch->rank);

			CHECK(domain_3_ne_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_3_ne_sw_patch->rank);

			CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_sw_ne_patch->rank);

			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_sw_nw_patch->rank);
			CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_se_nw_patch->rank);

			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_sw_se_patch->rank);
			CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_nw_se_patch->rank);

			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_sw_sw_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_se_sw_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_nw_sw_patch->rank);
			CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_ne_sw_patch->rank);

			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_se_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_se_ne_patch->rank);

			CHECK(domain_2_se_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_se_nw_patch->rank);

			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_sw_se_patch->rank);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_se_se_patch->rank);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_nw_se_patch->rank);
			CHECK(domain_2_se_nw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_ne_se_patch->rank);

			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_se_sw_patch->rank);
			CHECK(domain_2_se_ne_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_ne_sw_patch->rank);

			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_nw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_nw_ne_patch->rank);

			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_sw_nw_patch->rank);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_se_nw_patch->rank);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_nw_nw_patch->rank);
			CHECK(domain_2_nw_se_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_ne_nw_patch->rank);

			CHECK(domain_2_nw_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_nw_se_patch->rank);

			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_nw_sw_patch->rank);
			CHECK(domain_2_nw_ne_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_ne_sw_patch->rank);

			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_sw_ne_patch->rank);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_se_ne_patch->rank);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_nw_ne_patch->rank);
			CHECK(domain_2_ne_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_ne_ne_patch->rank);

			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_se_nw_patch->rank);
			CHECK(domain_2_ne_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_ne_nw_patch->rank);

			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_nw_se_patch->rank);
			CHECK(domain_2_ne_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_ne_se_patch->rank);

			CHECK(domain_2_ne_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_ne_sw_patch->rank);

			CHECK(domain_1_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_1_ne_patch->rank);

			CHECK(domain_1_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_1_nw_patch->rank);

			CHECK(domain_1_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_1_se_patch->rank);

			CHECK(domain_1_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_1_sw_patch->rank);
		}
	}
}
