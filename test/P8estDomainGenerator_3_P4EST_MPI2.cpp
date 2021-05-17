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

#include <ThunderEgg/P8estDomainGenerator.h>

#include <p8est_extended.h>
#include <p8est_mesh.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

namespace
{
std::vector<PatchInfo<3>> GetAllPatchesOnRank0(std::shared_ptr<const Domain<3>> domain)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	std::vector<PatchInfo<3>> all_patches;

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
				for (int i = 0; i < 3; i++) {
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
std::string GetAllPatchesJSONString(const std::vector<PatchInfo<3>> &patches)
{
	nlohmann::json patches_j = patches;
	return patches_j.dump(1);
}
} // namespace
TEST_CASE("P8estDomainGenerator 2x2x2 refined bsw", "[P8estDomainGenerator]")
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	p8est_connectivity_t *conn = p8est_connectivity_new_unitcube();

	p8est_t *p8est = p8est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

	p8est_refine(
	p8est, false,
	[](p8est_t *p8est, p4est_topidx_t witch_tree, p8est_quadrant_t *quadrant) -> int { return 1; },
	nullptr);
	p8est_refine(
	p8est, false,
	[](p8est_t *p8est, p4est_topidx_t witch_tree, p8est_quadrant_t *quadrant) -> int { return quadrant->x == 0 && quadrant->y == 0 && quadrant->z == 0; },
	nullptr);

	p8est_partition(p8est, true, nullptr);

	int    nx              = GENERATE(5, 10);
	int    ny              = GENERATE(5, 10);
	int    nz              = GENERATE(5, 10);
	double scale_x         = GENERATE(0.5, 1.0);
	double scale_y         = GENERATE(0.5, 1.0);
	double scale_z         = GENERATE(0.5, 1.0);
	int    num_ghost_cells = GENERATE(0, 1, 2);

	INFO("nx: " << nx);
	INFO("ny: " << ny);
	INFO("nz: " << nz);
	INFO("scale_x: " << scale_x);
	INFO("scale_y: " << scale_y);
	INFO("scale_z: " << scale_z);
	INFO("num_ghost_cells: " << num_ghost_cells);

	P8estDomainGenerator::BlockMapFunc bmf
	= [&](int block_no, double unit_x, double unit_y, double unit_z, double &x, double &y, double &z) {
		  x = scale_x * unit_x;
		  y = scale_y * unit_y;
		  z = scale_z * unit_z;
	  };

	P8estDomainGenerator dg(p8est, {nx, ny, nz}, num_ghost_cells, bmf);

	auto domain_2 = dg.getFinestDomain();
	auto domain_1 = dg.getCoarserDomain();
	auto domain_0 = dg.getCoarserDomain();

	SECTION("correct number of patches")
	{
		CHECK(domain_2->getNumGlobalPatches() == 15);
		CHECK(domain_1->getNumGlobalPatches() == 8);
		CHECK(domain_0->getNumGlobalPatches() == 1);
	}
	SECTION("patches have correct ns")
	{
		for (auto patch : domain_2->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
			CHECK(patch.ns[2] == nz);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
			CHECK(patch.ns[2] == nz);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.ns[0] == nx);
			CHECK(patch.ns[1] == ny);
			CHECK(patch.ns[2] == nz);
		}
	}

	SECTION("patches have ranks set")
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
	SECTION("patches have num_ghost_cells set")
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

	std::vector<PatchInfo<3>> domain_2_patches = GetAllPatchesOnRank0(domain_2);
	std::vector<PatchInfo<3>> domain_1_patches = GetAllPatchesOnRank0(domain_1);
	std::vector<PatchInfo<3>> domain_0_patches = GetAllPatchesOnRank0(domain_0);

	INFO("nx: " << nx);
	INFO("ny: " << ny);
	INFO("scale_x: " << scale_x);
	INFO("scale_y: " << scale_y);
	INFO("num_ghost_cells: " << num_ghost_cells);

	const PatchInfo<3> *domain_2_bsw_bsw_patch = nullptr;
	const PatchInfo<3> *domain_2_bsw_bse_patch = nullptr;
	const PatchInfo<3> *domain_2_bsw_bnw_patch = nullptr;
	const PatchInfo<3> *domain_2_bsw_bne_patch = nullptr;
	const PatchInfo<3> *domain_2_bsw_tsw_patch = nullptr;
	const PatchInfo<3> *domain_2_bsw_tse_patch = nullptr;
	const PatchInfo<3> *domain_2_bsw_tnw_patch = nullptr;
	const PatchInfo<3> *domain_2_bsw_tne_patch = nullptr;

	const PatchInfo<3> *domain_2_bse_patch = nullptr;
	const PatchInfo<3> *domain_2_bnw_patch = nullptr;
	const PatchInfo<3> *domain_2_bne_patch = nullptr;
	const PatchInfo<3> *domain_2_tsw_patch = nullptr;
	const PatchInfo<3> *domain_2_tse_patch = nullptr;
	const PatchInfo<3> *domain_2_tnw_patch = nullptr;
	const PatchInfo<3> *domain_2_tne_patch = nullptr;

	const PatchInfo<3> *domain_1_bsw_patch = nullptr;
	const PatchInfo<3> *domain_1_bse_patch = nullptr;
	const PatchInfo<3> *domain_1_bnw_patch = nullptr;
	const PatchInfo<3> *domain_1_bne_patch = nullptr;
	const PatchInfo<3> *domain_1_tsw_patch = nullptr;
	const PatchInfo<3> *domain_1_tse_patch = nullptr;
	const PatchInfo<3> *domain_1_tnw_patch = nullptr;
	const PatchInfo<3> *domain_1_tne_patch = nullptr;

	const PatchInfo<3> *domain_0_coarser_patch = nullptr;

	if (rank == 0) {
		for (PatchInfo<3> &patch : domain_2_patches) {
			double x = patch.starts[0];
			double y = patch.starts[1];
			double z = patch.starts[2];

			if (x == Catch::Approx(0) && y == Catch::Approx(0) && z == Catch::Approx(0 * scale_z)) {
				domain_2_bsw_bsw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0) && z == Catch::Approx(0 * scale_z)) {
				domain_2_bsw_bse_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.25 * scale_y) && z == Catch::Approx(0 * scale_z)) {
				domain_2_bsw_bnw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.25 * scale_y) && z == Catch::Approx(0 * scale_z)) {
				domain_2_bsw_bne_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0) && z == Catch::Approx(0.25 * scale_z)) {
				domain_2_bsw_tsw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0) && z == Catch::Approx(0.25 * scale_z)) {
				domain_2_bsw_tse_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.25 * scale_y) && z == Catch::Approx(0.25 * scale_z)) {
				domain_2_bsw_tnw_patch = &patch;
			}
			if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.25 * scale_y) && z == Catch::Approx(0.25 * scale_z)) {
				domain_2_bsw_tne_patch = &patch;
			}

			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0) && z == Catch::Approx(0 * scale_z)) {
				domain_2_bse_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y) && z == Catch::Approx(0 * scale_z)) {
				domain_2_bnw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y) && z == Catch::Approx(0 * scale_z)) {
				domain_2_bne_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0) && z == Catch::Approx(0.5 * scale_z)) {
				domain_2_tsw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0) && z == Catch::Approx(0.5 * scale_z)) {
				domain_2_tse_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y) && z == Catch::Approx(0.5 * scale_z)) {
				domain_2_tnw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y) && z == Catch::Approx(0.5 * scale_z)) {
				domain_2_tne_patch = &patch;
			}
		}

		for (PatchInfo<3> &patch : domain_1_patches) {
			double x = patch.starts[0];
			double y = patch.starts[1];
			double z = patch.starts[2];
			if (x == Catch::Approx(0) && y == Catch::Approx(0) && z == Catch::Approx(0 * scale_z)) {
				domain_1_bsw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0) && z == Catch::Approx(0 * scale_z)) {
				domain_1_bse_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y) && z == Catch::Approx(0 * scale_z)) {
				domain_1_bnw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y) && z == Catch::Approx(0 * scale_z)) {
				domain_1_bne_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0) && z == Catch::Approx(0.5 * scale_z)) {
				domain_1_tsw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0) && z == Catch::Approx(0.5 * scale_z)) {
				domain_1_tse_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y) && z == Catch::Approx(0.5 * scale_z)) {
				domain_1_tnw_patch = &patch;
			}
			if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y) && z == Catch::Approx(0.5 * scale_z)) {
				domain_1_tne_patch = &patch;
			}
		}

		domain_0_coarser_patch = &domain_0_patches[0];

		REQUIRE(domain_2_bsw_bsw_patch != nullptr);
		REQUIRE(domain_2_bsw_bse_patch != nullptr);
		REQUIRE(domain_2_bsw_bnw_patch != nullptr);
		REQUIRE(domain_2_bsw_bne_patch != nullptr);
		REQUIRE(domain_2_bsw_tsw_patch != nullptr);
		REQUIRE(domain_2_bsw_tse_patch != nullptr);
		REQUIRE(domain_2_bsw_tnw_patch != nullptr);
		REQUIRE(domain_2_bsw_tne_patch != nullptr);

		REQUIRE(domain_2_bse_patch != nullptr);
		REQUIRE(domain_2_bnw_patch != nullptr);
		REQUIRE(domain_2_bne_patch != nullptr);
		REQUIRE(domain_2_tsw_patch != nullptr);
		REQUIRE(domain_2_tse_patch != nullptr);
		REQUIRE(domain_2_tnw_patch != nullptr);
		REQUIRE(domain_2_tne_patch != nullptr);

		REQUIRE(domain_1_bsw_patch != nullptr);
		REQUIRE(domain_1_bse_patch != nullptr);
		REQUIRE(domain_1_bnw_patch != nullptr);
		REQUIRE(domain_1_bne_patch != nullptr);
		REQUIRE(domain_1_tsw_patch != nullptr);
		REQUIRE(domain_1_tse_patch != nullptr);
		REQUIRE(domain_1_tnw_patch != nullptr);
		REQUIRE(domain_1_tne_patch != nullptr);

		CHECK(domain_0_coarser_patch->starts[0] == Catch::Approx(0.0));
		CHECK(domain_0_coarser_patch->starts[1] == Catch::Approx(0.0));
		CHECK(domain_0_coarser_patch->starts[2] == Catch::Approx(0.0));
	}

	SECTION("patches have correct spacings")
	{
		if (rank == 0) {
			REQUIRE(domain_2_bsw_bsw_patch->spacings[0] == Catch::Approx(0.25 * scale_x / nx));
			REQUIRE(domain_2_bsw_bsw_patch->spacings[1] == Catch::Approx(0.25 * scale_y / ny));
			REQUIRE(domain_2_bsw_bsw_patch->spacings[2] == Catch::Approx(0.25 * scale_z / nz));
			REQUIRE(domain_2_bsw_bse_patch->spacings[0] == Catch::Approx(0.25 * scale_x / nx));
			REQUIRE(domain_2_bsw_bse_patch->spacings[1] == Catch::Approx(0.25 * scale_y / ny));
			REQUIRE(domain_2_bsw_bse_patch->spacings[2] == Catch::Approx(0.25 * scale_z / nz));
			REQUIRE(domain_2_bsw_bnw_patch->spacings[0] == Catch::Approx(0.25 * scale_x / nx));
			REQUIRE(domain_2_bsw_bnw_patch->spacings[1] == Catch::Approx(0.25 * scale_y / ny));
			REQUIRE(domain_2_bsw_bnw_patch->spacings[2] == Catch::Approx(0.25 * scale_z / nz));
			REQUIRE(domain_2_bsw_bne_patch->spacings[0] == Catch::Approx(0.25 * scale_x / nx));
			REQUIRE(domain_2_bsw_bne_patch->spacings[1] == Catch::Approx(0.25 * scale_y / ny));
			REQUIRE(domain_2_bsw_bne_patch->spacings[2] == Catch::Approx(0.25 * scale_z / nz));
			REQUIRE(domain_2_bsw_tsw_patch->spacings[0] == Catch::Approx(0.25 * scale_x / nx));
			REQUIRE(domain_2_bsw_tsw_patch->spacings[1] == Catch::Approx(0.25 * scale_y / ny));
			REQUIRE(domain_2_bsw_tsw_patch->spacings[2] == Catch::Approx(0.25 * scale_z / nz));
			REQUIRE(domain_2_bsw_tse_patch->spacings[0] == Catch::Approx(0.25 * scale_x / nx));
			REQUIRE(domain_2_bsw_tse_patch->spacings[1] == Catch::Approx(0.25 * scale_y / ny));
			REQUIRE(domain_2_bsw_tse_patch->spacings[2] == Catch::Approx(0.25 * scale_z / nz));
			REQUIRE(domain_2_bsw_tnw_patch->spacings[0] == Catch::Approx(0.25 * scale_x / nx));
			REQUIRE(domain_2_bsw_tnw_patch->spacings[1] == Catch::Approx(0.25 * scale_y / ny));
			REQUIRE(domain_2_bsw_tnw_patch->spacings[2] == Catch::Approx(0.25 * scale_z / nz));
			REQUIRE(domain_2_bsw_tne_patch->spacings[0] == Catch::Approx(0.25 * scale_x / nx));
			REQUIRE(domain_2_bsw_tne_patch->spacings[1] == Catch::Approx(0.25 * scale_y / ny));
			REQUIRE(domain_2_bsw_tne_patch->spacings[2] == Catch::Approx(0.25 * scale_z / nz));

			REQUIRE(domain_2_bse_patch->spacings[0] == Catch::Approx(0.5 * scale_x / nx));
			REQUIRE(domain_2_bse_patch->spacings[1] == Catch::Approx(0.5 * scale_y / ny));
			REQUIRE(domain_2_bse_patch->spacings[2] == Catch::Approx(0.5 * scale_z / nz));
			REQUIRE(domain_2_bnw_patch->spacings[0] == Catch::Approx(0.5 * scale_x / nx));
			REQUIRE(domain_2_bnw_patch->spacings[1] == Catch::Approx(0.5 * scale_y / ny));
			REQUIRE(domain_2_bnw_patch->spacings[2] == Catch::Approx(0.5 * scale_z / nz));
			REQUIRE(domain_2_bne_patch->spacings[0] == Catch::Approx(0.5 * scale_x / nx));
			REQUIRE(domain_2_bne_patch->spacings[1] == Catch::Approx(0.5 * scale_y / ny));
			REQUIRE(domain_2_bne_patch->spacings[2] == Catch::Approx(0.5 * scale_z / nz));
			REQUIRE(domain_2_tsw_patch->spacings[0] == Catch::Approx(0.5 * scale_x / nx));
			REQUIRE(domain_2_tsw_patch->spacings[1] == Catch::Approx(0.5 * scale_y / ny));
			REQUIRE(domain_2_tsw_patch->spacings[2] == Catch::Approx(0.5 * scale_z / nz));
			REQUIRE(domain_2_tse_patch->spacings[0] == Catch::Approx(0.5 * scale_x / nx));
			REQUIRE(domain_2_tse_patch->spacings[1] == Catch::Approx(0.5 * scale_y / ny));
			REQUIRE(domain_2_tse_patch->spacings[2] == Catch::Approx(0.5 * scale_z / nz));
			REQUIRE(domain_2_tnw_patch->spacings[0] == Catch::Approx(0.5 * scale_x / nx));
			REQUIRE(domain_2_tnw_patch->spacings[1] == Catch::Approx(0.5 * scale_y / ny));
			REQUIRE(domain_2_tnw_patch->spacings[2] == Catch::Approx(0.5 * scale_z / nz));
			REQUIRE(domain_2_tne_patch->spacings[0] == Catch::Approx(0.5 * scale_x / nx));
			REQUIRE(domain_2_tne_patch->spacings[1] == Catch::Approx(0.5 * scale_y / ny));
			REQUIRE(domain_2_tne_patch->spacings[2] == Catch::Approx(0.5 * scale_z / nz));
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
			CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
			CHECK(patch.spacings[2] == Catch::Approx(scale_z * 0.5 / nz));
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.spacings[0] == Catch::Approx(scale_x * 1.0 / nx));
			CHECK(patch.spacings[1] == Catch::Approx(scale_y * 1.0 / ny));
			CHECK(patch.spacings[2] == Catch::Approx(scale_z * 1.0 / nz));
		}
	}

	SECTION("patches have refine_level set")
	{
		if (rank == 0) {
			REQUIRE(domain_2_bsw_bsw_patch->refine_level == 2);
			REQUIRE(domain_2_bsw_bse_patch->refine_level == 2);
			REQUIRE(domain_2_bsw_bnw_patch->refine_level == 2);
			REQUIRE(domain_2_bsw_bne_patch->refine_level == 2);
			REQUIRE(domain_2_bsw_tsw_patch->refine_level == 2);
			REQUIRE(domain_2_bsw_tse_patch->refine_level == 2);
			REQUIRE(domain_2_bsw_tnw_patch->refine_level == 2);
			REQUIRE(domain_2_bsw_tne_patch->refine_level == 2);

			REQUIRE(domain_2_bse_patch->refine_level == 1);
			REQUIRE(domain_2_bnw_patch->refine_level == 1);
			REQUIRE(domain_2_bne_patch->refine_level == 1);
			REQUIRE(domain_2_tsw_patch->refine_level == 1);
			REQUIRE(domain_2_tse_patch->refine_level == 1);
			REQUIRE(domain_2_tnw_patch->refine_level == 1);
			REQUIRE(domain_2_tne_patch->refine_level == 1);
		}

		for (auto patch : domain_1->getPatchInfoVector()) {
			CHECK(patch.refine_level == 1);
		}

		for (auto patch : domain_0->getPatchInfoVector()) {
			CHECK(patch.refine_level == 0);
		}
	}

	SECTION("parent ids are set correctly")
	{
		if (rank == 0) {
			CHECK(domain_2_bsw_bsw_patch->parent_id == domain_1_bsw_patch->id);
			CHECK(domain_2_bsw_bse_patch->parent_id == domain_1_bsw_patch->id);
			CHECK(domain_2_bsw_bnw_patch->parent_id == domain_1_bsw_patch->id);
			CHECK(domain_2_bsw_bne_patch->parent_id == domain_1_bsw_patch->id);
			CHECK(domain_2_bsw_tsw_patch->parent_id == domain_1_bsw_patch->id);
			CHECK(domain_2_bsw_tse_patch->parent_id == domain_1_bsw_patch->id);
			CHECK(domain_2_bsw_tnw_patch->parent_id == domain_1_bsw_patch->id);
			CHECK(domain_2_bsw_tne_patch->parent_id == domain_1_bsw_patch->id);

			CHECK(domain_2_bse_patch->parent_id == domain_1_bse_patch->id);
			CHECK(domain_2_bnw_patch->parent_id == domain_1_bnw_patch->id);
			CHECK(domain_2_bne_patch->parent_id == domain_1_bne_patch->id);
			CHECK(domain_2_tsw_patch->parent_id == domain_1_tsw_patch->id);
			CHECK(domain_2_tse_patch->parent_id == domain_1_tse_patch->id);
			CHECK(domain_2_tnw_patch->parent_id == domain_1_tnw_patch->id);
			CHECK(domain_2_tne_patch->parent_id == domain_1_tne_patch->id);

			CHECK(domain_1_bsw_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_bse_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_bnw_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_bne_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_tsw_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_tse_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_tnw_patch->parent_id == domain_0_coarser_patch->id);
			CHECK(domain_1_tne_patch->parent_id == domain_0_coarser_patch->id);

			CHECK(domain_0_coarser_patch->parent_id == -1);
		}
	}
	SECTION("child ids are set correctly")
	{
		if (rank == 0) {
			for (int i = 0; i < 8; i++) {
				INFO("i: " << i);
				CHECK(domain_2_bsw_bsw_patch->child_ids[i] == -1);
				CHECK(domain_2_bsw_bse_patch->child_ids[i] == -1);
				CHECK(domain_2_bsw_bnw_patch->child_ids[i] == -1);
				CHECK(domain_2_bsw_bne_patch->child_ids[i] == -1);
				CHECK(domain_2_bsw_tsw_patch->child_ids[i] == -1);
				CHECK(domain_2_bsw_tse_patch->child_ids[i] == -1);
				CHECK(domain_2_bsw_tnw_patch->child_ids[i] == -1);
				CHECK(domain_2_bsw_tne_patch->child_ids[i] == -1);

				CHECK(domain_2_bse_patch->child_ids[i] == -1);
				CHECK(domain_2_bnw_patch->child_ids[i] == -1);
				CHECK(domain_2_bne_patch->child_ids[i] == -1);
				CHECK(domain_2_tsw_patch->child_ids[i] == -1);
				CHECK(domain_2_tse_patch->child_ids[i] == -1);
				CHECK(domain_2_tnw_patch->child_ids[i] == -1);
				CHECK(domain_2_tne_patch->child_ids[i] == -1);
			}

			CHECK(domain_1_bsw_patch->child_ids[0] == domain_2_bsw_bsw_patch->id);
			CHECK(domain_1_bsw_patch->child_ids[1] == domain_2_bsw_bse_patch->id);
			CHECK(domain_1_bsw_patch->child_ids[2] == domain_2_bsw_bnw_patch->id);
			CHECK(domain_1_bsw_patch->child_ids[3] == domain_2_bsw_bne_patch->id);
			CHECK(domain_1_bsw_patch->child_ids[4] == domain_2_bsw_tsw_patch->id);
			CHECK(domain_1_bsw_patch->child_ids[5] == domain_2_bsw_tse_patch->id);
			CHECK(domain_1_bsw_patch->child_ids[6] == domain_2_bsw_tnw_patch->id);
			CHECK(domain_1_bsw_patch->child_ids[7] == domain_2_bsw_tne_patch->id);

			CHECK(domain_1_bse_patch->child_ids[0] == domain_2_bse_patch->id);
			CHECK(domain_1_bnw_patch->child_ids[0] == domain_2_bnw_patch->id);
			CHECK(domain_1_bne_patch->child_ids[0] == domain_2_bne_patch->id);
			CHECK(domain_1_tsw_patch->child_ids[0] == domain_2_tsw_patch->id);
			CHECK(domain_1_tse_patch->child_ids[0] == domain_2_tse_patch->id);
			CHECK(domain_1_tnw_patch->child_ids[0] == domain_2_tnw_patch->id);
			CHECK(domain_1_tne_patch->child_ids[0] == domain_2_tne_patch->id);

			for (int i = 1; i < 8; i++) {
				INFO("i: " << i);
				CHECK(domain_1_bse_patch->child_ids[i] == -1);
				CHECK(domain_1_bnw_patch->child_ids[i] == -1);
				CHECK(domain_1_bne_patch->child_ids[i] == -1);
				CHECK(domain_1_tsw_patch->child_ids[i] == -1);
				CHECK(domain_1_tse_patch->child_ids[i] == -1);
				CHECK(domain_1_tnw_patch->child_ids[i] == -1);
				CHECK(domain_1_tne_patch->child_ids[i] == -1);
			}

			CHECK(domain_0_coarser_patch->child_ids[0] == domain_1_bsw_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[1] == domain_1_bse_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[2] == domain_1_bnw_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[3] == domain_1_bne_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[4] == domain_1_tsw_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[5] == domain_1_tse_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[6] == domain_1_tnw_patch->id);
			CHECK(domain_0_coarser_patch->child_ids[7] == domain_1_tne_patch->id);
		}
	}
	SECTION("orth on parent is set correctly")
	{
		if (rank == 0) {
			CHECK(domain_2_bsw_bsw_patch->orth_on_parent == Orthant<3>::bsw());
			CHECK(domain_2_bsw_bse_patch->orth_on_parent == Orthant<3>::bse());
			CHECK(domain_2_bsw_bnw_patch->orth_on_parent == Orthant<3>::bnw());
			CHECK(domain_2_bsw_bne_patch->orth_on_parent == Orthant<3>::bne());
			CHECK(domain_2_bsw_tsw_patch->orth_on_parent == Orthant<3>::tsw());
			CHECK(domain_2_bsw_tse_patch->orth_on_parent == Orthant<3>::tse());
			CHECK(domain_2_bsw_tnw_patch->orth_on_parent == Orthant<3>::tnw());
			CHECK(domain_2_bsw_tne_patch->orth_on_parent == Orthant<3>::tne());

			CHECK(domain_2_bse_patch->orth_on_parent == Orthant<3>::null());
			CHECK(domain_2_bnw_patch->orth_on_parent == Orthant<3>::null());
			CHECK(domain_2_bne_patch->orth_on_parent == Orthant<3>::null());
			CHECK(domain_2_tsw_patch->orth_on_parent == Orthant<3>::null());
			CHECK(domain_2_tse_patch->orth_on_parent == Orthant<3>::null());
			CHECK(domain_2_tnw_patch->orth_on_parent == Orthant<3>::null());
			CHECK(domain_2_tne_patch->orth_on_parent == Orthant<3>::null());

			CHECK(domain_1_bsw_patch->orth_on_parent == Orthant<3>::bsw());
			CHECK(domain_1_bse_patch->orth_on_parent == Orthant<3>::bse());
			CHECK(domain_1_bnw_patch->orth_on_parent == Orthant<3>::bnw());
			CHECK(domain_1_bne_patch->orth_on_parent == Orthant<3>::bne());
			CHECK(domain_1_tsw_patch->orth_on_parent == Orthant<3>::tsw());
			CHECK(domain_1_tse_patch->orth_on_parent == Orthant<3>::tse());
			CHECK(domain_1_tnw_patch->orth_on_parent == Orthant<3>::tnw());
			CHECK(domain_1_tne_patch->orth_on_parent == Orthant<3>::tne());

			CHECK(domain_0_coarser_patch->orth_on_parent == Orthant<3>::null());
		}
	}
	SECTION("parent ranks are set correctly")
	{
		if (rank == 0) {
			CHECK(domain_2_bsw_bsw_patch->parent_rank == domain_1_bsw_patch->rank);
			CHECK(domain_2_bsw_bse_patch->parent_rank == domain_1_bsw_patch->rank);
			CHECK(domain_2_bsw_bnw_patch->parent_rank == domain_1_bsw_patch->rank);
			CHECK(domain_2_bsw_bne_patch->parent_rank == domain_1_bsw_patch->rank);
			CHECK(domain_2_bsw_tsw_patch->parent_rank == domain_1_bsw_patch->rank);
			CHECK(domain_2_bsw_tse_patch->parent_rank == domain_1_bsw_patch->rank);
			CHECK(domain_2_bsw_tnw_patch->parent_rank == domain_1_bsw_patch->rank);
			CHECK(domain_2_bsw_tne_patch->parent_rank == domain_1_bsw_patch->rank);

			CHECK(domain_2_bse_patch->parent_rank == domain_1_bse_patch->rank);
			CHECK(domain_2_bnw_patch->parent_rank == domain_1_bnw_patch->rank);
			CHECK(domain_2_bne_patch->parent_rank == domain_1_bne_patch->rank);
			CHECK(domain_2_tsw_patch->parent_rank == domain_1_tsw_patch->rank);
			CHECK(domain_2_tse_patch->parent_rank == domain_1_tse_patch->rank);
			CHECK(domain_2_tnw_patch->parent_rank == domain_1_tnw_patch->rank);
			CHECK(domain_2_tne_patch->parent_rank == domain_1_tne_patch->rank);

			CHECK(domain_1_bsw_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_bse_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_bnw_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_bne_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_tsw_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_tse_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_tnw_patch->parent_rank == domain_0_coarser_patch->rank);
			CHECK(domain_1_tne_patch->parent_rank == domain_0_coarser_patch->rank);

			CHECK(domain_0_coarser_patch->parent_rank == -1);
		}
	}
	SECTION("child ranks are set correctly")
	{
		if (rank == 0) {
			for (int i = 0; i < 8; i++) {
				INFO("i: " << i);
				CHECK(domain_2_bsw_bsw_patch->child_ranks[i] == -1);
				CHECK(domain_2_bsw_bse_patch->child_ranks[i] == -1);
				CHECK(domain_2_bsw_bnw_patch->child_ranks[i] == -1);
				CHECK(domain_2_bsw_bne_patch->child_ranks[i] == -1);
				CHECK(domain_2_bsw_tsw_patch->child_ranks[i] == -1);
				CHECK(domain_2_bsw_tse_patch->child_ranks[i] == -1);
				CHECK(domain_2_bsw_tnw_patch->child_ranks[i] == -1);
				CHECK(domain_2_bsw_tne_patch->child_ranks[i] == -1);

				CHECK(domain_2_bse_patch->child_ranks[i] == -1);
				CHECK(domain_2_bnw_patch->child_ranks[i] == -1);
				CHECK(domain_2_bne_patch->child_ranks[i] == -1);
				CHECK(domain_2_tsw_patch->child_ranks[i] == -1);
				CHECK(domain_2_tse_patch->child_ranks[i] == -1);
				CHECK(domain_2_tnw_patch->child_ranks[i] == -1);
				CHECK(domain_2_tne_patch->child_ranks[i] == -1);
			}

			CHECK(domain_1_bsw_patch->child_ranks[0] == domain_2_bsw_bsw_patch->rank);
			CHECK(domain_1_bsw_patch->child_ranks[1] == domain_2_bsw_bse_patch->rank);
			CHECK(domain_1_bsw_patch->child_ranks[2] == domain_2_bsw_bnw_patch->rank);
			CHECK(domain_1_bsw_patch->child_ranks[3] == domain_2_bsw_bne_patch->rank);
			CHECK(domain_1_bsw_patch->child_ranks[4] == domain_2_bsw_tsw_patch->rank);
			CHECK(domain_1_bsw_patch->child_ranks[5] == domain_2_bsw_tse_patch->rank);
			CHECK(domain_1_bsw_patch->child_ranks[6] == domain_2_bsw_tnw_patch->rank);
			CHECK(domain_1_bsw_patch->child_ranks[7] == domain_2_bsw_tne_patch->rank);

			CHECK(domain_1_bse_patch->child_ranks[0] == domain_2_bse_patch->rank);
			CHECK(domain_1_bnw_patch->child_ranks[0] == domain_2_bnw_patch->rank);
			CHECK(domain_1_bne_patch->child_ranks[0] == domain_2_bne_patch->rank);
			CHECK(domain_1_tsw_patch->child_ranks[0] == domain_2_tsw_patch->rank);
			CHECK(domain_1_tse_patch->child_ranks[0] == domain_2_tse_patch->rank);
			CHECK(domain_1_tnw_patch->child_ranks[0] == domain_2_tnw_patch->rank);
			CHECK(domain_1_tne_patch->child_ranks[0] == domain_2_tne_patch->rank);

			for (int i = 1; i < 8; i++) {
				INFO("i: " << i);
				CHECK(domain_1_bse_patch->child_ranks[i] == -1);
				CHECK(domain_1_bnw_patch->child_ranks[i] == -1);
				CHECK(domain_1_bne_patch->child_ranks[i] == -1);
				CHECK(domain_1_tsw_patch->child_ranks[i] == -1);
				CHECK(domain_1_tse_patch->child_ranks[i] == -1);
				CHECK(domain_1_tnw_patch->child_ranks[i] == -1);
				CHECK(domain_1_tne_patch->child_ranks[i] == -1);
			}

			CHECK(domain_0_coarser_patch->child_ranks[0] == domain_1_bsw_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[1] == domain_1_bse_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[2] == domain_1_bnw_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[3] == domain_1_bne_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[4] == domain_1_tsw_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[5] == domain_1_tse_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[6] == domain_1_tnw_patch->rank);
			CHECK(domain_0_coarser_patch->child_ranks[7] == domain_1_tne_patch->rank);
		}
	}
	SECTION("correct sides have nbr_infos")
	{
		if (rank == 0) {
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Side<3>::west()) == false);
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Side<3>::south()) == false);
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Side<3>::bottom()) == false);
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_2_bsw_bse_patch->hasNbr(Side<3>::west()) == true);
			CHECK(domain_2_bsw_bse_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_2_bsw_bse_patch->hasNbr(Side<3>::south()) == false);
			CHECK(domain_2_bsw_bse_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_2_bsw_bse_patch->hasNbr(Side<3>::bottom()) == false);
			CHECK(domain_2_bsw_bse_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_2_bsw_bnw_patch->hasNbr(Side<3>::west()) == false);
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Side<3>::south()) == true);
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Side<3>::bottom()) == false);
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_2_bsw_bne_patch->hasNbr(Side<3>::west()) == true);
			CHECK(domain_2_bsw_bne_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_2_bsw_bne_patch->hasNbr(Side<3>::south()) == true);
			CHECK(domain_2_bsw_bne_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_2_bsw_bne_patch->hasNbr(Side<3>::bottom()) == false);
			CHECK(domain_2_bsw_bne_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_2_bsw_tsw_patch->hasNbr(Side<3>::west()) == false);
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Side<3>::south()) == false);
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Side<3>::bottom()) == true);
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_2_bsw_tse_patch->hasNbr(Side<3>::west()) == true);
			CHECK(domain_2_bsw_tse_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_2_bsw_tse_patch->hasNbr(Side<3>::south()) == false);
			CHECK(domain_2_bsw_tse_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_2_bsw_tse_patch->hasNbr(Side<3>::bottom()) == true);
			CHECK(domain_2_bsw_tse_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_2_bsw_tnw_patch->hasNbr(Side<3>::west()) == false);
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Side<3>::south()) == true);
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Side<3>::bottom()) == true);
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_2_bsw_tne_patch->hasNbr(Side<3>::west()) == true);
			CHECK(domain_2_bsw_tne_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_2_bsw_tne_patch->hasNbr(Side<3>::south()) == true);
			CHECK(domain_2_bsw_tne_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_2_bsw_tne_patch->hasNbr(Side<3>::bottom()) == true);
			CHECK(domain_2_bsw_tne_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_2_bse_patch->hasNbr(Side<3>::west()) == true);
			CHECK(domain_2_bse_patch->hasNbr(Side<3>::east()) == false);
			CHECK(domain_2_bse_patch->hasNbr(Side<3>::south()) == false);
			CHECK(domain_2_bse_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_2_bse_patch->hasNbr(Side<3>::bottom()) == false);
			CHECK(domain_2_bse_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_2_bnw_patch->hasNbr(Side<3>::west()) == false);
			CHECK(domain_2_bnw_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_2_bnw_patch->hasNbr(Side<3>::south()) == true);
			CHECK(domain_2_bnw_patch->hasNbr(Side<3>::north()) == false);
			CHECK(domain_2_bnw_patch->hasNbr(Side<3>::bottom()) == false);
			CHECK(domain_2_bnw_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_2_bne_patch->hasNbr(Side<3>::west()) == true);
			CHECK(domain_2_bne_patch->hasNbr(Side<3>::east()) == false);
			CHECK(domain_2_bne_patch->hasNbr(Side<3>::south()) == true);
			CHECK(domain_2_bne_patch->hasNbr(Side<3>::north()) == false);
			CHECK(domain_2_bne_patch->hasNbr(Side<3>::bottom()) == false);
			CHECK(domain_2_bne_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_2_tsw_patch->hasNbr(Side<3>::west()) == false);
			CHECK(domain_2_tsw_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_2_tsw_patch->hasNbr(Side<3>::south()) == false);
			CHECK(domain_2_tsw_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_2_tsw_patch->hasNbr(Side<3>::bottom()) == true);
			CHECK(domain_2_tsw_patch->hasNbr(Side<3>::top()) == false);

			CHECK(domain_2_tse_patch->hasNbr(Side<3>::west()) == true);
			CHECK(domain_2_tse_patch->hasNbr(Side<3>::east()) == false);
			CHECK(domain_2_tse_patch->hasNbr(Side<3>::south()) == false);
			CHECK(domain_2_tse_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_2_tse_patch->hasNbr(Side<3>::bottom()) == true);
			CHECK(domain_2_tse_patch->hasNbr(Side<3>::top()) == false);

			CHECK(domain_2_tnw_patch->hasNbr(Side<3>::west()) == false);
			CHECK(domain_2_tnw_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_2_tnw_patch->hasNbr(Side<3>::south()) == true);
			CHECK(domain_2_tnw_patch->hasNbr(Side<3>::north()) == false);
			CHECK(domain_2_tnw_patch->hasNbr(Side<3>::bottom()) == true);
			CHECK(domain_2_tnw_patch->hasNbr(Side<3>::top()) == false);

			CHECK(domain_2_tne_patch->hasNbr(Side<3>::west()) == true);
			CHECK(domain_2_tne_patch->hasNbr(Side<3>::east()) == false);
			CHECK(domain_2_tne_patch->hasNbr(Side<3>::south()) == true);
			CHECK(domain_2_tne_patch->hasNbr(Side<3>::north()) == false);
			CHECK(domain_2_tne_patch->hasNbr(Side<3>::bottom()) == true);
			CHECK(domain_2_tne_patch->hasNbr(Side<3>::top()) == false);

			CHECK(domain_1_bsw_patch->hasNbr(Side<3>::west()) == false);
			CHECK(domain_1_bsw_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_1_bsw_patch->hasNbr(Side<3>::south()) == false);
			CHECK(domain_1_bsw_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_1_bsw_patch->hasNbr(Side<3>::bottom()) == false);
			CHECK(domain_1_bsw_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_1_bse_patch->hasNbr(Side<3>::west()) == true);
			CHECK(domain_1_bse_patch->hasNbr(Side<3>::east()) == false);
			CHECK(domain_1_bse_patch->hasNbr(Side<3>::south()) == false);
			CHECK(domain_1_bse_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_1_bse_patch->hasNbr(Side<3>::bottom()) == false);
			CHECK(domain_1_bse_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_1_bnw_patch->hasNbr(Side<3>::west()) == false);
			CHECK(domain_1_bnw_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_1_bnw_patch->hasNbr(Side<3>::south()) == true);
			CHECK(domain_1_bnw_patch->hasNbr(Side<3>::north()) == false);
			CHECK(domain_1_bnw_patch->hasNbr(Side<3>::bottom()) == false);
			CHECK(domain_1_bnw_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_1_bne_patch->hasNbr(Side<3>::west()) == true);
			CHECK(domain_1_bne_patch->hasNbr(Side<3>::east()) == false);
			CHECK(domain_1_bne_patch->hasNbr(Side<3>::south()) == true);
			CHECK(domain_1_bne_patch->hasNbr(Side<3>::north()) == false);
			CHECK(domain_1_bne_patch->hasNbr(Side<3>::bottom()) == false);
			CHECK(domain_1_bne_patch->hasNbr(Side<3>::top()) == true);

			CHECK(domain_1_tsw_patch->hasNbr(Side<3>::west()) == false);
			CHECK(domain_1_tsw_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_1_tsw_patch->hasNbr(Side<3>::south()) == false);
			CHECK(domain_1_tsw_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_1_tsw_patch->hasNbr(Side<3>::bottom()) == true);
			CHECK(domain_1_tsw_patch->hasNbr(Side<3>::top()) == false);

			CHECK(domain_1_tse_patch->hasNbr(Side<3>::west()) == true);
			CHECK(domain_1_tse_patch->hasNbr(Side<3>::east()) == false);
			CHECK(domain_1_tse_patch->hasNbr(Side<3>::south()) == false);
			CHECK(domain_1_tse_patch->hasNbr(Side<3>::north()) == true);
			CHECK(domain_1_tse_patch->hasNbr(Side<3>::bottom()) == true);
			CHECK(domain_1_tse_patch->hasNbr(Side<3>::top()) == false);

			CHECK(domain_1_tnw_patch->hasNbr(Side<3>::west()) == false);
			CHECK(domain_1_tnw_patch->hasNbr(Side<3>::east()) == true);
			CHECK(domain_1_tnw_patch->hasNbr(Side<3>::south()) == true);
			CHECK(domain_1_tnw_patch->hasNbr(Side<3>::north()) == false);
			CHECK(domain_1_tnw_patch->hasNbr(Side<3>::bottom()) == true);
			CHECK(domain_1_tnw_patch->hasNbr(Side<3>::top()) == false);

			CHECK(domain_1_tne_patch->hasNbr(Side<3>::west()) == true);
			CHECK(domain_1_tne_patch->hasNbr(Side<3>::east()) == false);
			CHECK(domain_1_tne_patch->hasNbr(Side<3>::south()) == true);
			CHECK(domain_1_tne_patch->hasNbr(Side<3>::north()) == false);
			CHECK(domain_1_tne_patch->hasNbr(Side<3>::bottom()) == true);
			CHECK(domain_1_tne_patch->hasNbr(Side<3>::top()) == false);

			CHECK(domain_0_coarser_patch->hasNbr(Side<3>::west()) == false);
			CHECK(domain_0_coarser_patch->hasNbr(Side<3>::east()) == false);
			CHECK(domain_0_coarser_patch->hasNbr(Side<3>::south()) == false);
			CHECK(domain_0_coarser_patch->hasNbr(Side<3>::north()) == false);
			CHECK(domain_0_coarser_patch->hasNbr(Side<3>::bottom()) == false);
			CHECK(domain_0_coarser_patch->hasNbr(Side<3>::top()) == false);
		}
	}
	SECTION("correct sides have nbr_infos")
	{
		if (rank == 0) {
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_2_bsw_bse_patch->id);
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Side<3>::north()).id == domain_2_bsw_bnw_patch->id);
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Side<3>::top()).id == domain_2_bsw_tsw_patch->id);

			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Side<3>::west()).id == domain_2_bsw_bsw_patch->id);
			CHECK(domain_2_bsw_bse_patch->getCoarseNbrInfo(Side<3>::east()).id == domain_2_bse_patch->id);
			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Side<3>::north()).id == domain_2_bsw_bne_patch->id);
			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Side<3>::top()).id == domain_2_bsw_tse_patch->id);

			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_2_bsw_bne_patch->id);
			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Side<3>::south()).id == domain_2_bsw_bsw_patch->id);
			CHECK(domain_2_bsw_bnw_patch->getCoarseNbrInfo(Side<3>::north()).id == domain_2_bnw_patch->id);
			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Side<3>::top()).id == domain_2_bsw_tnw_patch->id);

			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Side<3>::west()).id == domain_2_bsw_bnw_patch->id);
			CHECK(domain_2_bsw_bne_patch->getCoarseNbrInfo(Side<3>::east()).id == domain_2_bse_patch->id);
			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Side<3>::south()).id == domain_2_bsw_bse_patch->id);
			CHECK(domain_2_bsw_bne_patch->getCoarseNbrInfo(Side<3>::north()).id == domain_2_bnw_patch->id);
			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Side<3>::top()).id == domain_2_bsw_tne_patch->id);

			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_2_bsw_tse_patch->id);
			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Side<3>::north()).id == domain_2_bsw_tnw_patch->id);
			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_2_bsw_bsw_patch->id);
			CHECK(domain_2_bsw_tsw_patch->getCoarseNbrInfo(Side<3>::top()).id == domain_2_tsw_patch->id);

			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Side<3>::west()).id == domain_2_bsw_tsw_patch->id);
			CHECK(domain_2_bsw_tse_patch->getCoarseNbrInfo(Side<3>::east()).id == domain_2_bse_patch->id);
			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Side<3>::north()).id == domain_2_bsw_tne_patch->id);
			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_2_bsw_bse_patch->id);
			CHECK(domain_2_bsw_tse_patch->getCoarseNbrInfo(Side<3>::top()).id == domain_2_tsw_patch->id);

			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_2_bsw_tne_patch->id);
			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Side<3>::south()).id == domain_2_bsw_tsw_patch->id);
			CHECK(domain_2_bsw_tnw_patch->getCoarseNbrInfo(Side<3>::north()).id == domain_2_bnw_patch->id);
			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_2_bsw_bnw_patch->id);
			CHECK(domain_2_bsw_tnw_patch->getCoarseNbrInfo(Side<3>::top()).id == domain_2_tsw_patch->id);

			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Side<3>::west()).id == domain_2_bsw_tnw_patch->id);
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Side<3>::east()).id == domain_2_bse_patch->id);
			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Side<3>::south()).id == domain_2_bsw_tse_patch->id);
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Side<3>::north()).id == domain_2_bnw_patch->id);
			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_2_bsw_bne_patch->id);
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Side<3>::top()).id == domain_2_tsw_patch->id);

			CHECK(domain_2_bse_patch->getFineNbrInfo(Side<3>::west()).ids[0] == domain_2_bsw_bse_patch->id);
			CHECK(domain_2_bse_patch->getFineNbrInfo(Side<3>::west()).ids[1] == domain_2_bsw_bne_patch->id);
			CHECK(domain_2_bse_patch->getFineNbrInfo(Side<3>::west()).ids[2] == domain_2_bsw_tse_patch->id);
			CHECK(domain_2_bse_patch->getFineNbrInfo(Side<3>::west()).ids[3] == domain_2_bsw_tne_patch->id);
			CHECK(domain_2_bse_patch->getNormalNbrInfo(Side<3>::north()).id == domain_2_bne_patch->id);
			CHECK(domain_2_bse_patch->getNormalNbrInfo(Side<3>::top()).id == domain_2_tse_patch->id);

			CHECK(domain_2_bnw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_2_bne_patch->id);
			CHECK(domain_2_bnw_patch->getFineNbrInfo(Side<3>::south()).ids[0] == domain_2_bsw_bnw_patch->id);
			CHECK(domain_2_bnw_patch->getFineNbrInfo(Side<3>::south()).ids[1] == domain_2_bsw_bne_patch->id);
			CHECK(domain_2_bnw_patch->getFineNbrInfo(Side<3>::south()).ids[2] == domain_2_bsw_tnw_patch->id);
			CHECK(domain_2_bnw_patch->getFineNbrInfo(Side<3>::south()).ids[3] == domain_2_bsw_tne_patch->id);
			CHECK(domain_2_bnw_patch->getNormalNbrInfo(Side<3>::top()).id == domain_2_tnw_patch->id);

			CHECK(domain_2_bne_patch->getNormalNbrInfo(Side<3>::west()).id == domain_2_bnw_patch->id);
			CHECK(domain_2_bne_patch->getNormalNbrInfo(Side<3>::south()).id == domain_2_bse_patch->id);
			CHECK(domain_2_bne_patch->getNormalNbrInfo(Side<3>::top()).id == domain_2_tne_patch->id);

			CHECK(domain_2_tsw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_2_tse_patch->id);
			CHECK(domain_2_tsw_patch->getNormalNbrInfo(Side<3>::north()).id == domain_2_tnw_patch->id);
			CHECK(domain_2_tsw_patch->getFineNbrInfo(Side<3>::bottom()).ids[0] == domain_2_bsw_tsw_patch->id);
			CHECK(domain_2_tsw_patch->getFineNbrInfo(Side<3>::bottom()).ids[1] == domain_2_bsw_tse_patch->id);
			CHECK(domain_2_tsw_patch->getFineNbrInfo(Side<3>::bottom()).ids[2] == domain_2_bsw_tnw_patch->id);
			CHECK(domain_2_tsw_patch->getFineNbrInfo(Side<3>::bottom()).ids[3] == domain_2_bsw_tne_patch->id);

			CHECK(domain_2_tse_patch->getNormalNbrInfo(Side<3>::west()).id == domain_2_tsw_patch->id);
			CHECK(domain_2_tse_patch->getNormalNbrInfo(Side<3>::north()).id == domain_2_tne_patch->id);
			CHECK(domain_2_tse_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_2_bse_patch->id);

			CHECK(domain_2_tnw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_2_tne_patch->id);
			CHECK(domain_2_tnw_patch->getNormalNbrInfo(Side<3>::south()).id == domain_2_tsw_patch->id);
			CHECK(domain_2_tnw_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_2_bnw_patch->id);

			CHECK(domain_2_tne_patch->getNormalNbrInfo(Side<3>::west()).id == domain_2_tnw_patch->id);
			CHECK(domain_2_tne_patch->getNormalNbrInfo(Side<3>::south()).id == domain_2_tse_patch->id);
			CHECK(domain_2_tne_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_2_bne_patch->id);

			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_1_bse_patch->id);
			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Side<3>::north()).id == domain_1_bnw_patch->id);
			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Side<3>::top()).id == domain_1_tsw_patch->id);

			CHECK(domain_1_bse_patch->getNormalNbrInfo(Side<3>::west()).id == domain_1_bsw_patch->id);
			CHECK(domain_1_bse_patch->getNormalNbrInfo(Side<3>::north()).id == domain_1_bne_patch->id);
			CHECK(domain_1_bse_patch->getNormalNbrInfo(Side<3>::top()).id == domain_1_tse_patch->id);

			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_1_bne_patch->id);
			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Side<3>::south()).id == domain_1_bsw_patch->id);
			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Side<3>::top()).id == domain_1_tnw_patch->id);

			CHECK(domain_1_bne_patch->getNormalNbrInfo(Side<3>::west()).id == domain_1_bnw_patch->id);
			CHECK(domain_1_bne_patch->getNormalNbrInfo(Side<3>::south()).id == domain_1_bse_patch->id);
			CHECK(domain_1_bne_patch->getNormalNbrInfo(Side<3>::top()).id == domain_1_tne_patch->id);

			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_1_tse_patch->id);
			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Side<3>::north()).id == domain_1_tnw_patch->id);
			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_1_bsw_patch->id);

			CHECK(domain_1_tse_patch->getNormalNbrInfo(Side<3>::west()).id == domain_1_tsw_patch->id);
			CHECK(domain_1_tse_patch->getNormalNbrInfo(Side<3>::north()).id == domain_1_tne_patch->id);
			CHECK(domain_1_tse_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_1_bse_patch->id);

			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_1_tne_patch->id);
			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Side<3>::south()).id == domain_1_tsw_patch->id);
			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_1_bnw_patch->id);

			CHECK(domain_1_tne_patch->getNormalNbrInfo(Side<3>::west()).id == domain_1_tnw_patch->id);
			CHECK(domain_1_tne_patch->getNormalNbrInfo(Side<3>::south()).id == domain_1_tse_patch->id);
			CHECK(domain_1_tne_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_1_bne_patch->id);
		}
	}

	SECTION("nbr_info ranks are correct")
	{
		if (rank == 0) {
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_2_bsw_bse_patch->rank);
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_2_bsw_bnw_patch->rank);
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_2_bsw_tsw_patch->rank);

			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_2_bsw_bsw_patch->rank);
			CHECK(domain_2_bsw_bse_patch->getCoarseNbrInfo(Side<3>::east()).rank == domain_2_bse_patch->rank);
			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_2_bsw_bne_patch->rank);
			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_2_bsw_tse_patch->rank);

			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_2_bsw_bne_patch->rank);
			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_2_bsw_bsw_patch->rank);
			CHECK(domain_2_bsw_bnw_patch->getCoarseNbrInfo(Side<3>::north()).rank == domain_2_bnw_patch->rank);
			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_2_bsw_tnw_patch->rank);

			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_2_bsw_bnw_patch->rank);
			CHECK(domain_2_bsw_bne_patch->getCoarseNbrInfo(Side<3>::east()).rank == domain_2_bse_patch->rank);
			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_2_bsw_bse_patch->rank);
			CHECK(domain_2_bsw_bne_patch->getCoarseNbrInfo(Side<3>::north()).rank == domain_2_bnw_patch->rank);
			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_2_bsw_tne_patch->rank);

			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_2_bsw_tse_patch->rank);
			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_2_bsw_tnw_patch->rank);
			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_2_bsw_bsw_patch->rank);
			CHECK(domain_2_bsw_tsw_patch->getCoarseNbrInfo(Side<3>::top()).rank == domain_2_tsw_patch->rank);

			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_2_bsw_tsw_patch->rank);
			CHECK(domain_2_bsw_tse_patch->getCoarseNbrInfo(Side<3>::east()).rank == domain_2_bse_patch->rank);
			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_2_bsw_tne_patch->rank);
			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_2_bsw_bse_patch->rank);
			CHECK(domain_2_bsw_tse_patch->getCoarseNbrInfo(Side<3>::top()).rank == domain_2_tsw_patch->rank);

			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_2_bsw_tne_patch->rank);
			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_2_bsw_tsw_patch->rank);
			CHECK(domain_2_bsw_tnw_patch->getCoarseNbrInfo(Side<3>::north()).rank == domain_2_bnw_patch->rank);
			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_2_bsw_bnw_patch->rank);
			CHECK(domain_2_bsw_tnw_patch->getCoarseNbrInfo(Side<3>::top()).rank == domain_2_tsw_patch->rank);

			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_2_bsw_tnw_patch->rank);
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Side<3>::east()).rank == domain_2_bse_patch->rank);
			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_2_bsw_tse_patch->rank);
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Side<3>::north()).rank == domain_2_bnw_patch->rank);
			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_2_bsw_bne_patch->rank);
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Side<3>::top()).rank == domain_2_tsw_patch->rank);

			CHECK(domain_2_bse_patch->getFineNbrInfo(Side<3>::west()).ranks[0] == domain_2_bsw_bse_patch->rank);
			CHECK(domain_2_bse_patch->getFineNbrInfo(Side<3>::west()).ranks[1] == domain_2_bsw_bne_patch->rank);
			CHECK(domain_2_bse_patch->getFineNbrInfo(Side<3>::west()).ranks[2] == domain_2_bsw_tse_patch->rank);
			CHECK(domain_2_bse_patch->getFineNbrInfo(Side<3>::west()).ranks[3] == domain_2_bsw_tne_patch->rank);
			CHECK(domain_2_bse_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_2_bne_patch->rank);
			CHECK(domain_2_bse_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_2_tse_patch->rank);

			CHECK(domain_2_bnw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_2_bne_patch->rank);
			CHECK(domain_2_bnw_patch->getFineNbrInfo(Side<3>::south()).ranks[0] == domain_2_bsw_bnw_patch->rank);
			CHECK(domain_2_bnw_patch->getFineNbrInfo(Side<3>::south()).ranks[1] == domain_2_bsw_bne_patch->rank);
			CHECK(domain_2_bnw_patch->getFineNbrInfo(Side<3>::south()).ranks[2] == domain_2_bsw_tnw_patch->rank);
			CHECK(domain_2_bnw_patch->getFineNbrInfo(Side<3>::south()).ranks[3] == domain_2_bsw_tne_patch->rank);
			CHECK(domain_2_bnw_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_2_tnw_patch->rank);

			CHECK(domain_2_bne_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_2_bnw_patch->rank);
			CHECK(domain_2_bne_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_2_bse_patch->rank);
			CHECK(domain_2_bne_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_2_tne_patch->rank);

			CHECK(domain_2_tsw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_2_tse_patch->rank);
			CHECK(domain_2_tsw_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_2_tnw_patch->rank);
			CHECK(domain_2_tsw_patch->getFineNbrInfo(Side<3>::bottom()).ranks[0] == domain_2_bsw_tsw_patch->rank);
			CHECK(domain_2_tsw_patch->getFineNbrInfo(Side<3>::bottom()).ranks[1] == domain_2_bsw_tse_patch->rank);
			CHECK(domain_2_tsw_patch->getFineNbrInfo(Side<3>::bottom()).ranks[2] == domain_2_bsw_tnw_patch->rank);
			CHECK(domain_2_tsw_patch->getFineNbrInfo(Side<3>::bottom()).ranks[3] == domain_2_bsw_tne_patch->rank);

			CHECK(domain_2_tse_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_2_tsw_patch->rank);
			CHECK(domain_2_tse_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_2_tne_patch->rank);
			CHECK(domain_2_tse_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_2_bse_patch->rank);

			CHECK(domain_2_tnw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_2_tne_patch->rank);
			CHECK(domain_2_tnw_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_2_tsw_patch->rank);
			CHECK(domain_2_tnw_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_2_bnw_patch->rank);

			CHECK(domain_2_tne_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_2_tnw_patch->rank);
			CHECK(domain_2_tne_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_2_tse_patch->rank);
			CHECK(domain_2_tne_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_2_bne_patch->rank);

			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_1_bse_patch->rank);
			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_1_bnw_patch->rank);
			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_1_tsw_patch->rank);

			CHECK(domain_1_bse_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_1_bsw_patch->rank);
			CHECK(domain_1_bse_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_1_bne_patch->rank);
			CHECK(domain_1_bse_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_1_tse_patch->rank);

			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_1_bne_patch->rank);
			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_1_bsw_patch->rank);
			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_1_tnw_patch->rank);

			CHECK(domain_1_bne_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_1_bnw_patch->rank);
			CHECK(domain_1_bne_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_1_bse_patch->rank);
			CHECK(domain_1_bne_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_1_tne_patch->rank);

			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_1_tse_patch->rank);
			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_1_tnw_patch->rank);
			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_1_bsw_patch->rank);

			CHECK(domain_1_tse_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_1_tsw_patch->rank);
			CHECK(domain_1_tse_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_1_tne_patch->rank);
			CHECK(domain_1_tse_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_1_bse_patch->rank);

			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_1_tne_patch->rank);
			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_1_tsw_patch->rank);
			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_1_bnw_patch->rank);

			CHECK(domain_1_tne_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_1_tnw_patch->rank);
			CHECK(domain_1_tne_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_1_tse_patch->rank);
			CHECK(domain_1_tne_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_1_bne_patch->rank);
		}
	}
	SECTION("orth_on_coarse")
	{
		if (rank == 0) {
			CHECK(domain_2_bsw_bse_patch->getCoarseNbrInfo(Side<3>::east()).orth_on_coarse == Orthant<2>::sw());
			CHECK(domain_2_bsw_bne_patch->getCoarseNbrInfo(Side<3>::east()).orth_on_coarse == Orthant<2>::se());
			CHECK(domain_2_bsw_tse_patch->getCoarseNbrInfo(Side<3>::east()).orth_on_coarse == Orthant<2>::nw());
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Side<3>::east()).orth_on_coarse == Orthant<2>::ne());

			CHECK(domain_2_bsw_bnw_patch->getCoarseNbrInfo(Side<3>::north()).orth_on_coarse == Orthant<2>::sw());
			CHECK(domain_2_bsw_bne_patch->getCoarseNbrInfo(Side<3>::north()).orth_on_coarse == Orthant<2>::se());
			CHECK(domain_2_bsw_tnw_patch->getCoarseNbrInfo(Side<3>::north()).orth_on_coarse == Orthant<2>::nw());
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Side<3>::north()).orth_on_coarse == Orthant<2>::ne());

			CHECK(domain_2_bsw_tsw_patch->getCoarseNbrInfo(Side<3>::top()).orth_on_coarse == Orthant<2>::sw());
			CHECK(domain_2_bsw_tse_patch->getCoarseNbrInfo(Side<3>::top()).orth_on_coarse == Orthant<2>::se());
			CHECK(domain_2_bsw_tnw_patch->getCoarseNbrInfo(Side<3>::top()).orth_on_coarse == Orthant<2>::nw());
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Side<3>::top()).orth_on_coarse == Orthant<2>::ne());
		}
	}
	SECTION("edge_nbr_info ids are correct")
	{
		if (rank == 0) {
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Edge::tn()).id == domain_2_bsw_tnw_patch->id);
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Edge::te()).id == domain_2_bsw_tse_patch->id);
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Edge::ne()).id == domain_2_bsw_bne_patch->id);
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_bsw_bse_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Edge::tn()).id == domain_2_bsw_tne_patch->id);
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::be()));
			CHECK(domain_2_bsw_bse_patch->hasNbr(Edge::tw()));
			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Edge::tw()).id == domain_2_bsw_tsw_patch->id);
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::se()));
			CHECK(domain_2_bsw_bse_patch->hasNbr(Edge::nw()));
			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Edge::nw()).id == domain_2_bsw_bnw_patch->id);

			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Edge::ts()).id == domain_2_bsw_tsw_patch->id);
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Edge::te()).id == domain_2_bsw_tne_patch->id);
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Edge::se()));
			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Edge::se()).id == domain_2_bsw_bse_patch->id);
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_bsw_bne_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Edge::ts()).id == domain_2_bsw_tse_patch->id);
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::be()));
			CHECK(domain_2_bsw_bne_patch->hasNbr(Edge::tw()));
			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Edge::tw()).id == domain_2_bsw_tnw_patch->id);
			CHECK(domain_2_bsw_bne_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Edge::sw()).id == domain_2_bsw_bsw_patch->id);
			CHECK(domain_2_bsw_bne_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bsw_bne_patch->getCoarseNbrInfo(Edge::ne()).id == domain_2_bne_patch->id);
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Edge::bn()).id == domain_2_bsw_bnw_patch->id);
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Edge::be()));
			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Edge::be()).id == domain_2_bsw_bse_patch->id);
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Edge::ne()).id == domain_2_bsw_tne_patch->id);
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bsw_tse_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Edge::bn()).id == domain_2_bsw_bne_patch->id);
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_bsw_tse_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Edge::bw()).id == domain_2_bsw_bsw_patch->id);
			CHECK(domain_2_bsw_tse_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bsw_tse_patch->getCoarseNbrInfo(Edge::te()).id == domain_2_tse_patch->id);
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::se()));
			CHECK(domain_2_bsw_tse_patch->hasNbr(Edge::nw()));
			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Edge::nw()).id == domain_2_bsw_tnw_patch->id);

			CHECK(domain_2_bsw_tnw_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Edge::bs()).id == domain_2_bsw_bsw_patch->id);
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bsw_tnw_patch->getCoarseNbrInfo(Edge::tn()).id == domain_2_tnw_patch->id);
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Edge::be()));
			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Edge::be()).id == domain_2_bsw_bne_patch->id);
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Edge::se()));
			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Edge::se()).id == domain_2_bsw_tse_patch->id);
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::nw()));

			CHECK(domain_2_bsw_tne_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Edge::bs()).id == domain_2_bsw_bse_patch->id);
			CHECK(domain_2_bsw_tne_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Edge::tn()).id == domain_2_tnw_patch->id);
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_bsw_tne_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Edge::bw()).id == domain_2_bsw_bnw_patch->id);
			CHECK(domain_2_bsw_tne_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Edge::te()).id == domain_2_tse_patch->id);
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Edge::tw()));
			CHECK(domain_2_bsw_tne_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Edge::sw()).id == domain_2_bsw_tsw_patch->id);
			CHECK(domain_2_bsw_tne_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Edge::ne()).id == domain_2_bne_patch->id);
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_bse_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bse_patch->getNormalNbrInfo(Edge::tn()).id == domain_2_tne_patch->id);
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::be()));
			CHECK(domain_2_bse_patch->hasNbr(Edge::tw()));
			CHECK(domain_2_bse_patch->getNormalNbrInfo(Edge::tw()).id == domain_2_tsw_patch->id);
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::se()));
			CHECK(domain_2_bse_patch->hasNbr(Edge::nw()));
			CHECK(domain_2_bse_patch->getNormalNbrInfo(Edge::nw()).id == domain_2_bnw_patch->id);

			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_bnw_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_bnw_patch->getNormalNbrInfo(Edge::ts()).id == domain_2_tsw_patch->id);
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_bnw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bnw_patch->getNormalNbrInfo(Edge::te()).id == domain_2_tne_patch->id);
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bnw_patch->hasNbr(Edge::se()));
			CHECK(domain_2_bnw_patch->getNormalNbrInfo(Edge::se()).id == domain_2_bse_patch->id);
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_bne_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_bne_patch->getNormalNbrInfo(Edge::ts()).id == domain_2_tse_patch->id);
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::be()));
			CHECK(domain_2_bne_patch->hasNbr(Edge::tw()));
			CHECK(domain_2_bne_patch->getNormalNbrInfo(Edge::tw()).id == domain_2_tnw_patch->id);
			CHECK(domain_2_bne_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_bne_patch->getFineNbrInfo(Edge::sw()).ids[0] == domain_2_bsw_bne_patch->id);
			CHECK(domain_2_bne_patch->getFineNbrInfo(Edge::sw()).ids[1] == domain_2_bsw_tne_patch->id);
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_tsw_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_tsw_patch->getNormalNbrInfo(Edge::bn()).id == domain_2_bnw_patch->id);
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_tsw_patch->hasNbr(Edge::be()));
			CHECK(domain_2_tsw_patch->getNormalNbrInfo(Edge::be()).id == domain_2_bse_patch->id);
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_tsw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_tsw_patch->getNormalNbrInfo(Edge::ne()).id == domain_2_tne_patch->id);
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_tse_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_tse_patch->getNormalNbrInfo(Edge::bn()).id == domain_2_bne_patch->id);
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_tse_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_tse_patch->getFineNbrInfo(Edge::bw()).ids[0] == domain_2_bsw_tse_patch->id);
			CHECK(domain_2_tse_patch->getFineNbrInfo(Edge::bw()).ids[1] == domain_2_bsw_tne_patch->id);
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::se()));
			CHECK(domain_2_tse_patch->hasNbr(Edge::nw()));
			CHECK(domain_2_tse_patch->getNormalNbrInfo(Edge::nw()).id == domain_2_tnw_patch->id);

			CHECK(domain_2_tnw_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_tnw_patch->getFineNbrInfo(Edge::bs()).ids[0] == domain_2_bsw_tnw_patch->id);
			CHECK(domain_2_tnw_patch->getFineNbrInfo(Edge::bs()).ids[1] == domain_2_bsw_tne_patch->id);
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_tnw_patch->hasNbr(Edge::be()));
			CHECK(domain_2_tnw_patch->getNormalNbrInfo(Edge::be()).id == domain_2_bne_patch->id);
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_tnw_patch->hasNbr(Edge::se()));
			CHECK(domain_2_tnw_patch->getNormalNbrInfo(Edge::se()).id == domain_2_tse_patch->id);
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::nw()));

			CHECK(domain_2_tne_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_tne_patch->getNormalNbrInfo(Edge::bs()).id == domain_2_bse_patch->id);
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_tne_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_tne_patch->getNormalNbrInfo(Edge::bw()).id == domain_2_bnw_patch->id);
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::tw()));
			CHECK(domain_2_tne_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_tne_patch->getNormalNbrInfo(Edge::sw()).id == domain_2_tsw_patch->id);
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::bs()));
			CHECK(domain_1_bsw_patch->hasNbr(Edge::tn()));
			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Edge::tn()).id == domain_1_tnw_patch->id);
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::bw()));
			CHECK(domain_1_bsw_patch->hasNbr(Edge::te()));
			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Edge::te()).id == domain_1_tse_patch->id);
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::sw()));
			CHECK(domain_1_bsw_patch->hasNbr(Edge::ne()));
			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Edge::ne()).id == domain_1_bne_patch->id);
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::bs()));
			CHECK(domain_1_bse_patch->hasNbr(Edge::tn()));
			CHECK(domain_1_bse_patch->getNormalNbrInfo(Edge::tn()).id == domain_1_tne_patch->id);
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::be()));
			CHECK(domain_1_bse_patch->hasNbr(Edge::tw()));
			CHECK(domain_1_bse_patch->getNormalNbrInfo(Edge::tw()).id == domain_1_tsw_patch->id);
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::se()));
			CHECK(domain_1_bse_patch->hasNbr(Edge::nw()));
			CHECK(domain_1_bse_patch->getNormalNbrInfo(Edge::nw()).id == domain_1_bnw_patch->id);

			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::bn()));
			CHECK(domain_1_bnw_patch->hasNbr(Edge::ts()));
			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Edge::ts()).id == domain_1_tsw_patch->id);
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::bw()));
			CHECK(domain_1_bnw_patch->hasNbr(Edge::te()));
			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Edge::te()).id == domain_1_tne_patch->id);
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::ne()));
			CHECK(domain_1_bnw_patch->hasNbr(Edge::se()));
			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Edge::se()).id == domain_1_bse_patch->id);
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::bn()));
			CHECK(domain_1_bne_patch->hasNbr(Edge::ts()));
			CHECK(domain_1_bne_patch->getNormalNbrInfo(Edge::ts()).id == domain_1_tse_patch->id);
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::be()));
			CHECK(domain_1_bne_patch->hasNbr(Edge::tw()));
			CHECK(domain_1_bne_patch->getNormalNbrInfo(Edge::tw()).id == domain_1_tnw_patch->id);
			CHECK(domain_1_bne_patch->hasNbr(Edge::sw()));
			CHECK(domain_1_bne_patch->getNormalNbrInfo(Edge::sw()).id == domain_1_bsw_patch->id);
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::tn()));
			CHECK(domain_1_tsw_patch->hasNbr(Edge::bn()));
			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Edge::bn()).id == domain_1_bnw_patch->id);
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::te()));
			CHECK(domain_1_tsw_patch->hasNbr(Edge::be()));
			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Edge::be()).id == domain_1_bse_patch->id);
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::sw()));
			CHECK(domain_1_tsw_patch->hasNbr(Edge::ne()));
			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Edge::ne()).id == domain_1_tne_patch->id);
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::tn()));
			CHECK(domain_1_tse_patch->hasNbr(Edge::bn()));
			CHECK(domain_1_tse_patch->getNormalNbrInfo(Edge::bn()).id == domain_1_bne_patch->id);
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::ts()));
			CHECK(domain_1_tse_patch->hasNbr(Edge::bw()));
			CHECK(domain_1_tse_patch->getNormalNbrInfo(Edge::bw()).id == domain_1_bsw_patch->id);
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::se()));
			CHECK(domain_1_tse_patch->hasNbr(Edge::nw()));
			CHECK(domain_1_tse_patch->getNormalNbrInfo(Edge::nw()).id == domain_1_tnw_patch->id);

			CHECK(domain_1_tnw_patch->hasNbr(Edge::bs()));
			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Edge::bs()).id == domain_1_bsw_patch->id);
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::te()));
			CHECK(domain_1_tnw_patch->hasNbr(Edge::be()));
			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Edge::be()).id == domain_1_bne_patch->id);
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::ne()));
			CHECK(domain_1_tnw_patch->hasNbr(Edge::se()));
			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Edge::se()).id == domain_1_tse_patch->id);
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::nw()));

			CHECK(domain_1_tne_patch->hasNbr(Edge::bs()));
			CHECK(domain_1_tne_patch->getNormalNbrInfo(Edge::bs()).id == domain_1_bse_patch->id);
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::ts()));
			CHECK(domain_1_tne_patch->hasNbr(Edge::bw()));
			CHECK(domain_1_tne_patch->getNormalNbrInfo(Edge::bw()).id == domain_1_bnw_patch->id);
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::tw()));
			CHECK(domain_1_tne_patch->hasNbr(Edge::sw()));
			CHECK(domain_1_tne_patch->getNormalNbrInfo(Edge::sw()).id == domain_1_tsw_patch->id);
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::nw()));
		}
	}
	SECTION("edge_nbr_info ranks are correct")
	{
		if (rank == 0) {
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Edge::tn()).rank == domain_2_bsw_tnw_patch->rank);
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Edge::te()).rank == domain_2_bsw_tse_patch->rank);
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Edge::ne()).rank == domain_2_bsw_bne_patch->rank);
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_bsw_bse_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Edge::tn()).rank == domain_2_bsw_tne_patch->rank);
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::be()));
			CHECK(domain_2_bsw_bse_patch->hasNbr(Edge::tw()));
			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Edge::tw()).rank == domain_2_bsw_tsw_patch->rank);
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Edge::se()));
			CHECK(domain_2_bsw_bse_patch->hasNbr(Edge::nw()));
			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Edge::nw()).rank == domain_2_bsw_bnw_patch->rank);

			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Edge::ts()).rank == domain_2_bsw_tsw_patch->rank);
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Edge::te()).rank == domain_2_bsw_tne_patch->rank);
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Edge::se()));
			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Edge::se()).rank == domain_2_bsw_bse_patch->rank);
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_bsw_bne_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Edge::ts()).rank == domain_2_bsw_tse_patch->rank);
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::be()));
			CHECK(domain_2_bsw_bne_patch->hasNbr(Edge::tw()));
			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Edge::tw()).rank == domain_2_bsw_tnw_patch->rank);
			CHECK(domain_2_bsw_bne_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Edge::sw()).rank == domain_2_bsw_bsw_patch->rank);
			CHECK(domain_2_bsw_bne_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bsw_bne_patch->getCoarseNbrInfo(Edge::ne()).rank == domain_2_bne_patch->rank);
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Edge::bn()).rank == domain_2_bsw_bnw_patch->rank);
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Edge::be()));
			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Edge::be()).rank == domain_2_bsw_bse_patch->rank);
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Edge::ne()).rank == domain_2_bsw_tne_patch->rank);
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bsw_tse_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Edge::bn()).rank == domain_2_bsw_bne_patch->rank);
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_bsw_tse_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Edge::bw()).rank == domain_2_bsw_bsw_patch->rank);
			CHECK(domain_2_bsw_tse_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bsw_tse_patch->getCoarseNbrInfo(Edge::te()).rank == domain_2_tse_patch->rank);
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Edge::se()));
			CHECK(domain_2_bsw_tse_patch->hasNbr(Edge::nw()));
			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Edge::nw()).rank == domain_2_bsw_tnw_patch->rank);

			CHECK(domain_2_bsw_tnw_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Edge::bs()).rank == domain_2_bsw_bsw_patch->rank);
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bsw_tnw_patch->getCoarseNbrInfo(Edge::tn()).rank == domain_2_tnw_patch->rank);
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Edge::be()));
			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Edge::be()).rank == domain_2_bsw_bne_patch->rank);
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Edge::se()));
			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Edge::se()).rank == domain_2_bsw_tse_patch->rank);
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Edge::nw()));

			CHECK(domain_2_bsw_tne_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Edge::bs()).rank == domain_2_bsw_bse_patch->rank);
			CHECK(domain_2_bsw_tne_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Edge::tn()).rank == domain_2_tnw_patch->rank);
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_bsw_tne_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Edge::bw()).rank == domain_2_bsw_bnw_patch->rank);
			CHECK(domain_2_bsw_tne_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Edge::te()).rank == domain_2_tse_patch->rank);
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Edge::tw()));
			CHECK(domain_2_bsw_tne_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Edge::sw()).rank == domain_2_bsw_tse_patch->rank);
			CHECK(domain_2_bsw_tne_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Edge::ne()).rank == domain_2_bne_patch->rank);
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_bse_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_bse_patch->getNormalNbrInfo(Edge::tn()).rank == domain_2_tne_patch->rank);
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::be()));
			CHECK(domain_2_bse_patch->hasNbr(Edge::tw()));
			CHECK(domain_2_bse_patch->getNormalNbrInfo(Edge::tw()).rank == domain_2_tsw_patch->rank);
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Edge::se()));
			CHECK(domain_2_bse_patch->hasNbr(Edge::nw()));
			CHECK(domain_2_bse_patch->getNormalNbrInfo(Edge::nw()).rank == domain_2_bnw_patch->rank);

			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_bnw_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_bnw_patch->getNormalNbrInfo(Edge::ts()).rank == domain_2_tsw_patch->rank);
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_bnw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_bnw_patch->getNormalNbrInfo(Edge::te()).rank == domain_2_tne_patch->rank);
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_bnw_patch->hasNbr(Edge::se()));
			CHECK(domain_2_bnw_patch->getNormalNbrInfo(Edge::se()).rank == domain_2_bse_patch->rank);
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_bne_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_bne_patch->getNormalNbrInfo(Edge::ts()).rank == domain_2_tse_patch->rank);
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::be()));
			CHECK(domain_2_bne_patch->hasNbr(Edge::tw()));
			CHECK(domain_2_bne_patch->getNormalNbrInfo(Edge::tw()).rank == domain_2_tnw_patch->rank);
			CHECK(domain_2_bne_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_bne_patch->getFineNbrInfo(Edge::sw()).ranks[0] == domain_2_bsw_bne_patch->rank);
			CHECK(domain_2_bne_patch->getFineNbrInfo(Edge::sw()).ranks[1] == domain_2_bsw_tne_patch->rank);
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_tsw_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_tsw_patch->getNormalNbrInfo(Edge::bn()).rank == domain_2_bnw_patch->rank);
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_tsw_patch->hasNbr(Edge::be()));
			CHECK(domain_2_tsw_patch->getNormalNbrInfo(Edge::be()).rank == domain_2_bse_patch->rank);
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_tsw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_tsw_patch->getNormalNbrInfo(Edge::ne()).rank == domain_2_tne_patch->rank);
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::tn()));
			CHECK(domain_2_tse_patch->hasNbr(Edge::bn()));
			CHECK(domain_2_tse_patch->getNormalNbrInfo(Edge::bn()).rank == domain_2_bne_patch->rank);
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_tse_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_tse_patch->getFineNbrInfo(Edge::bw()).ranks[0] == domain_2_bsw_tse_patch->rank);
			CHECK(domain_2_tse_patch->getFineNbrInfo(Edge::bw()).ranks[1] == domain_2_bsw_tne_patch->rank);
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Edge::se()));
			CHECK(domain_2_tse_patch->hasNbr(Edge::nw()));
			CHECK(domain_2_tse_patch->getNormalNbrInfo(Edge::nw()).rank == domain_2_tnw_patch->rank);

			CHECK(domain_2_tnw_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_tnw_patch->getFineNbrInfo(Edge::bs()).ranks[0] == domain_2_bsw_tnw_patch->rank);
			CHECK(domain_2_tnw_patch->getFineNbrInfo(Edge::bs()).ranks[1] == domain_2_bsw_tne_patch->rank);
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::te()));
			CHECK(domain_2_tnw_patch->hasNbr(Edge::be()));
			CHECK(domain_2_tnw_patch->getNormalNbrInfo(Edge::be()).rank == domain_2_bne_patch->rank);
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::ne()));
			CHECK(domain_2_tnw_patch->hasNbr(Edge::se()));
			CHECK(domain_2_tnw_patch->getNormalNbrInfo(Edge::se()).rank == domain_2_tse_patch->rank);
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Edge::nw()));

			CHECK(domain_2_tne_patch->hasNbr(Edge::bs()));
			CHECK(domain_2_tne_patch->getNormalNbrInfo(Edge::bs()).rank == domain_2_bse_patch->rank);
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::ts()));
			CHECK(domain_2_tne_patch->hasNbr(Edge::bw()));
			CHECK(domain_2_tne_patch->getNormalNbrInfo(Edge::bw()).rank == domain_2_bnw_patch->rank);
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::tw()));
			CHECK(domain_2_tne_patch->hasNbr(Edge::sw()));
			CHECK(domain_2_tne_patch->getNormalNbrInfo(Edge::sw()).rank == domain_2_tsw_patch->rank);
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::bs()));
			CHECK(domain_1_bsw_patch->hasNbr(Edge::tn()));
			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Edge::tn()).rank == domain_1_tnw_patch->rank);
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::bw()));
			CHECK(domain_1_bsw_patch->hasNbr(Edge::te()));
			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Edge::te()).rank == domain_1_tse_patch->rank);
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::sw()));
			CHECK(domain_1_bsw_patch->hasNbr(Edge::ne()));
			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Edge::ne()).rank == domain_1_bne_patch->rank);
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::bs()));
			CHECK(domain_1_bse_patch->hasNbr(Edge::tn()));
			CHECK(domain_1_bse_patch->getNormalNbrInfo(Edge::tn()).rank == domain_1_tne_patch->rank);
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::be()));
			CHECK(domain_1_bse_patch->hasNbr(Edge::tw()));
			CHECK(domain_1_bse_patch->getNormalNbrInfo(Edge::tw()).rank == domain_1_tsw_patch->rank);
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Edge::se()));
			CHECK(domain_1_bse_patch->hasNbr(Edge::nw()));
			CHECK(domain_1_bse_patch->getNormalNbrInfo(Edge::nw()).rank == domain_1_bnw_patch->rank);

			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::bn()));
			CHECK(domain_1_bnw_patch->hasNbr(Edge::ts()));
			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Edge::ts()).rank == domain_1_tsw_patch->rank);
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::bw()));
			CHECK(domain_1_bnw_patch->hasNbr(Edge::te()));
			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Edge::te()).rank == domain_1_tne_patch->rank);
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::ne()));
			CHECK(domain_1_bnw_patch->hasNbr(Edge::se()));
			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Edge::se()).rank == domain_1_bse_patch->rank);
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::bn()));
			CHECK(domain_1_bne_patch->hasNbr(Edge::ts()));
			CHECK(domain_1_bne_patch->getNormalNbrInfo(Edge::ts()).rank == domain_1_tse_patch->rank);
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::be()));
			CHECK(domain_1_bne_patch->hasNbr(Edge::tw()));
			CHECK(domain_1_bne_patch->getNormalNbrInfo(Edge::tw()).rank == domain_1_tnw_patch->rank);
			CHECK(domain_1_bne_patch->hasNbr(Edge::sw()));
			CHECK(domain_1_bne_patch->getNormalNbrInfo(Edge::sw()).rank == domain_1_bsw_patch->rank);
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::tn()));
			CHECK(domain_1_tsw_patch->hasNbr(Edge::bn()));
			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Edge::bn()).rank == domain_1_bnw_patch->rank);
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::te()));
			CHECK(domain_1_tsw_patch->hasNbr(Edge::be()));
			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Edge::be()).rank == domain_1_bse_patch->rank);
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::sw()));
			CHECK(domain_1_tsw_patch->hasNbr(Edge::ne()));
			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Edge::ne()).rank == domain_1_tne_patch->rank);
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::tn()));
			CHECK(domain_1_tse_patch->hasNbr(Edge::bn()));
			CHECK(domain_1_tse_patch->getNormalNbrInfo(Edge::bn()).rank == domain_1_bne_patch->rank);
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::ts()));
			CHECK(domain_1_tse_patch->hasNbr(Edge::bw()));
			CHECK(domain_1_tse_patch->getNormalNbrInfo(Edge::bw()).rank == domain_1_bsw_patch->rank);
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Edge::se()));
			CHECK(domain_1_tse_patch->hasNbr(Edge::nw()));
			CHECK(domain_1_tse_patch->getNormalNbrInfo(Edge::nw()).rank == domain_1_tnw_patch->rank);

			CHECK(domain_1_tnw_patch->hasNbr(Edge::bs()));
			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Edge::bs()).rank == domain_1_bsw_patch->rank);
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::te()));
			CHECK(domain_1_tnw_patch->hasNbr(Edge::be()));
			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Edge::be()).rank == domain_1_bne_patch->rank);
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::ne()));
			CHECK(domain_1_tnw_patch->hasNbr(Edge::se()));
			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Edge::se()).rank == domain_1_tse_patch->rank);
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Edge::nw()));

			CHECK(domain_1_tne_patch->hasNbr(Edge::bs()));
			CHECK(domain_1_tne_patch->getNormalNbrInfo(Edge::bs()).rank == domain_1_bse_patch->rank);
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::ts()));
			CHECK(domain_1_tne_patch->hasNbr(Edge::bw()));
			CHECK(domain_1_tne_patch->getNormalNbrInfo(Edge::bw()).rank == domain_1_bnw_patch->rank);
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::tw()));
			CHECK(domain_1_tne_patch->hasNbr(Edge::sw()));
			CHECK(domain_1_tne_patch->getNormalNbrInfo(Edge::sw()).rank == domain_1_tsw_patch->rank);
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Edge::nw()));

			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::bs()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::tn()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::bn()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::ts()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::bw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::te()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::be()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::tw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::sw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::ne()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::se()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Edge::nw()));
		}
	}
	SECTION("edge_nbr_info ranks are correct")
	{
		if (rank == 0) {
			CHECK(domain_2_bsw_tnw_patch->getCoarseNbrInfo(Edge::tn()).orth_on_coarse == Orthant<1>::lower());
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Edge::tn()).orth_on_coarse == Orthant<1>::upper());

			CHECK(domain_2_bsw_tse_patch->getCoarseNbrInfo(Edge::te()).orth_on_coarse == Orthant<1>::lower());
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Edge::te()).orth_on_coarse == Orthant<1>::upper());

			CHECK(domain_2_bsw_bne_patch->getCoarseNbrInfo(Edge::ne()).orth_on_coarse == Orthant<1>::lower());
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Edge::ne()).orth_on_coarse == Orthant<1>::upper());
		}
	}
	SECTION("corner_nbr_info ids are correct")
	{
		if (rank == 0) {
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::tnw()));
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::tne()));
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Corner<3>::tne()).id == domain_2_bsw_tne_patch->id);

			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::tse()));
			CHECK(domain_2_bsw_bse_patch->hasNbr(Corner<3>::tnw()));
			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Corner<3>::tnw()).id == domain_2_bsw_tnw_patch->id);
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::tsw()));
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::tse()));
			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Corner<3>::tse()).id == domain_2_bsw_tse_patch->id);
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::bne()));
			CHECK(domain_2_bsw_bne_patch->hasNbr(Corner<3>::tsw()));
			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Corner<3>::tsw()).id == domain_2_bsw_tsw_patch->id);
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::bnw()));
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::bne()));
			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Corner<3>::bne()).id == domain_2_bsw_bne_patch->id);
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::bse()));
			CHECK(domain_2_bsw_tse_patch->hasNbr(Corner<3>::bnw()));
			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Corner<3>::bnw()).id == domain_2_bsw_bnw_patch->id);
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::bsw()));
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::bse()));
			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Corner<3>::bse()).id == domain_2_bsw_bse_patch->id);
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::tne()));

			CHECK(domain_2_bsw_tne_patch->hasNbr(Corner<3>::bsw()));
			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Corner<3>::bsw()).id == domain_2_bsw_bsw_patch->id);
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Corner<3>::tnw()));
			CHECK(domain_2_bsw_tne_patch->hasNbr(Corner<3>::tne()));
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Corner<3>::tne()).id == domain_2_tne_patch->id);

			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::tse()));
			CHECK(domain_2_bse_patch->hasNbr(Corner<3>::tnw()));
			CHECK(domain_2_bse_patch->getNormalNbrInfo(Corner<3>::tnw()).id == domain_2_tnw_patch->id);
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::tsw()));
			CHECK(domain_2_bnw_patch->hasNbr(Corner<3>::tse()));
			CHECK(domain_2_bnw_patch->getNormalNbrInfo(Corner<3>::tse()).id == domain_2_tse_patch->id);
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::bne()));
			CHECK(domain_2_bne_patch->hasNbr(Corner<3>::tsw()));
			CHECK(domain_2_bne_patch->getNormalNbrInfo(Corner<3>::tsw()).id == domain_2_tsw_patch->id);
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::bnw()));
			CHECK(domain_2_tsw_patch->hasNbr(Corner<3>::bne()));
			CHECK(domain_2_tsw_patch->getNormalNbrInfo(Corner<3>::bne()).id == domain_2_bne_patch->id);
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::bse()));
			CHECK(domain_2_tse_patch->hasNbr(Corner<3>::bnw()));
			CHECK(domain_2_tse_patch->getNormalNbrInfo(Corner<3>::bnw()).id == domain_2_bnw_patch->id);
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::bsw()));
			CHECK(domain_2_tnw_patch->hasNbr(Corner<3>::bse()));
			CHECK(domain_2_tnw_patch->getNormalNbrInfo(Corner<3>::bse()).id == domain_2_bse_patch->id);
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::tne()));

			CHECK(domain_2_tne_patch->hasNbr(Corner<3>::bsw()));
			CHECK(domain_2_tne_patch->getFineNbrInfo(Corner<3>::bsw()).ids[0] == domain_2_bsw_tne_patch->id);
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::tnw()));
			CHECK(domain_1_bsw_patch->hasNbr(Corner<3>::tne()));
			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Corner<3>::tne()).id == domain_1_tne_patch->id);

			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::tse()));
			CHECK(domain_1_bse_patch->hasNbr(Corner<3>::tnw()));
			CHECK(domain_1_bse_patch->getNormalNbrInfo(Corner<3>::tnw()).id == domain_1_tnw_patch->id);
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::tsw()));
			CHECK(domain_1_bnw_patch->hasNbr(Corner<3>::tse()));
			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Corner<3>::tse()).id == domain_1_tse_patch->id);
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::bne()));
			CHECK(domain_1_bne_patch->hasNbr(Corner<3>::tsw()));
			CHECK(domain_1_bne_patch->getNormalNbrInfo(Corner<3>::tsw()).id == domain_1_tsw_patch->id);
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::bnw()));
			CHECK(domain_1_tsw_patch->hasNbr(Corner<3>::bne()));
			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Corner<3>::bne()).id == domain_1_bne_patch->id);
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::bse()));
			CHECK(domain_1_tse_patch->hasNbr(Corner<3>::bnw()));
			CHECK(domain_1_tse_patch->getNormalNbrInfo(Corner<3>::bnw()).id == domain_1_bnw_patch->id);
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::bsw()));
			CHECK(domain_1_tnw_patch->hasNbr(Corner<3>::bse()));
			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Corner<3>::bse()).id == domain_1_bse_patch->id);
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::tne()));

			CHECK(domain_1_tne_patch->hasNbr(Corner<3>::bsw()));
			CHECK(domain_1_tne_patch->getNormalNbrInfo(Corner<3>::bsw()).id == domain_1_bsw_patch->id);
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::tne()));
		}
	}
	SECTION("corner_nbr_info ranks are correct")
	{
		if (rank == 0) {
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::tnw()));
			CHECK(domain_2_bsw_bsw_patch->hasNbr(Corner<3>::tne()));
			CHECK(domain_2_bsw_bsw_patch->getNormalNbrInfo(Corner<3>::tne()).rank == domain_2_bsw_tne_patch->rank);

			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::tse()));
			CHECK(domain_2_bsw_bse_patch->hasNbr(Corner<3>::tnw()));
			CHECK(domain_2_bsw_bse_patch->getNormalNbrInfo(Corner<3>::tnw()).rank == domain_2_bsw_tnw_patch->rank);
			CHECK_FALSE(domain_2_bsw_bse_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::tsw()));
			CHECK(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::tse()));
			CHECK(domain_2_bsw_bnw_patch->getNormalNbrInfo(Corner<3>::tse()).rank == domain_2_bsw_tse_patch->rank);
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bsw_bnw_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::bne()));
			CHECK(domain_2_bsw_bne_patch->hasNbr(Corner<3>::tsw()));
			CHECK(domain_2_bsw_bne_patch->getNormalNbrInfo(Corner<3>::tsw()).rank == domain_2_bsw_tsw_patch->rank);
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bsw_bne_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::bnw()));
			CHECK(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::bne()));
			CHECK(domain_2_bsw_tsw_patch->getNormalNbrInfo(Corner<3>::bne()).rank == domain_2_bsw_bne_patch->rank);
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bsw_tsw_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::bse()));
			CHECK(domain_2_bsw_tse_patch->hasNbr(Corner<3>::bnw()));
			CHECK(domain_2_bsw_tse_patch->getNormalNbrInfo(Corner<3>::bnw()).rank == domain_2_bsw_bnw_patch->rank);
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bsw_tse_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::bsw()));
			CHECK(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::bse()));
			CHECK(domain_2_bsw_tnw_patch->getNormalNbrInfo(Corner<3>::bse()).rank == domain_2_bsw_bse_patch->rank);
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bsw_tnw_patch->hasNbr(Corner<3>::tne()));

			CHECK(domain_2_bsw_tne_patch->hasNbr(Corner<3>::bsw()));
			CHECK(domain_2_bsw_tne_patch->getNormalNbrInfo(Corner<3>::bsw()).rank == domain_2_bsw_bsw_patch->rank);
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bsw_tne_patch->hasNbr(Corner<3>::tnw()));
			CHECK(domain_2_bsw_tne_patch->hasNbr(Corner<3>::tne()));
			CHECK(domain_2_bsw_tne_patch->getCoarseNbrInfo(Corner<3>::tne()).rank == domain_2_tne_patch->rank);

			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::tse()));
			CHECK(domain_2_bse_patch->hasNbr(Corner<3>::tnw()));
			CHECK(domain_2_bse_patch->getNormalNbrInfo(Corner<3>::tnw()).rank == domain_2_tnw_patch->rank);
			CHECK_FALSE(domain_2_bse_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::tsw()));
			CHECK(domain_2_bnw_patch->hasNbr(Corner<3>::tse()));
			CHECK(domain_2_bnw_patch->getNormalNbrInfo(Corner<3>::tse()).rank == domain_2_tse_patch->rank);
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bnw_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::bne()));
			CHECK(domain_2_bne_patch->hasNbr(Corner<3>::tsw()));
			CHECK(domain_2_bne_patch->getNormalNbrInfo(Corner<3>::tsw()).rank == domain_2_tsw_patch->rank);
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_bne_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::bnw()));
			CHECK(domain_2_tsw_patch->hasNbr(Corner<3>::bne()));
			CHECK(domain_2_tsw_patch->getNormalNbrInfo(Corner<3>::bne()).rank == domain_2_bne_patch->rank);
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_tsw_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::bse()));
			CHECK(domain_2_tse_patch->hasNbr(Corner<3>::bnw()));
			CHECK(domain_2_tse_patch->getNormalNbrInfo(Corner<3>::bnw()).rank == domain_2_bnw_patch->rank);
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_tse_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::bsw()));
			CHECK(domain_2_tnw_patch->hasNbr(Corner<3>::bse()));
			CHECK(domain_2_tnw_patch->getNormalNbrInfo(Corner<3>::bse()).rank == domain_2_bse_patch->rank);
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_tnw_patch->hasNbr(Corner<3>::tne()));

			CHECK(domain_2_tne_patch->hasNbr(Corner<3>::bsw()));
			CHECK(domain_2_tne_patch->getFineNbrInfo(Corner<3>::bsw()).ranks[0] == domain_2_bsw_tne_patch->rank);
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_2_tne_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_1_bsw_patch->hasNbr(Corner<3>::tnw()));
			CHECK(domain_1_bsw_patch->hasNbr(Corner<3>::tne()));
			CHECK(domain_1_bsw_patch->getNormalNbrInfo(Corner<3>::tne()).rank == domain_1_tne_patch->rank);

			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::tse()));
			CHECK(domain_1_bse_patch->hasNbr(Corner<3>::tnw()));
			CHECK(domain_1_bse_patch->getNormalNbrInfo(Corner<3>::tnw()).rank == domain_1_tnw_patch->rank);
			CHECK_FALSE(domain_1_bse_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::tsw()));
			CHECK(domain_1_bnw_patch->hasNbr(Corner<3>::tse()));
			CHECK(domain_1_bnw_patch->getNormalNbrInfo(Corner<3>::tse()).rank == domain_1_tse_patch->rank);
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_1_bnw_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::bne()));
			CHECK(domain_1_bne_patch->hasNbr(Corner<3>::tsw()));
			CHECK(domain_1_bne_patch->getNormalNbrInfo(Corner<3>::tsw()).rank == domain_1_tsw_patch->rank);
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_1_bne_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::bnw()));
			CHECK(domain_1_tsw_patch->hasNbr(Corner<3>::bne()));
			CHECK(domain_1_tsw_patch->getNormalNbrInfo(Corner<3>::bne()).rank == domain_1_bne_patch->rank);
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_1_tsw_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::bse()));
			CHECK(domain_1_tse_patch->hasNbr(Corner<3>::bnw()));
			CHECK(domain_1_tse_patch->getNormalNbrInfo(Corner<3>::bnw()).rank == domain_1_bnw_patch->rank);
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_1_tse_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::bsw()));
			CHECK(domain_1_tnw_patch->hasNbr(Corner<3>::bse()));
			CHECK(domain_1_tnw_patch->getNormalNbrInfo(Corner<3>::bse()).rank == domain_1_bse_patch->rank);
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_1_tnw_patch->hasNbr(Corner<3>::tne()));

			CHECK(domain_1_tne_patch->hasNbr(Corner<3>::bsw()));
			CHECK(domain_1_tne_patch->getNormalNbrInfo(Corner<3>::bsw()).rank == domain_1_bsw_patch->rank);
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_1_tne_patch->hasNbr(Corner<3>::tne()));

			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::bsw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::bse()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::bnw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::bne()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::tsw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::tse()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::tnw()));
			CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<3>::tne()));
		}
	}
}
