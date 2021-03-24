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

#include <ThunderEgg/P4estDomGen.h>

#include <p4est.h>
#include <p4est_mesh.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

#define MESHES                                                                                     \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"

TEST_CASE("SinglePatch", "[p4estDomGen]")
{
	p4est_connectivity_t *conn  = p4est_connectivity_new_unitsquare();
	p4est_t *             p4est = p4est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);
	IsNeumannFunc<2>      inf   = [](Side<2> s, const std::array<double, 2> &lower,
                              const std::array<double, 2> &upper) { return false; };
	P4estDomGen::BlockMapFunc bmf
	= [](int block_no, double unit_x, double unit_y, double &x, double &y) {
		  x = unit_x;
		  y = unit_y;
	  };
	P4estDomGen dg(p4est, {10, 10}, 1, inf, bmf);
	auto        domain = dg.getFinestDomain();
	CHECK(domain->getNumGlobalPatches() == 1);
	auto patch = domain->getPatchInfoVector()[0];
	CHECK_FALSE(patch->hasNbr(Side<2>::west()));
	CHECK_FALSE(patch->hasNbr(Side<2>::east()));
	CHECK_FALSE(patch->hasNbr(Side<2>::south()));
	CHECK_FALSE(patch->hasNbr(Side<2>::north()));
	CHECK(patch->spacings[0] == Catch::Approx(1.0 / 10));
	CHECK(patch->spacings[1] == Catch::Approx(1.0 / 10));
	CHECK(patch->starts[0] == Catch::Approx(0));
	CHECK(patch->starts[1] == Catch::Approx(0));
	CHECK(patch->ns[0] == 10);
	CHECK(patch->ns[1] == 10);

	CHECK_FALSE(dg.hasCoarserDomain());
}
TEST_CASE("P4estDomGen 2x2 Uniform", "[p4estDomGen]")
{
	p4est_connectivity_t *conn = p4est_connectivity_new_unitsquare();

	p4est_t *p4est = p4est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

	p4est_refine(
	p4est, false,
	[](p4est_t *p4est, p4est_topidx_t witch_tree, p4est_quadrant_t *quadrant) -> int { return 1; },
	nullptr);

	IsNeumannFunc<2> inf = [](Side<2> s, const std::array<double, 2> &lower,
	                          const std::array<double, 2> &upper) { return false; };

	int    nx              = GENERATE(5, 10);
	int    ny              = GENERATE(5, 10);
	double scale_x         = GENERATE(0.5, 1.0);
	double scale_y         = GENERATE(0.5, 1.0);
	int    num_ghost_cells = GENERATE(0, 1, 2);

	INFO("nx: " << nx);
	INFO("ny: " << ny);
	INFO("scale_x: " << scale_x);
	INFO("scale_x: " << scale_y);
	INFO("num_ghost_cells: " << num_ghost_cells);

	P4estDomGen::BlockMapFunc bmf
	= [&](int block_no, double unit_x, double unit_y, double &x, double &y) {
		  x = scale_x * unit_x;
		  y = scale_y * unit_y;
	  };

	P4estDomGen dg(p4est, {nx, ny}, num_ghost_cells, inf, bmf);

	auto domain         = dg.getFinestDomain();
	auto coarser_domain = dg.getCoarserDomain();

	SECTION("patches have correct spacings")
	{
		for (auto patch : domain->getPatchInfoVector()) {
			CHECK(patch->spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
			CHECK(patch->spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
		}

		auto patch = coarser_domain->getPatchInfoVector()[0];

		CHECK(patch->spacings[0] == Catch::Approx(scale_x * 1.0 / nx));
		CHECK(patch->spacings[1] == Catch::Approx(scale_y * 1.0 / ny));
	}
	SECTION("patches have correct ns")
	{
		for (auto patch : domain->getPatchInfoVector()) {
			CHECK(patch->ns[0] == nx);
			CHECK(patch->ns[1] == ny);
		}

		auto patch = coarser_domain->getPatchInfoVector()[0];

		CHECK(patch->ns[0] == nx);
		CHECK(patch->ns[1] == ny);
	}
	SECTION("patches have ranks set")
	{
		for (auto patch : domain->getPatchInfoVector()) {
			CHECK(patch->rank == 0);
		}

		auto patch = coarser_domain->getPatchInfoVector()[0];

		CHECK(patch->rank == 0);
	}
	SECTION("patches have refine_level set")
	{
		for (auto patch : domain->getPatchInfoVector()) {
			CHECK(patch->refine_level == 1);
		}

		auto patch = coarser_domain->getPatchInfoVector()[0];

		CHECK(patch->refine_level == 0);
	}
	SECTION("patches have num_ghost_cells set")
	{
		for (auto patch : domain->getPatchInfoVector()) {
			CHECK(patch->num_ghost_cells == num_ghost_cells);
		}

		auto patch = coarser_domain->getPatchInfoVector()[0];

		CHECK(patch->num_ghost_cells == num_ghost_cells);
	}
	std::shared_ptr<const PatchInfo<2>> sw_patch;
	std::shared_ptr<const PatchInfo<2>> se_patch;
	std::shared_ptr<const PatchInfo<2>> nw_patch;
	std::shared_ptr<const PatchInfo<2>> ne_patch;
	std::shared_ptr<const PatchInfo<2>> coarser_patch;

	for (auto patch : domain->getPatchInfoVector()) {
		double x = patch->starts[0];
		double y = patch->starts[1];
		if (x == Catch::Approx(0) && y == Catch::Approx(0)) {
			sw_patch = patch;
		}
		if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0)) {
			se_patch = patch;
		}
		if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y)) {
			nw_patch = patch;
		}
		if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
			ne_patch = patch;
		}
	}

	coarser_patch = coarser_domain->getPatchInfoVector()[0];

	SECTION("patches have starts set")
	{
		REQUIRE(sw_patch != nullptr);
		REQUIRE(se_patch != nullptr);
		REQUIRE(nw_patch != nullptr);
		REQUIRE(ne_patch != nullptr);

		CHECK(coarser_patch->starts[0] == Catch::Approx(0.0));
		CHECK(coarser_patch->starts[1] == Catch::Approx(0.0));
	}
	SECTION("parent ids are set correctly")
	{
		CHECK(sw_patch->parent_id == coarser_patch->id);
		CHECK(se_patch->parent_id == coarser_patch->id);
		CHECK(nw_patch->parent_id == coarser_patch->id);
		CHECK(ne_patch->parent_id == coarser_patch->id);

		CHECK(coarser_patch->parent_id == -1);
	}
	SECTION("child ids are set correctly")
	{
		for (int i = 0; i < 4; i++) {
			CHECK(sw_patch->child_ids[i] == -1);
			CHECK(se_patch->child_ids[i] == -1);
			CHECK(nw_patch->child_ids[i] == -1);
			CHECK(ne_patch->child_ids[i] == -1);
		}

		CHECK(coarser_patch->child_ids[0] == sw_patch->id);
		CHECK(coarser_patch->child_ids[1] == se_patch->id);
		CHECK(coarser_patch->child_ids[2] == nw_patch->id);
		CHECK(coarser_patch->child_ids[3] == ne_patch->id);
	}
	SECTION("orth on parent is set correctly")
	{
		CHECK(sw_patch->orth_on_parent == Orthant<2>::sw());
		CHECK(se_patch->orth_on_parent == Orthant<2>::se());
		CHECK(nw_patch->orth_on_parent == Orthant<2>::nw());
		CHECK(ne_patch->orth_on_parent == Orthant<2>::ne());

		CHECK(coarser_patch->orth_on_parent == Orthant<2>::null());
	}
	SECTION("parent ranks are set correctly")
	{
		CHECK(sw_patch->parent_rank == coarser_patch->rank);
		CHECK(se_patch->parent_rank == coarser_patch->rank);
		CHECK(nw_patch->parent_rank == coarser_patch->rank);
		CHECK(ne_patch->parent_rank == coarser_patch->rank);

		CHECK(coarser_patch->parent_rank == -1);
	}
	SECTION("child ranks are set correctly")
	{
		for (int i = 0; i < 4; i++) {
			CHECK(sw_patch->child_ranks[i] == -1);
			CHECK(se_patch->child_ranks[i] == -1);
			CHECK(nw_patch->child_ranks[i] == -1);
			CHECK(ne_patch->child_ranks[i] == -1);
		}

		CHECK(coarser_patch->child_ranks[0] == sw_patch->rank);
		CHECK(coarser_patch->child_ranks[1] == se_patch->rank);
		CHECK(coarser_patch->child_ranks[2] == nw_patch->rank);
		CHECK(coarser_patch->child_ranks[3] == ne_patch->rank);
	}
	SECTION("nbr_info ids are correct")
	{
		CHECK_FALSE(sw_patch->hasNbr(Side<2>::west()));
		CHECK(sw_patch->getNormalNbrInfo(Side<2>::east()).id == se_patch->id);
		CHECK_FALSE(sw_patch->hasNbr(Side<2>::south()));
		CHECK(sw_patch->getNormalNbrInfo(Side<2>::north()).id == nw_patch->id);

		CHECK(se_patch->getNormalNbrInfo(Side<2>::west()).id == sw_patch->id);
		CHECK_FALSE(se_patch->hasNbr(Side<2>::east()));
		CHECK_FALSE(se_patch->hasNbr(Side<2>::south()));
		CHECK(se_patch->getNormalNbrInfo(Side<2>::north()).id == ne_patch->id);

		CHECK_FALSE(nw_patch->hasNbr(Side<2>::west()));
		CHECK(nw_patch->getNormalNbrInfo(Side<2>::east()).id == ne_patch->id);
		CHECK(nw_patch->getNormalNbrInfo(Side<2>::south()).id == sw_patch->id);
		CHECK_FALSE(nw_patch->hasNbr(Side<2>::north()));

		CHECK(ne_patch->getNormalNbrInfo(Side<2>::west()).id == nw_patch->id);
		CHECK_FALSE(ne_patch->hasNbr(Side<2>::east()));
		CHECK(ne_patch->getNormalNbrInfo(Side<2>::south()).id == se_patch->id);
		CHECK_FALSE(ne_patch->hasNbr(Side<2>::north()));

		CHECK_FALSE(coarser_patch->hasNbr(Side<2>::west()));
		CHECK_FALSE(coarser_patch->hasNbr(Side<2>::east()));
		CHECK_FALSE(coarser_patch->hasNbr(Side<2>::south()));
		CHECK_FALSE(coarser_patch->hasNbr(Side<2>::north()));
	}
	SECTION("nbr_info ranks are correct")
	{
		CHECK_FALSE(sw_patch->hasNbr(Side<2>::west()));
		CHECK(sw_patch->getNormalNbrInfo(Side<2>::east()).rank == se_patch->rank);
		CHECK_FALSE(sw_patch->hasNbr(Side<2>::south()));
		CHECK(sw_patch->getNormalNbrInfo(Side<2>::north()).rank == nw_patch->rank);

		CHECK(se_patch->getNormalNbrInfo(Side<2>::west()).rank == sw_patch->rank);
		CHECK_FALSE(se_patch->hasNbr(Side<2>::east()));
		CHECK_FALSE(se_patch->hasNbr(Side<2>::south()));
		CHECK(se_patch->getNormalNbrInfo(Side<2>::north()).rank == ne_patch->rank);

		CHECK_FALSE(nw_patch->hasNbr(Side<2>::west()));
		CHECK(nw_patch->getNormalNbrInfo(Side<2>::east()).rank == ne_patch->rank);
		CHECK(nw_patch->getNormalNbrInfo(Side<2>::south()).rank == sw_patch->rank);
		CHECK_FALSE(nw_patch->hasNbr(Side<2>::north()));

		CHECK(ne_patch->getNormalNbrInfo(Side<2>::west()).rank == nw_patch->rank);
		CHECK_FALSE(ne_patch->hasNbr(Side<2>::east()));
		CHECK(ne_patch->getNormalNbrInfo(Side<2>::south()).rank == se_patch->rank);
		CHECK_FALSE(ne_patch->hasNbr(Side<2>::north()));

		CHECK_FALSE(coarser_patch->hasNbr(Side<2>::west()));
		CHECK_FALSE(coarser_patch->hasNbr(Side<2>::east()));
		CHECK_FALSE(coarser_patch->hasNbr(Side<2>::south()));
		CHECK_FALSE(coarser_patch->hasNbr(Side<2>::north()));
	}
}
TEST_CASE("2x1 brick", "[p4estDomGen]")
{
	p4est_connectivity_t *conn  = p4est_connectivity_new_brick(2, 1, false, false);
	p4est_t *             p4est = p4est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);
	IsNeumannFunc<2>      inf   = [](Side<2> s, const std::array<double, 2> &lower,
                              const std::array<double, 2> &upper) { return false; };
	P4estDomGen::BlockMapFunc bmf
	= [](int block_no, double unit_x, double unit_y, double &x, double &y) {
		  x = unit_x;
		  y = unit_y;
	  };
	P4estDomGen dg(p4est, {10, 10}, 1, inf, bmf);
	auto        domain = dg.getFinestDomain();
	CHECK(domain->getNumGlobalPatches() == 2);
	auto patch1 = domain->getPatchInfoVector()[0];
	CHECK_FALSE(patch1->hasNbr(Side<2>::west()));
	CHECK(patch1->hasNbr(Side<2>::east()));
	CHECK_FALSE(patch1->hasNbr(Side<2>::south()));
	CHECK_FALSE(patch1->hasNbr(Side<2>::north()));
	CHECK(patch1->spacings[0] == Catch::Approx(1.0 / 10));
	CHECK(patch1->spacings[1] == Catch::Approx(1.0 / 10));
	CHECK(patch1->starts[0] == Catch::Approx(0));
	CHECK(patch1->starts[1] == Catch::Approx(0));
	CHECK(patch1->ns[0] == 10);
	CHECK(patch1->ns[1] == 10);
	auto patch2 = domain->getPatchInfoVector()[1];
	CHECK(patch2->hasNbr(Side<2>::west()));
	CHECK_FALSE(patch2->hasNbr(Side<2>::east()));
	CHECK_FALSE(patch2->hasNbr(Side<2>::south()));
	CHECK_FALSE(patch2->hasNbr(Side<2>::north()));
	CHECK(patch2->spacings[0] == Catch::Approx(1.0 / 10));
	CHECK(patch2->spacings[1] == Catch::Approx(1.0 / 10));
	CHECK(patch2->starts[0] == Catch::Approx(0));
	CHECK(patch2->starts[1] == Catch::Approx(0));
	CHECK(patch2->ns[0] == 10);
	CHECK(patch2->ns[1] == 10);

	CHECK_FALSE(dg.hasCoarserDomain());
}