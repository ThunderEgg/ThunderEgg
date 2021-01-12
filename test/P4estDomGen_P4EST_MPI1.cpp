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
#include <ThunderEgg/P4estDomGen.h>
#include <p4est.h>
#include <p4est_mesh.h>
using namespace std;
using namespace ThunderEgg;
#define MESHES                                                                                     \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
TEST_CASE("SinglePatch", "[p4estDomGen]")
{
	//
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
	CHECK(patch->spacings[0] == Approx(1.0 / 10));
	CHECK(patch->spacings[1] == Approx(1.0 / 10));
	CHECK(patch->starts[0] == Approx(0));
	CHECK(patch->starts[1] == Approx(0));
	CHECK(patch->ns[0] == 10);
	CHECK(patch->ns[1] == 10);

	CHECK_FALSE(dg.hasCoarserDomain());
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
	CHECK(patch1->spacings[0] == Approx(1.0 / 10));
	CHECK(patch1->spacings[1] == Approx(1.0 / 10));
	CHECK(patch1->starts[0] == Approx(0));
	CHECK(patch1->starts[1] == Approx(0));
	CHECK(patch1->ns[0] == 10);
	CHECK(patch1->ns[1] == 10);
	auto patch2 = domain->getPatchInfoVector()[1];
	CHECK(patch2->hasNbr(Side<2>::west()));
	CHECK_FALSE(patch2->hasNbr(Side<2>::east()));
	CHECK_FALSE(patch2->hasNbr(Side<2>::south()));
	CHECK_FALSE(patch2->hasNbr(Side<2>::north()));
	CHECK(patch2->spacings[0] == Approx(1.0 / 10));
	CHECK(patch2->spacings[1] == Approx(1.0 / 10));
	CHECK(patch2->starts[0] == Approx(0));
	CHECK(patch2->starts[1] == Approx(0));
	CHECK(patch2->ns[0] == 10);
	CHECK(patch2->ns[1] == 10);

	CHECK_FALSE(dg.hasCoarserDomain());
}