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
#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/Experimental/DomGen.h>
#include <ThunderEgg/Experimental/OctTree.h>
#include <ThunderEgg/GMG/LinearRestrictor.h>
#include <ThunderEgg/ValVector.h>
using namespace std;
using namespace ThunderEgg;
const string uniform_mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";
const string refined_mesh_file = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json";
TEST_CASE("Linear Test LinearRestrictor", "[GMG::LinearRestrictor]")
{
	auto mesh_file = GENERATE(as<std::string>{}, uniform_mesh_file, refined_mesh_file);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();

	auto fine_vec        = ValVector<2>::GetNewVector(d_fine, 1);
	auto coarse_vec      = ValVector<2>::GetNewVector(d_coarse, 1);
	auto coarse_expected = ValVector<2>::GetNewVector(d_coarse, 1);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 1 + ((x * 0.3) + y);
	};

	DomainTools::SetValuesWithGhost<2>(d_fine, fine_vec, f);
	DomainTools::SetValuesWithGhost<2>(d_coarse, coarse_expected, f);

	auto restrictor = std::make_shared<GMG::LinearRestrictor<2>>(d_fine, d_coarse, 1, true);

	restrictor->restrict(fine_vec, coarse_vec);

	for (auto pinfo : d_coarse->getPatchInfoVector()) {
		INFO("Patch:          " << pinfo->id);
		INFO("x:              " << pinfo->starts[0]);
		INFO("y:              " << pinfo->starts[1]);
		INFO("nx:             " << pinfo->ns[0]);
		INFO("ny:             " << pinfo->ns[1]);
		INFO("parent_orth:    " << pinfo->orth_on_parent);
		LocalData<2> vec_ld      = coarse_vec->getLocalData(0, pinfo->local_index);
		LocalData<2> expected_ld = coarse_expected->getLocalData(0, pinfo->local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			REQUIRE(vec_ld[coord] == Approx(expected_ld[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			INFO("side:      " << s);
			if (!pinfo->hasNbr(s)) {
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 1> &coord) {
					               INFO("coord:  " << coord[0]);
					               CHECK(vec_ghost[coord] == Approx(expected_ghost[coord]));
				               });
			} else {
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 1> &coord) {
					               INFO("coord:  " << coord[0]);
					               CHECK(vec_ghost[coord] == 0);
				               });
			}
		}
	}
}
TEST_CASE("Linear Test LinearRestrictor two components", "[GMG::LinearRestrictor]")
{
	auto mesh_file = GENERATE(as<std::string>{}, uniform_mesh_file, refined_mesh_file);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();

	auto fine_vec        = ValVector<2>::GetNewVector(d_fine, 2);
	auto coarse_vec      = ValVector<2>::GetNewVector(d_coarse, 2);
	auto coarse_expected = ValVector<2>::GetNewVector(d_coarse, 2);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 1 + ((x * 0.3) + y);
	};
	auto g = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 9 + ((x * 0.9) + y * 4);
	};

	DomainTools::SetValuesWithGhost<2>(d_fine, fine_vec, f, g);
	DomainTools::SetValuesWithGhost<2>(d_coarse, coarse_expected, f, g);

	auto restrictor = std::make_shared<GMG::LinearRestrictor<2>>(d_fine, d_coarse, 2, true);

	restrictor->restrict(fine_vec, coarse_vec);

	for (auto pinfo : d_coarse->getPatchInfoVector()) {
		INFO("Patch:          " << pinfo->id);
		INFO("x:              " << pinfo->starts[0]);
		INFO("y:              " << pinfo->starts[1]);
		INFO("nx:             " << pinfo->ns[0]);
		INFO("ny:             " << pinfo->ns[1]);
		INFO("parent_orth:    " << pinfo->orth_on_parent);
		LocalData<2> vec_ld       = coarse_vec->getLocalData(0, pinfo->local_index);
		LocalData<2> expected_ld  = coarse_expected->getLocalData(0, pinfo->local_index);
		LocalData<2> vec_ld2      = coarse_vec->getLocalData(1, pinfo->local_index);
		LocalData<2> expected_ld2 = coarse_expected->getLocalData(1, pinfo->local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			REQUIRE(vec_ld[coord] == Approx(expected_ld[coord]));
			REQUIRE(vec_ld2[coord] == Approx(expected_ld2[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost       = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost  = expected_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> vec_ghost2      = vec_ld2.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost2 = expected_ld2.getGhostSliceOnSide(s, 1);
			INFO("side:      " << s);
			if (!pinfo->hasNbr(s)) {
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 1> &coord) {
					               INFO("coord:  " << coord[0]);
					               CHECK(vec_ghost[coord] == Approx(expected_ghost[coord]));
					               CHECK(vec_ghost2[coord] == Approx(expected_ghost2[coord]));
				               });
			} else {
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 1> &coord) {
					               INFO("coord:  " << coord[0]);
					               CHECK(vec_ghost[coord] == 0);
					               CHECK(vec_ghost2[coord] == 0);
				               });
			}
		}
	}
}
TEST_CASE("Linear Test LinearRestrictor dont extrapolate bound ghosts", "[GMG::LinearRestrictor]")
{
	auto mesh_file = GENERATE(as<std::string>{}, uniform_mesh_file, refined_mesh_file);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();

	auto fine_vec        = ValVector<2>::GetNewVector(d_fine, 1);
	auto coarse_vec      = ValVector<2>::GetNewVector(d_coarse, 1);
	auto coarse_expected = ValVector<2>::GetNewVector(d_coarse, 1);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 1 + ((x * 0.3) + y);
	};

	DomainTools::SetValuesWithGhost<2>(d_fine, fine_vec, f);
	DomainTools::SetValuesWithGhost<2>(d_coarse, coarse_expected, f);

	auto restrictor = std::make_shared<GMG::LinearRestrictor<2>>(d_fine, d_coarse, 1, false);

	restrictor->restrict(fine_vec, coarse_vec);

	for (auto pinfo : d_coarse->getPatchInfoVector()) {
		INFO("Patch:          " << pinfo->id);
		INFO("x:              " << pinfo->starts[0]);
		INFO("y:              " << pinfo->starts[1]);
		INFO("nx:             " << pinfo->ns[0]);
		INFO("ny:             " << pinfo->ns[1]);
		INFO("parent_orth:    " << pinfo->orth_on_parent);
		LocalData<2> vec_ld      = coarse_vec->getLocalData(0, pinfo->local_index);
		LocalData<2> expected_ld = coarse_expected->getLocalData(0, pinfo->local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			REQUIRE(vec_ld[coord] == Approx(expected_ld[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			INFO("side:      " << s);
			nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
			               [&](const array<int, 1> &coord) {
				               INFO("coord:  " << coord[0]);
				               CHECK(vec_ghost[coord] == 0);
			               });
		}
	}
}
TEST_CASE("Linear Test LinearRestrictor two components dont extrapolate boundary ghosts",
          "[GMG::LinearRestrictor]")
{
	auto mesh_file = GENERATE(as<std::string>{}, uniform_mesh_file, refined_mesh_file);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();

	auto fine_vec        = ValVector<2>::GetNewVector(d_fine, 2);
	auto coarse_vec      = ValVector<2>::GetNewVector(d_coarse, 2);
	auto coarse_expected = ValVector<2>::GetNewVector(d_coarse, 2);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 1 + ((x * 0.3) + y);
	};
	auto g = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 9 + ((x * 0.9) + y * 4);
	};

	DomainTools::SetValuesWithGhost<2>(d_fine, fine_vec, f, g);
	DomainTools::SetValuesWithGhost<2>(d_coarse, coarse_expected, f, g);

	auto restrictor = std::make_shared<GMG::LinearRestrictor<2>>(d_fine, d_coarse, 2, false);

	restrictor->restrict(fine_vec, coarse_vec);

	for (auto pinfo : d_coarse->getPatchInfoVector()) {
		INFO("Patch:          " << pinfo->id);
		INFO("x:              " << pinfo->starts[0]);
		INFO("y:              " << pinfo->starts[1]);
		INFO("nx:             " << pinfo->ns[0]);
		INFO("ny:             " << pinfo->ns[1]);
		INFO("parent_orth:    " << pinfo->orth_on_parent);
		LocalData<2> vec_ld       = coarse_vec->getLocalData(0, pinfo->local_index);
		LocalData<2> expected_ld  = coarse_expected->getLocalData(0, pinfo->local_index);
		LocalData<2> vec_ld2      = coarse_vec->getLocalData(1, pinfo->local_index);
		LocalData<2> expected_ld2 = coarse_expected->getLocalData(1, pinfo->local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			REQUIRE(vec_ld[coord] == Approx(expected_ld[coord]));
			REQUIRE(vec_ld2[coord] == Approx(expected_ld2[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost       = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost  = expected_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> vec_ghost2      = vec_ld2.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost2 = expected_ld2.getGhostSliceOnSide(s, 1);
			INFO("side:      " << s);
			nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
			               [&](const array<int, 1> &coord) {
				               INFO("coord:  " << coord[0]);
				               CHECK(vec_ghost[coord] == 0);
				               CHECK(vec_ghost2[coord] == 0);
			               });
		}
	}
}