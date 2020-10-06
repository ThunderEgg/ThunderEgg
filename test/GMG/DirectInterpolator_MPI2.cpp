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
#include <ThunderEgg/GMG/DirectInterpolator.h>
#include <ThunderEgg/ValVector.h>
using namespace std;
using namespace ThunderEgg;
const string mesh_file = "mesh_inputs/2d_uniform_quad_mpi2.json";
TEST_CASE("Test DirectInterpolator", "[GMG::DirectInterpolator]")
{
	auto                  num_components = GENERATE(1, 2, 3);
	auto                  nx             = GENERATE(2, 10);
	auto                  ny             = GENERATE(2, 10);
	int                   num_ghost      = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();

	auto coarse_vec    = ValVector<2>::GetNewVector(d_coarse, num_components);
	auto fine_vec      = ValVector<2>::GetNewVector(d_fine, num_components);
	auto fine_expected = ValVector<2>::GetNewVector(d_fine, num_components);

	// set coarse vector
	for (auto pinfo : d_coarse->getPatchInfoVector()) {
		auto lds = coarse_vec->getLocalDatas(pinfo->local_index);
		for (int c = 0; c < num_components; c++) {
			nested_loop<2>(lds[c].getStart(), lds[c].getEnd(), [&](const array<int, 2> &coord) {
				lds[c][coord] = 1 + pinfo->id * nx * ny + coord[0] + coord[1] * nx + c;
			});
		}
	}

	// set expected finer vector vector
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		auto lds = fine_expected->getLocalDatas(pinfo->local_index);

		Orthant<2>         orth = pinfo->orth_on_parent;
		std::array<int, 2> starts;
		for (size_t i = 0; i < 2; i++) {
			starts[i] = orth.isOnSide(Side<2>(2 * i)) ? 0 : lds[0].getLengths()[i];
		}

		for (int c = 0; c < num_components; c++) {
			nested_loop<2>(lds[c].getStart(), lds[c].getEnd(), [&](const array<int, 2> &coord) {
				lds[c][coord] = 1 + pinfo->parent_id * nx * ny + (coord[0] + starts[0]) / 2
				                + (coord[1] + starts[1]) / 2 * nx + c;
			});
		}
	}

	auto interpolator
	= std::make_shared<GMG::DirectInterpolator<2>>(d_coarse, d_fine, num_components);

	interpolator->interpolate(coarse_vec, fine_vec);

	for (auto pinfo : d_fine->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		auto vec_lds      = fine_vec->getLocalDatas(pinfo->local_index);
		auto expected_lds = fine_expected->getLocalDatas(pinfo->local_index);
		for (int c = 0; c < num_components; c++) {
			nested_loop<2>(vec_lds[c].getStart(), vec_lds[c].getEnd(),
			               [&](const array<int, 2> &coord) {
				               REQUIRE(vec_lds[c][coord] == Approx(expected_lds[c][coord]));
			               });
		}
	}
}
TEST_CASE("Linear Test DirectInterpolator with values already set", "[GMG::DirectInterpolator]")
{
	auto                  num_components = GENERATE(1, 2, 3);
	auto                  nx             = GENERATE(2, 10);
	auto                  ny             = GENERATE(2, 10);
	int                   num_ghost      = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();

	auto coarse_vec    = ValVector<2>::GetNewVector(d_coarse, num_components);
	auto fine_vec      = ValVector<2>::GetNewVector(d_fine, num_components);
	auto fine_expected = ValVector<2>::GetNewVector(d_fine, num_components);

	// set coarse vector
	for (auto pinfo : d_coarse->getPatchInfoVector()) {
		auto lds = coarse_vec->getLocalDatas(pinfo->local_index);
		for (int c = 0; c < num_components; c++) {
			nested_loop<2>(lds[c].getStart(), lds[c].getEnd(), [&](const array<int, 2> &coord) {
				lds[c][coord] = 1 + pinfo->id * nx * ny + coord[0] + coord[1] * nx + c;
			});
		}
	}

	// set expected finer vector vector
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		auto lds = fine_expected->getLocalDatas(pinfo->local_index);

		Orthant<2>         orth = pinfo->orth_on_parent;
		std::array<int, 2> starts;
		for (size_t i = 0; i < 2; i++) {
			starts[i] = orth.isOnSide(Side<2>(2 * i)) ? 0 : lds[0].getLengths()[i];
		}

		for (int c = 0; c < num_components; c++) {
			nested_loop<2>(lds[c].getStart(), lds[c].getEnd(), [&](const array<int, 2> &coord) {
				lds[c][coord] = 2 + pinfo->parent_id * nx * ny + (coord[0] + starts[0]) / 2
				                + (coord[1] + starts[1]) / 2 * nx + c;
			});
		}
	}

	fine_vec->set(1.0);

	auto interpolator
	= std::make_shared<GMG::DirectInterpolator<2>>(d_coarse, d_fine, num_components);

	interpolator->interpolate(coarse_vec, fine_vec);

	for (auto pinfo : d_fine->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		auto vec_lds      = fine_vec->getLocalDatas(pinfo->local_index);
		auto expected_lds = fine_expected->getLocalDatas(pinfo->local_index);
		for (int c = 0; c < num_components; c++) {
			nested_loop<2>(vec_lds[c].getStart(), vec_lds[c].getEnd(),
			               [&](const array<int, 2> &coord) {
				               REQUIRE(vec_lds[c][coord] == Approx(expected_lds[c][coord]));
			               });
		}
	}
}