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
#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/DirectInterpolator.h>
#include <ThunderEgg/ValVector.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

const string mesh_file         = "mesh_inputs/2d_uniform_4x4_mpi1.json";
const string refined_mesh_file = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json";
TEST_CASE("Test DirectInterpolator on uniform 4x4", "[GMG::DirectInterpolator]")
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
		PatchView<double, 2> view = coarse_vec->getPatchView(pinfo.local_index);
		loop_over_interior_indexes<3>(view, [&](const array<int, 3> &coord) {
			view[coord] = 1 + pinfo.id * nx * ny + coord[0] + coord[1] * nx + coord[2];
		});
	}

	// set expected finer vector vector
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		PatchView<double, 2> view = fine_expected->getPatchView(pinfo.local_index);

		Orthant<2>         orth = pinfo.orth_on_parent;
		std::array<int, 2> starts;
		for (size_t i = 0; i < 2; i++) {
			starts[i] = orth.isOnSide(Side<2>(2 * i)) ? 0 : (view.getEnd()[i] + 1);
		}

		loop_over_interior_indexes<3>(view, [&](const array<int, 3> &coord) {
			view[coord] = 1 + pinfo.parent_id * nx * ny + (coord[0] + starts[0]) / 2
			              + (coord[1] + starts[1]) / 2 * nx + coord[2];
		});
	}

	auto interpolator
	= std::make_shared<GMG::DirectInterpolator<2>>(d_coarse, d_fine, num_components);

	interpolator->interpolate(*coarse_vec, *fine_vec);

	for (auto pinfo : d_fine->getPatchInfoVector()) {
		INFO("Patch: " << pinfo.id);
		INFO("x:     " << pinfo.starts[0]);
		INFO("y:     " << pinfo.starts[1]);
		INFO("nx:    " << pinfo.ns[0]);
		INFO("c:     " << pinfo.ns[1]);
		PatchView<double, 2> vec_view      = fine_vec->getPatchView(pinfo.local_index);
		PatchView<double, 2> expected_view = fine_expected->getPatchView(pinfo.local_index);
		loop_over_interior_indexes<3>(vec_view,
		                              [&](const array<int, 3> &coord) {
			                              REQUIRE(vec_view[coord] == Catch::Approx(expected_view[coord]));
		                              });
	}
}
TEST_CASE("Linear Test DirectInterpolator with values already set on uniform 4x4",
          "[GMG::DirectInterpolator]")
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
		PatchView<double, 2> view = coarse_vec->getPatchView(pinfo.local_index);
		loop_over_interior_indexes<3>(view, [&](const array<int, 3> &coord) {
			view[coord] = 1 + pinfo.id * nx * ny + coord[0] + coord[1] * nx + coord[2];
		});
	}

	// set expected finer vector vector
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		PatchView<double, 2> view = fine_expected->getPatchView(pinfo.local_index);

		Orthant<2>         orth = pinfo.orth_on_parent;
		std::array<int, 2> starts;
		for (size_t i = 0; i < 2; i++) {
			starts[i] = orth.isOnSide(Side<2>(2 * i)) ? 0 : (view.getEnd()[i] + 1);
		}

		loop_over_interior_indexes<3>(view, [&](const array<int, 3> &coord) {
			view[coord] = 2 + pinfo.parent_id * nx * ny + (coord[0] + starts[0]) / 2
			              + (coord[1] + starts[1]) / 2 * nx + coord[2];
		});
	}

	fine_vec->set(1.0);

	auto interpolator
	= std::make_shared<GMG::DirectInterpolator<2>>(d_coarse, d_fine, num_components);

	interpolator->interpolate(*coarse_vec, *fine_vec);

	for (auto pinfo : d_fine->getPatchInfoVector()) {
		INFO("Patch: " << pinfo.id);
		INFO("x:     " << pinfo.starts[0]);
		INFO("y:     " << pinfo.starts[1]);
		INFO("nx:    " << pinfo.ns[0]);
		INFO("ny:    " << pinfo.ns[1]);
		PatchView<double, 2> vec_view      = fine_vec->getPatchView(pinfo.local_index);
		PatchView<double, 2> expected_view = fine_expected->getPatchView(pinfo.local_index);
		loop_over_interior_indexes<3>(vec_view,
		                              [&](const array<int, 3> &coord) {
			                              REQUIRE(vec_view[coord] == Catch::Approx(expected_view[coord]));
		                              });
	}
}
TEST_CASE("Test DirectInterpolator on refined 2x2", "[GMG::DirectInterpolator]")
{
	auto                  num_components = GENERATE(1, 2, 3);
	auto                  nx             = GENERATE(2, 10);
	auto                  ny             = GENERATE(2, 10);
	int                   num_ghost      = 1;
	DomainReader<2>       domain_reader(refined_mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();

	auto coarse_vec    = ValVector<2>::GetNewVector(d_coarse, num_components);
	auto fine_vec      = ValVector<2>::GetNewVector(d_fine, num_components);
	auto fine_expected = ValVector<2>::GetNewVector(d_fine, num_components);

	// set coarse vector
	for (auto pinfo : d_coarse->getPatchInfoVector()) {
		PatchView<double, 2> view = coarse_vec->getPatchView(pinfo.local_index);
		loop_over_interior_indexes<3>(view, [&](const array<int, 3> &coord) {
			view[coord] = 1 + pinfo.id * nx * ny + coord[0] + coord[1] * nx + coord[2];
		});
	}

	// set expected finer vector vector
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		PatchView<double, 2> view = fine_expected->getPatchView(pinfo.local_index);

		if (pinfo.hasCoarseParent()) {
			Orthant<2>         orth = pinfo.orth_on_parent;
			std::array<int, 2> starts;
			for (size_t i = 0; i < 2; i++) {
				starts[i] = orth.isOnSide(Side<2>(2 * i)) ? 0 : (view.getEnd()[i] + 1);
			}

			loop_over_interior_indexes<3>(view, [&](const array<int, 3> &coord) {
				view[coord] = 1 + pinfo.parent_id * nx * ny + (coord[0] + starts[0]) / 2
				              + (coord[1] + starts[1]) / 2 * nx + coord[2];
			});
		} else {
			loop_over_interior_indexes<3>(view, [&](const array<int, 3> &coord) {
				view[coord] = 1 + pinfo.id * nx * ny + coord[0] + coord[1] * nx + coord[2];
			});
		}
	}

	auto interpolator
	= std::make_shared<GMG::DirectInterpolator<2>>(d_coarse, d_fine, num_components);

	interpolator->interpolate(*coarse_vec, *fine_vec);

	for (auto pinfo : d_fine->getPatchInfoVector()) {
		INFO("Patch: " << pinfo.id);
		INFO("x:     " << pinfo.starts[0]);
		INFO("y:     " << pinfo.starts[1]);
		INFO("nx:    " << pinfo.ns[0]);
		INFO("ny:    " << pinfo.ns[1]);
		PatchView<double, 2> vec_view      = fine_vec->getPatchView(pinfo.local_index);
		PatchView<double, 2> expected_view = fine_expected->getPatchView(pinfo.local_index);
		loop_over_interior_indexes<3>(vec_view,
		                              [&](const array<int, 3> &coord) {
			                              REQUIRE(vec_view[coord] == Catch::Approx(expected_view[coord]));
		                              });
	}
}