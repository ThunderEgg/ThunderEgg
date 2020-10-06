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

#include "Vector_MOCKS.h"
#include "catch.hpp"
#include "utils/DomainReader.h"
using namespace std;
using namespace ThunderEgg;
TEST_CASE("Vector<3> twoNorm", "[Vector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	MockVector<3> vec(MPI_COMM_WORLD, num_components, num_local_patches, num_ghost_cells, ns);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);
	INFO("num_components:    " << num_components);

	for (int i = 0; i < vec.data.size(); i++) {
		double x    = (i + 0.5) / vec.data.size();
		vec.data[i] = 10 - (x - 0.75) * (x - 0.75);
	}
	vec.setWithGhost(1);
	vec.shift(28);

	double expected_norm = 0;
	for (int i = 0; i < vec.getNumLocalPatches(); i++) {
		for (int c = 0; c < vec.getNumComponents(); c++) {
			auto ld = vec.getLocalData(c, i);
			nested_loop<3>(ld.getStart(), ld.getEnd(), [&](std::array<int, 3> &coord) {
				expected_norm += ld[coord] * ld[coord];
			});
		}
	}
	double global_expected_norm;
	MPI_Allreduce(&expected_norm, &global_expected_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	global_expected_norm = sqrt(global_expected_norm);

	CHECK(vec.twoNorm() == Approx(global_expected_norm));
}
TEST_CASE("Vector<3> infNorm", "[Vector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	MockVector<3> vec(MPI_COMM_WORLD, num_components, num_local_patches, num_ghost_cells, ns);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);
	INFO("num_components:    " << num_components);

	for (int i = 0; i < vec.data.size(); i++) {
		double x    = (i + 0.5) / vec.data.size();
		vec.data[i] = 10 - (x - 0.75) * (x - 0.75);
	}
	vec.setWithGhost(1);
	vec.shift(28);

	double expected_norm = 0;
	for (int i = 0; i < vec.getNumLocalPatches(); i++) {
		for (int c = 0; c < vec.getNumComponents(); c++) {
			auto ld = vec.getLocalData(c, i);
			nested_loop<3>(ld.getStart(), ld.getEnd(), [&](std::array<int, 3> &coord) {
				expected_norm = max(abs(ld[coord]), expected_norm);
			});
		}
	}
	double global_expected_norm;
	MPI_Allreduce(&expected_norm, &global_expected_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	CHECK(vec.infNorm() == global_expected_norm);
}
TEST_CASE("Vector<3> dot", "[Vector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	auto a = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
	                                    num_ghost_cells, ns);
	auto b = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
	                                    num_ghost_cells, ns);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);
	INFO("num_components:    " << num_components);

	for (int i = 0; i < a->data.size(); i++) {
		double x   = (i + 0.5) / a->data.size();
		a->data[i] = 10 - (x - 0.75) * (x - 0.75);
	}

	for (int i = 0; i < b->data.size(); i++) {
		double x   = (i + 0.5) / b->data.size();
		b->data[i] = (x - 0.5) * (x - 0.5);
	}

	double expected_value = 0;
	for (int i = 0; i < a->getNumLocalPatches(); i++) {
		for (int c = 0; c < a->getNumComponents(); c++) {
			auto a_ld = a->getLocalData(c, i);
			auto b_ld = b->getLocalData(c, i);
			nested_loop<3>(b_ld.getStart(), b_ld.getEnd(), [&](std::array<int, 3> &coord) {
				expected_value += b_ld[coord] * a_ld[coord];
			});
		}
	}
	double global_expected_value;
	MPI_Allreduce(&expected_value, &global_expected_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	CHECK(a->dot(b) == Approx(global_expected_value));
}