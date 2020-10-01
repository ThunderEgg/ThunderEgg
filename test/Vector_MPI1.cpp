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
#define MESHES                                                                                     \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";
TEST_CASE("MockVector<3> getMPIComm", "[MockVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	int           num_ghost_cells   = GENERATE(0, 1, 5);
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

	CHECK(vec.getMPIComm() == MPI_COMM_WORLD);
}
TEST_CASE("Vector<3> getNumComponents", "[Vector]")
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

	CHECK(vec.getNumComponents() == num_components);
}
TEST_CASE("Vector<3> getNumLocalPatches", "[Vector]")
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

	CHECK(vec.getNumLocalPatches() == num_local_patches);
}
TEST_CASE("Vector<3> getNumLocalCells", "[Vector]")
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

	CHECK(vec.getNumLocalCells() == nx * ny * nz * num_local_patches * num_components);
}
TEST_CASE("Vector<3> getLocalDatas", "[Vector]")
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

	for (int i = 0; i < vec.getNumLocalPatches(); i++) {
		auto lds = vec.getLocalDatas(i);
		for (int c = 0; c < vec.getNumComponents(); c++) {
			auto ld = vec.getLocalData(c, i);
			CHECK(ld.getPtr() == lds[c].getPtr());
		}
	}
}
TEST_CASE("Vector<3> getLocalDatas const", "[Vector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	const MockVector<3> vec(MPI_COMM_WORLD, num_components, num_local_patches, num_ghost_cells, ns);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);
	INFO("num_components:    " << num_components);

	for (int i = 0; i < vec.getNumLocalPatches(); i++) {
		auto lds = vec.getLocalDatas(i);
		for (int c = 0; c < vec.getNumComponents(); c++) {
			auto ld = vec.getLocalData(c, i);
			CHECK(ld.getPtr() == lds[c].getPtr());
		}
	}
}
TEST_CASE("Vector<3> set", "[Vector]")
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

	vec.set(28);

	for (int i = 0; i < vec.getNumLocalPatches(); i++) {
		for (int c = 0; c < vec.getNumComponents(); c++) {
			auto ld = vec.getLocalData(c, i);
			nested_loop<3>(ld.getGhostStart(), ld.getGhostEnd(), [&](std::array<int, 3> &coord) {
				if (isGhost(coord, ns, num_ghost_cells)) {
					CHECK(ld[coord] == 0);
				} else {
					CHECK(ld[coord] == 28);
				}
			});
		}
	}
}
TEST_CASE("Vector<3> setWithGhost", "[Vector]")
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

	vec.setWithGhost(28);

	for (int i = 0; i < vec.getNumLocalPatches(); i++) {
		for (int c = 0; c < vec.getNumComponents(); c++) {
			auto ld = vec.getLocalData(c, i);
			nested_loop<3>(ld.getGhostStart(), ld.getGhostEnd(),
			               [&](std::array<int, 3> &coord) { CHECK(ld[coord] == 28); });
		}
	}
}
TEST_CASE("Vector<3> scale", "[Vector]")
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

	vec.setWithGhost(28);
	vec.scale(0.25);

	for (int i = 0; i < vec.getNumLocalPatches(); i++) {
		for (int c = 0; c < vec.getNumComponents(); c++) {
			auto ld = vec.getLocalData(c, i);
			nested_loop<3>(ld.getGhostStart(), ld.getGhostEnd(), [&](std::array<int, 3> &coord) {
				if (isGhost(coord, ns, num_ghost_cells)) {
					CHECK(ld[coord] == 28);
				} else {
					CHECK(ld[coord] == Approx(7));
				}
			});
		}
	}
}
TEST_CASE("Vector<3> shift", "[Vector]")
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

	vec.setWithGhost(1);
	vec.shift(28);

	for (int i = 0; i < vec.getNumLocalPatches(); i++) {
		for (int c = 0; c < vec.getNumComponents(); c++) {
			auto ld = vec.getLocalData(c, i);
			nested_loop<3>(ld.getGhostStart(), ld.getGhostEnd(), [&](std::array<int, 3> &coord) {
				if (isGhost(coord, ns, num_ghost_cells)) {
					CHECK(ld[coord] == 1);
				} else {
					CHECK(ld[coord] == Approx(29));
				}
			});
		}
	}
}
TEST_CASE("Vector<3> copy", "[Vector]")
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

	b->copy(a);

	for (int i = 0; i < a->getNumLocalPatches(); i++) {
		for (int c = 0; c < a->getNumComponents(); c++) {
			auto a_ld = a->getLocalData(c, i);
			auto b_ld = b->getLocalData(c, i);
			nested_loop<3>(a_ld.getGhostStart(), a_ld.getGhostEnd(),
			               [&](std::array<int, 3> &coord) {
				               if (isGhost(coord, ns, num_ghost_cells)) {
					               CHECK(b_ld[coord] == 0);
				               } else {
					               CHECK(b_ld[coord] == a_ld[coord]);
				               }
			               });
		}
	}
}
TEST_CASE("Vector<3> add", "[Vector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	auto a        = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                        num_ghost_cells, ns);
	auto b        = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                        num_ghost_cells, ns);
	auto b_copy   = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                             num_ghost_cells, ns);
	auto expected = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
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
		double x        = (i + 0.5) / b->data.size();
		b->data[i]      = (x - 0.5) * (x - 0.5);
		b_copy->data[i] = (x - 0.5) * (x - 0.5);
	}

	for (int i = 0; i < expected->data.size(); i++) {
		double x          = (i + 0.5) / expected->data.size();
		expected->data[i] = 10 - (x - 0.75) * (x - 0.75) + (x - 0.5) * (x - 0.5);
	}
	b->add(a);

	for (int i = 0; i < a->getNumLocalPatches(); i++) {
		for (int c = 0; c < a->getNumComponents(); c++) {
			auto expected_ld = expected->getLocalData(c, i);
			auto b_ld        = b->getLocalData(c, i);
			auto b_copy_ld   = b_copy->getLocalData(c, i);
			nested_loop<3>(b_ld.getGhostStart(), b_ld.getGhostEnd(),
			               [&](std::array<int, 3> &coord) {
				               if (isGhost(coord, ns, num_ghost_cells)) {
					               CHECK(b_ld[coord] == b_copy_ld[coord]);
				               } else {
					               CHECK(b_ld[coord] == Approx(expected_ld[coord]));
				               }
			               });
		}
	}
}
TEST_CASE("Vector<3> addScaled", "[Vector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	auto a        = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                        num_ghost_cells, ns);
	auto b        = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                        num_ghost_cells, ns);
	auto b_copy   = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                             num_ghost_cells, ns);
	auto expected = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
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
		double x        = (i + 0.5) / b->data.size();
		b->data[i]      = (x - 0.5) * (x - 0.5);
		b_copy->data[i] = (x - 0.5) * (x - 0.5);
	}

	for (int i = 0; i < expected->data.size(); i++) {
		double x          = (i + 0.5) / expected->data.size();
		expected->data[i] = 0.7 * (10 - (x - 0.75) * (x - 0.75)) + (x - 0.5) * (x - 0.5);
	}
	b->addScaled(0.7, a);

	for (int i = 0; i < a->getNumLocalPatches(); i++) {
		for (int c = 0; c < a->getNumComponents(); c++) {
			auto expected_ld = expected->getLocalData(c, i);
			auto b_ld        = b->getLocalData(c, i);
			auto b_copy_ld   = b_copy->getLocalData(c, i);
			nested_loop<3>(b_ld.getGhostStart(), b_ld.getGhostEnd(),
			               [&](std::array<int, 3> &coord) {
				               if (isGhost(coord, ns, num_ghost_cells)) {
					               CHECK(b_ld[coord] == b_copy_ld[coord]);
				               } else {
					               CHECK(b_ld[coord] == Approx(expected_ld[coord]));
				               }
			               });
		}
	}
}
TEST_CASE("Vector<3> scaleThenAdd", "[Vector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	auto a        = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                        num_ghost_cells, ns);
	auto b        = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                        num_ghost_cells, ns);
	auto b_copy   = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                             num_ghost_cells, ns);
	auto expected = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
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
		double x        = (i + 0.5) / b->data.size();
		b->data[i]      = (x - 0.5) * (x - 0.5);
		b_copy->data[i] = (x - 0.5) * (x - 0.5);
	}

	for (int i = 0; i < expected->data.size(); i++) {
		double x          = (i + 0.5) / expected->data.size();
		expected->data[i] = (10 - (x - 0.75) * (x - 0.75)) + 0.7 * (x - 0.5) * (x - 0.5);
	}
	b->scaleThenAdd(0.7, a);

	for (int i = 0; i < a->getNumLocalPatches(); i++) {
		for (int c = 0; c < a->getNumComponents(); c++) {
			auto expected_ld = expected->getLocalData(c, i);
			auto b_ld        = b->getLocalData(c, i);
			auto b_copy_ld   = b_copy->getLocalData(c, i);
			nested_loop<3>(b_ld.getGhostStart(), b_ld.getGhostEnd(),
			               [&](std::array<int, 3> &coord) {
				               if (isGhost(coord, ns, num_ghost_cells)) {
					               CHECK(b_ld[coord] == b_copy_ld[coord]);
				               } else {
					               CHECK(b_ld[coord] == Approx(expected_ld[coord]));
				               }
			               });
		}
	}
}
TEST_CASE("Vector<3> scaleThenAddScaled", "[Vector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	auto a        = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                        num_ghost_cells, ns);
	auto b        = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                        num_ghost_cells, ns);
	auto b_copy   = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                             num_ghost_cells, ns);
	auto expected = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
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
		double x        = (i + 0.5) / b->data.size();
		b->data[i]      = (x - 0.5) * (x - 0.5);
		b_copy->data[i] = (x - 0.5) * (x - 0.5);
	}

	for (int i = 0; i < expected->data.size(); i++) {
		double x          = (i + 0.5) / expected->data.size();
		expected->data[i] = -2 * (10 - (x - 0.75) * (x - 0.75)) + 0.7 * (x - 0.5) * (x - 0.5);
	}
	b->scaleThenAddScaled(0.7, -2, a);

	for (int i = 0; i < a->getNumLocalPatches(); i++) {
		for (int c = 0; c < a->getNumComponents(); c++) {
			auto expected_ld = expected->getLocalData(c, i);
			auto b_ld        = b->getLocalData(c, i);
			auto b_copy_ld   = b_copy->getLocalData(c, i);
			nested_loop<3>(b_ld.getGhostStart(), b_ld.getGhostEnd(),
			               [&](std::array<int, 3> &coord) {
				               if (isGhost(coord, ns, num_ghost_cells)) {
					               CHECK(b_ld[coord] == b_copy_ld[coord]);
				               } else {
					               CHECK(b_ld[coord] == Approx(expected_ld[coord]));
				               }
			               });
		}
	}
}
TEST_CASE("Vector<3> scaleThenAddScaled two vectors", "[Vector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	auto a        = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                        num_ghost_cells, ns);
	auto b        = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                        num_ghost_cells, ns);
	auto c        = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                        num_ghost_cells, ns);
	auto b_copy   = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
                                             num_ghost_cells, ns);
	auto expected = make_shared<MockVector<3>>(MPI_COMM_WORLD, num_components, num_local_patches,
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
		double x        = (i + 0.5) / b->data.size();
		b->data[i]      = (x - 0.5) * (x - 0.5);
		b_copy->data[i] = (x - 0.5) * (x - 0.5);
	}

	for (int i = 0; i < c->data.size(); i++) {
		double x   = (i + 0.5) / c->data.size();
		c->data[i] = 1 + (x - 0.25) * (x - 0.25);
	}

	for (int i = 0; i < expected->data.size(); i++) {
		double x          = (i + 0.5) / expected->data.size();
		expected->data[i] = -2 * (10 - (x - 0.75) * (x - 0.75)) + 0.7 * (x - 0.5) * (x - 0.5)
		                    + 9 * (1 + (x - 0.25) * (x - 0.25));
	}
	b->scaleThenAddScaled(0.7, -2, a, 9, c);

	for (int i = 0; i < a->getNumLocalPatches(); i++) {
		for (int c = 0; c < a->getNumComponents(); c++) {
			auto expected_ld = expected->getLocalData(c, i);
			auto b_ld        = b->getLocalData(c, i);
			auto b_copy_ld   = b_copy->getLocalData(c, i);
			nested_loop<3>(b_ld.getGhostStart(), b_ld.getGhostEnd(),
			               [&](std::array<int, 3> &coord) {
				               if (isGhost(coord, ns, num_ghost_cells)) {
					               CHECK(b_ld[coord] == b_copy_ld[coord]);
				               } else {
					               CHECK(b_ld[coord] == Approx(expected_ld[coord]));
				               }
			               });
		}
	}
}
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
	expected_norm = sqrt(expected_norm);

	CHECK(vec.twoNorm() == Approx(expected_norm));
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

	CHECK(vec.infNorm() == expected_norm);
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

	CHECK(a->dot(b) == Approx(expected_value));
}