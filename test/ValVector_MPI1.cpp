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

#include "utils/DomainReader.h"
#include <ThunderEgg/ValVector.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

#define MESHES \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";

TEST_CASE("ValVector<1> getNumGhostCells", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);
	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<1>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getNumGhostCells() == num_ghost_cells);
}
TEST_CASE("ValVector<1> getMPIComm", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<1>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getMPIComm() == MPI_COMM_WORLD);
}
TEST_CASE("ValVector<1> getNumLocalPatches", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<1>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getNumLocalPatches() == num_local_patches);
}
TEST_CASE("ValVector<1> getNumComponents", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<1>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getNumComponents() == num_components);
}
TEST_CASE("ValVector<1> getNumLocalCells", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<1>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getNumLocalCells() == nx * num_local_patches);
}
TEST_CASE("ValVector<1> getValArray", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);
	size_t        size              = (nx + 2 * num_ghost_cells) * num_local_patches * num_components;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<1>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getValArray().size() == size);
}
TEST_CASE("ValVector<2> getNumGhostCells", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<2>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getNumGhostCells() == num_ghost_cells);
}
TEST_CASE("ValVector<2> getMPIComm", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<2>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getMPIComm() == MPI_COMM_WORLD);
}
TEST_CASE("ValVector<2> getNumLocalPatches", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<2>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getNumLocalPatches() == num_local_patches);
}
TEST_CASE("ValVector<2> getNumComponents", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<2>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getNumComponents() == num_components);
}
TEST_CASE("ValVector<2> getNumLocalCells", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<2>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getNumLocalCells() == nx * ny * num_local_patches);
}
TEST_CASE("ValVector<2> getValArray", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);
	size_t        size
	= (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_local_patches * num_components;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<2>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getValArray().size() == size);
}
TEST_CASE("ValVector<3> getNumGhostCells", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<3>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getNumGhostCells() == num_ghost_cells);
}
TEST_CASE("ValVector<3> getMPIComm", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<3>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getMPIComm() == MPI_COMM_WORLD);
}
TEST_CASE("ValVector<3> getNumLocalPatches", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<3>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getNumLocalPatches() == num_local_patches);
}
TEST_CASE("ValVector<3> getNumComponents", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<3>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getNumComponents() == num_components);
}
TEST_CASE("ValVector<3> getNumLocalCells", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<3>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getNumLocalCells() == nx * ny * nz * num_local_patches);
}
TEST_CASE("ValVector<3> getValArray", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);
	size_t        size              = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells)
	              * (nz + 2 * num_ghost_cells) * num_local_patches * num_components;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<3>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	CHECK(val_vector->getValArray().size() == size);
}

TEST_CASE("ValVector<1> getComponentView.h", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);
	int           component_stride  = (nx + 2 * num_ghost_cells);
	int           patch_stride      = component_stride * num_components;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<1>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	double *view = &val_vector->getValArray()[0];
	for (int i = 0; i < num_local_patches; i++) {
		INFO("i:                 " << i);
		for (int c = 0; c < num_components; c++) {
			INFO("c:                 " << c);
			ComponentView<double, 1> ld = val_vector->getComponentView(c, i);
			CHECK(&ld[ld.getGhostStart()] == view + (patch_stride * i + c * component_stride));
			CHECK(&ld[ld.getGhostEnd()]
			      == view + (patch_stride * i + (c + 1) * component_stride) - 1);
		}
	}
}
TEST_CASE("ValVector<1> getComponentView const", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);
	int           component_stride  = (nx + 2 * num_ghost_cells);
	int           patch_stride      = component_stride * num_components;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector       = make_shared<ValVector<1>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
                                                num_local_patches);
	auto const_val_vector = std::const_pointer_cast<const ValVector<1>>(val_vector);

	double *view = &val_vector->getValArray()[0];
	for (int i = 0; i < num_local_patches; i++) {
		INFO("i:                 " << i);
		for (int c = 0; c < num_components; c++) {
			INFO("c:                 " << c);
			ComponentView<const double, 1> ld = const_val_vector->getComponentView(c, i);
			CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i + c * component_stride);
			CHECK(&ld[ld.getGhostEnd()]
			      == view + patch_stride * i + (c + 1) * component_stride - 1);
		}
	}
}
TEST_CASE("ValVector<2> getComponentView.h", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);
	int           component_stride  = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
	int           patch_stride      = component_stride * num_components;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<2>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	double *view = &val_vector->getValArray()[0];
	for (int i = 0; i < num_local_patches; i++) {
		INFO("i:                 " << i);
		for (int c = 0; c < num_components; c++) {
			INFO("c:                 " << c);
			ComponentView<double, 2> ld = val_vector->getComponentView(c, i);
			CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i + c * component_stride);
			CHECK(&ld[ld.getGhostEnd()]
			      == view + patch_stride * i + (c + 1) * component_stride - 1);
		}
	}
}
TEST_CASE("ValVector<2> getComponentView const", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);
	int           component_stride  = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
	int           patch_stride      = component_stride * num_components;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector       = make_shared<ValVector<2>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
                                                num_local_patches);
	auto const_val_vector = std::const_pointer_cast<const ValVector<2>>(val_vector);

	double *view = &val_vector->getValArray()[0];
	for (int i = 0; i < num_local_patches; i++) {
		INFO("i:                 " << i);
		for (int c = 0; c < num_components; c++) {
			INFO("c:                 " << c);
			ComponentView<const double, 2> ld = const_val_vector->getComponentView(c, i);
			CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i + c * component_stride);
			CHECK(&ld[ld.getGhostEnd()]
			      == view + patch_stride * i + (c + 1) * component_stride - 1);
		}
	}
}
TEST_CASE("ValVector<3> getComponentView.h", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);
	int           component_stride
	= (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells);
	int patch_stride = component_stride * num_components;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector = make_shared<ValVector<3>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
	                                            num_local_patches);

	double *view = &val_vector->getValArray()[0];
	for (int i = 0; i < num_local_patches; i++) {
		INFO("i:                 " << i);
		for (int c = 0; c < num_components; c++) {
			ComponentView<double, 3> ld = val_vector->getComponentView(c, i);
			CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i + c * component_stride);
			CHECK(&ld[ld.getGhostEnd()]
			      == view + patch_stride * i + (c + 1) * component_stride - 1);
		}
	}
}
TEST_CASE("ValVector<3> getComponentView const", "[ValVector]")
{
	int           num_components    = GENERATE(1, 2, 3);
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);
	int           component_stride
	= (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells);
	int patch_stride = component_stride * num_components;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	auto val_vector       = make_shared<ValVector<3>>(MPI_COMM_WORLD, ns, num_ghost_cells, num_components,
                                                num_local_patches);
	auto const_val_vector = std::const_pointer_cast<const ValVector<3>>(val_vector);

	double *view = &val_vector->getValArray()[0];
	for (int i = 0; i < num_local_patches; i++) {
		INFO("i:                 " << i);
		for (int c = 0; c < num_components; c++) {
			INFO("c:                 " << c);
			ComponentView<const double, 3> ld = const_val_vector->getComponentView(c, i);
			CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i + c * component_stride);
			CHECK(&ld[ld.getGhostEnd()]
			      == view + patch_stride * i + (c + 1) * component_stride - 1);
		}
	}
}
TEST_CASE("ValVector getNewVector works", "[ValVector]")
{
	int  num_components = GENERATE(1, 2, 3);
	auto mesh_file      = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int nx = GENERATE(1, 4, 5);
	int ny = GENERATE(1, 4, 5);
	INFO("nx:       " << nx);
	INFO("ny:       " << ny);
	int                   num_ghost = 0;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto val_vector = ValVector<2>::GetNewVector(d_fine, num_components);

	CHECK(val_vector->getNumGhostCells() == 0);
	CHECK(val_vector->getNumLocalCells() == d_fine->getNumLocalCells());
	CHECK(val_vector->getNumComponents() == num_components);
	CHECK(val_vector->getNumLocalPatches() == d_fine->getNumLocalPatches());
	CHECK(val_vector->getMPIComm() == MPI_COMM_WORLD);
	CHECK(val_vector->getComponentView(0, 0).getEnd()[0] + 1 == nx);
	CHECK(val_vector->getComponentView(0, 0).getEnd()[1] + 1 == ny);
}