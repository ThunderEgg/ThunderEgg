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
#include <ThunderEgg/PETSc/VecWrapper.h>
using namespace std;
using namespace ThunderEgg;
#define MESHES                                                                                     \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";
TEST_CASE("PETSc::VecWrapper<1> getNumGhostCells", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<2>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getNumGhostCells() == num_ghost_cells);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<1> getMPIComm", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);
	int           size              = (nx + 2 * num_ghost_cells) * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<1>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getMPIComm() == MPI_COMM_WORLD);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<1> getNumLocalPatches", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);
	int           size              = (nx + 2 * num_ghost_cells) * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<1>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getNumLocalPatches() == num_local_patches);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<1> getNumLocalCells", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);
	int           size              = (nx + 2 * num_ghost_cells) * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<1>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getNumLocalCells() == nx * num_local_patches);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<2> getNumGhostCells", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<2>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getNumGhostCells() == num_ghost_cells);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<2> getMPIComm", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<2>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getMPIComm() == MPI_COMM_WORLD);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<2> getNumLocalPatches", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<2>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getNumLocalPatches() == num_local_patches);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<2> getNumLocalCells", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<2>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getNumLocalCells() == nx * ny * num_local_patches);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<2> ", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells = GENERATE(0, 1, 5);
	int           nx              = GENERATE(1, 4, 5);
	int           ny              = GENERATE(1, 4, 5);
	array<int, 2> ns              = {nx, ny};
	int           size            = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<2>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getNumGhostCells() == num_ghost_cells);

	VecDestroy(&vec);
}

TEST_CASE("PETSc::VecWrapper<3> getNumGhostCells", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells)
	           * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<3>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getNumGhostCells() == num_ghost_cells);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<3> getMPIComm", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells)
	           * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<3>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getMPIComm() == MPI_COMM_WORLD);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<3> getNumLocalPatches", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells)
	           * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<3>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getNumLocalPatches() == num_local_patches);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<3> getNumLocalCells", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells)
	           * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<3>>(vec, ns, num_ghost_cells, false);

	CHECK(vec_wrapper->getNumLocalCells() == nx * ny * nz * num_local_patches);

	VecDestroy(&vec);
}

TEST_CASE("PETSc::VecWrapper<1> getLocalData", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);
	int           patch_stride      = (nx + 2 * num_ghost_cells);
	int           size              = (nx + 2 * num_ghost_cells) * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<1>>(vec, ns, num_ghost_cells, false);

	double *view;
	VecGetArray(vec, &view);
	for (int i = 0; i < num_local_patches; i++) {
		INFO("i:                 " << i);
		LocalData<1> ld = vec_wrapper->getLocalData(i);
		CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i);
		CHECK(&ld[ld.getGhostEnd()] == view + patch_stride * (i + 1) - 1);
	}
	VecRestoreArray(vec, &view);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<1> getLocalData const", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	array<int, 1> ns                = {nx};
	int           num_local_patches = GENERATE(1, 13);
	int           patch_stride      = (nx + 2 * num_ghost_cells);
	int           size              = (nx + 2 * num_ghost_cells) * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<1>>(vec, ns, num_ghost_cells, false);

	double *view;
	VecGetArray(vec, &view);
	for (int i = 0; i < num_local_patches; i++) {
		INFO("i:                 " << i);
		const LocalData<1> ld = vec_wrapper->getLocalData(i);
		CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i);
		CHECK(&ld[ld.getGhostEnd()] == view + patch_stride * (i + 1) - 1);
	}
	VecRestoreArray(vec, &view);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<2> getLocalData", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);
	int           patch_stride      = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<2>>(vec, ns, num_ghost_cells, false);

	double *view;
	VecGetArray(vec, &view);
	for (int i = 0; i < num_local_patches; i++) {
		INFO("i:                 " << i);
		LocalData<2> ld = vec_wrapper->getLocalData(i);
		CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i);
		CHECK(&ld[ld.getGhostEnd()] == view + patch_stride * (i + 1) - 1);
	}
	VecRestoreArray(vec, &view);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<2> getLocalData const", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	array<int, 2> ns                = {nx, ny};
	int           num_local_patches = GENERATE(1, 13);
	int           patch_stride      = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<2>>(vec, ns, num_ghost_cells, false);

	double *view;
	VecGetArray(vec, &view);
	for (int i = 0; i < num_local_patches; i++) {
		INFO("i:                 " << i);
		const LocalData<2> ld = vec_wrapper->getLocalData(i);
		CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i);
		CHECK(&ld[ld.getGhostEnd()] == view + patch_stride * (i + 1) - 1);
	}
	VecRestoreArray(vec, &view);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<3> getLocalData", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);
	int           patch_stride
	= (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells)
	           * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<3>>(vec, ns, num_ghost_cells, false);

	double *view;
	VecGetArray(vec, &view);
	for (int i = 0; i < num_local_patches; i++) {
		INFO("i:                 " << i);
		LocalData<3> ld = vec_wrapper->getLocalData(i);
		CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i);
		CHECK(&ld[ld.getGhostEnd()] == view + patch_stride * (i + 1) - 1);
	}
	VecRestoreArray(vec, &view);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper<3> getLocalData const", "[PETSc::VecWrapper]")
{
	auto          num_ghost_cells   = GENERATE(0, 1, 5);
	int           nx                = GENERATE(1, 4, 5);
	int           ny                = GENERATE(1, 4, 5);
	int           nz                = GENERATE(1, 4, 5);
	array<int, 3> ns                = {nx, ny, nz};
	int           num_local_patches = GENERATE(1, 13);
	int           patch_stride
	= (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells);
	int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells)
	           * num_local_patches;

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);
	INFO("nz:                " << nz);
	INFO("num_local_patches: " << num_local_patches);

	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, size, PETSC_DETERMINE, &vec);

	auto vec_wrapper = make_shared<PETSc::VecWrapper<3>>(vec, ns, num_ghost_cells, false);

	double *view;
	VecGetArray(vec, &view);
	for (int i = 0; i < num_local_patches; i++) {
		INFO("i:                 " << i);
		const LocalData<3> ld = vec_wrapper->getLocalData(i);
		CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i);
		CHECK(&ld[ld.getGhostEnd()] == view + patch_stride * (i + 1) - 1);
	}
	VecRestoreArray(vec, &view);

	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecWrapper getNewVector works", "[PETSc::VecWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int nx = GENERATE(1, 4, 5);
	int ny = GENERATE(1, 4, 5);
	INFO("nx:       " << nx);
	INFO("ny:       " << ny);
	int                   num_ghost = 0;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec_wrapper = PETSc::VecWrapper<2>::GetNewVector(d_fine);

	CHECK(vec_wrapper->getNumGhostCells() == 0);
	CHECK(vec_wrapper->getNumLocalCells() == d_fine->getNumLocalCells());
	CHECK(vec_wrapper->getNumLocalPatches() == d_fine->getNumLocalPatches());
	CHECK(vec_wrapper->getMPIComm() == MPI_COMM_WORLD);
	CHECK(vec_wrapper->getLocalData(0).getLengths()[0] == nx);
	CHECK(vec_wrapper->getLocalData(0).getLengths()[1] == ny);
}