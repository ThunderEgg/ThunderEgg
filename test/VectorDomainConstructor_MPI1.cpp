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
#include <ThunderEgg/Vector.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

#define MESHES \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"

TEST_CASE("Vector<2> getNumGhostCells domain constructor", "[Vector]")
{
	Communicator comm(MPI_COMM_WORLD);

	auto                  mesh_file       = GENERATE(as<std::string>{}, MESHES);
	int                   num_components  = GENERATE(1, 2, 3);
	auto                  num_ghost_cells = GENERATE(0, 1, 5);
	int                   nx              = GENERATE(1, 4, 5);
	int                   ny              = GENERATE(1, 4, 5);
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost_cells);
	shared_ptr<Domain<2>> domain = domain_reader.getFinerDomain();

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);

	Vector<2> vec(*domain, num_components);

	CHECK(vec.getNumGhostCells() == num_ghost_cells);
}
TEST_CASE("Vector<2> getMPIComm domain constructor", "[Vector]")
{
	Communicator comm(MPI_COMM_WORLD);

	auto                  mesh_file       = GENERATE(as<std::string>{}, MESHES);
	int                   num_components  = GENERATE(1, 2, 3);
	auto                  num_ghost_cells = GENERATE(0, 1, 5);
	int                   nx              = GENERATE(1, 4, 5);
	int                   ny              = GENERATE(1, 4, 5);
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost_cells);
	shared_ptr<Domain<2>> domain = domain_reader.getFinerDomain();

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);

	Vector<2> vec(*domain, num_components);

	int result;
	int err = MPI_Comm_compare(vec.getCommunicator().getMPIComm(), comm.getMPIComm(), &result);
	REQUIRE(err == MPI_SUCCESS);
	CHECK(result == MPI_CONGRUENT);
}
TEST_CASE("Vector<2> getNumLocalPatches domain constructor", "[Vector]")
{
	Communicator comm(MPI_COMM_WORLD);

	auto                  mesh_file       = GENERATE(as<std::string>{}, MESHES);
	int                   num_components  = GENERATE(1, 2, 3);
	auto                  num_ghost_cells = GENERATE(0, 1, 5);
	int                   nx              = GENERATE(1, 4, 5);
	int                   ny              = GENERATE(1, 4, 5);
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost_cells);
	shared_ptr<Domain<2>> domain = domain_reader.getFinerDomain();

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);

	Vector<2> vec(*domain, num_components);

	CHECK(vec.getNumLocalPatches() == domain->getNumLocalPatches());
}
TEST_CASE("Vector<2> getNumComponents domain constructor", "[Vector]")
{
	Communicator comm(MPI_COMM_WORLD);

	auto                  mesh_file       = GENERATE(as<std::string>{}, MESHES);
	int                   num_components  = GENERATE(1, 2, 3);
	auto                  num_ghost_cells = GENERATE(0, 1, 5);
	int                   nx              = GENERATE(1, 4, 5);
	int                   ny              = GENERATE(1, 4, 5);
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost_cells);
	shared_ptr<Domain<2>> domain = domain_reader.getFinerDomain();

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);

	Vector<2> vec(*domain, num_components);

	CHECK(vec.getNumComponents() == num_components);
}
TEST_CASE("Vector<2> getNumLocalCells domain constructor", "[Vector]")
{
	Communicator comm(MPI_COMM_WORLD);

	auto                  mesh_file       = GENERATE(as<std::string>{}, MESHES);
	int                   num_components  = GENERATE(1, 2, 3);
	auto                  num_ghost_cells = GENERATE(0, 1, 5);
	int                   nx              = GENERATE(1, 4, 5);
	int                   ny              = GENERATE(1, 4, 5);
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost_cells);
	shared_ptr<Domain<2>> domain = domain_reader.getFinerDomain();

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);

	Vector<2> vec(*domain, num_components);

	CHECK(vec.getNumLocalCells() == domain->getNumLocalCells());
}
TEST_CASE("Vector<2> getComponentView domain constructor", "[Vector]")
{
	Communicator comm(MPI_COMM_WORLD);

	auto                  mesh_file        = GENERATE(as<std::string>{}, MESHES);
	int                   num_components   = GENERATE(1, 2, 3);
	auto                  num_ghost_cells  = GENERATE(0, 1, 5);
	int                   nx               = GENERATE(1, 4, 5);
	int                   ny               = GENERATE(1, 4, 5);
	int                   component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
	int                   patch_stride     = component_stride * num_components;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost_cells);
	shared_ptr<Domain<2>> domain = domain_reader.getFinerDomain();

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);

	Vector<2> vec(*domain, num_components);

	double *view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
	for (int i = 0; i < domain->getNumLocalPatches(); i++) {
		INFO("i:                 " << i);
		for (int c = 0; c < num_components; c++) {
			INFO("c:                 " << c);
			ComponentView<double, 2> ld = vec.getComponentView(c, i);
			CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i + c * component_stride);
			CHECK(&ld[ld.getGhostEnd()]
			      == view + patch_stride * i + (c + 1) * component_stride - 1);
		}
	}
}
TEST_CASE("Vector<2> getComponentView const domain constructor", "[Vector]")
{
	Communicator comm(MPI_COMM_WORLD);

	auto                  mesh_file        = GENERATE(as<std::string>{}, MESHES);
	int                   num_components   = GENERATE(1, 2, 3);
	auto                  num_ghost_cells  = GENERATE(0, 1, 5);
	int                   nx               = GENERATE(1, 4, 5);
	int                   ny               = GENERATE(1, 4, 5);
	int                   component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
	int                   patch_stride     = component_stride * num_components;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost_cells);
	shared_ptr<Domain<2>> domain = domain_reader.getFinerDomain();

	INFO("num_ghost_cells:   " << num_ghost_cells);
	INFO("nx:                " << nx);
	INFO("ny:                " << ny);

	const Vector<2> vec(*domain, num_components);

	const double *view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
	for (int i = 0; i < domain->getNumLocalPatches(); i++) {
		INFO("i:                 " << i);
		for (int c = 0; c < num_components; c++) {
			INFO("c:                 " << c);
			ComponentView<const double, 2> ld = vec.getComponentView(c, i);
			CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i + c * component_stride);
			CHECK(&ld[ld.getGhostEnd()]
			      == view + patch_stride * i + (c + 1) * component_stride - 1);
		}
	}
}