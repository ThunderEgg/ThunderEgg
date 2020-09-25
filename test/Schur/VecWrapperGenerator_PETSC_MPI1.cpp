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
#include <ThunderEgg/Schur/VecWrapperGenerator.h>
#include <limits>
using namespace std;
using namespace ThunderEgg;
#define MESHES "mesh_inputs/2d_refined_east_1x2_mpi1.json", "mesh_inputs/2d_uniform_1x2_mpi1.json"
TEST_CASE("Schur::VecWrapperGenerator<2> works for various meshes", "[Schur::VecWrapperGenerator]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	DomainReader<2> domain_reader(mesh_file, {10, 10}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	Schur::VecWrapperGenerator<1> vg(iface_domain);

	auto vector = vg.getNewVector();

	CHECK(vector->getNumLocalPatches() == iface_domain->getNumLocalInterfaces());
	CHECK(vector->getNumLocalCells() == 10 * iface_domain->getNumLocalInterfaces());
}
TEST_CASE("Schur::VecWrapperGenerator<2> throws exception for non-square patches",
          "[Schur::VecWrapperGenerator]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto            nx = GENERATE(5, 7);
	auto            ny = GENERATE(6, 8);
	DomainReader<2> domain_reader(mesh_file, {nx, ny}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	CHECK_THROWS_AS(Schur::VecWrapperGenerator<1>(iface_domain), RuntimeError);
}