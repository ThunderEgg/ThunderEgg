/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
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
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/PETSc/MatWrapper.h>

#include <petscmat.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

#define MESHES \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";

TEST_CASE("PETSc::MatWrapper works with ValVector and 0.5I", "[PETSc::MatWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int             n         = 32;
	int             num_ghost = 1;
	DomainReader<2> domain_reader(mesh_file, {n, n}, num_ghost);
	Domain<2>       d_fine = domain_reader.getFinerDomain();

	auto gfun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return sinl(M_PI * y) * cosl(2 * M_PI * x);
	};

	Vector<2> x(d_fine, 1);
	DomainTools::SetValues<2>(d_fine, x, gfun);
	Vector<2> b(d_fine, 1);

	// create an Identity matrix
	Mat A;
	MatCreateAIJ(MPI_COMM_WORLD, d_fine.getNumLocalCells(), d_fine.getNumLocalCells(),
	             PETSC_DETERMINE, PETSC_DETERMINE, 1, nullptr, 1, nullptr, &A);
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	MatShift(A, 0.5);

	// create MatWrapper
	PETSc::MatWrapper<2> m_operator(A);
	m_operator.apply(x, b);

	for (auto pinfo : d_fine.getPatchInfoVector()) {
		INFO("Patch: " << pinfo.id);
		INFO("x:     " << pinfo.starts[0]);
		INFO("y:     " << pinfo.starts[1]);
		INFO("nx:    " << pinfo.ns[0]);
		INFO("ny:    " << pinfo.ns[1]);
		INFO("dx:    " << pinfo.spacings[0]);
		INFO("dy:    " << pinfo.spacings[1]);
		ComponentView<double, 2> x_ld = x.getComponentView(0, pinfo.local_index);
		ComponentView<double, 2> b_ld = b.getComponentView(0, pinfo.local_index);
		Loop::Nested<2>(x_ld.getStart(), x_ld.getEnd(), [&](const array<int, 2> &coord) {
			INFO("xi:    " << coord[0]);
			INFO("yi:    " << coord[1]);
			CHECK(0.5 * x_ld[coord] == Catch::Approx(b_ld[coord]));
		});
	}
	MatDestroy(&A);
}
TEST_CASE("PETSc::MatWrapper works with ValVector and 0.5I two components", "[PETSc::MatWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int             n         = 32;
	int             num_ghost = 1;
	DomainReader<2> domain_reader(mesh_file, {n, n}, num_ghost);
	Domain<2>       d_fine = domain_reader.getFinerDomain();

	auto gfun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return sinl(M_PI * y) * cosl(2 * M_PI * x);
	};
	auto ffun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return x + y;
	};

	Vector<2> x(d_fine, 2);
	DomainTools::SetValues<2>(d_fine, x, gfun, ffun);
	Vector<2> b(d_fine, 2);

	// create an Identity matrix
	Mat A;
	MatCreateAIJ(MPI_COMM_WORLD, d_fine.getNumLocalCells() * 2, d_fine.getNumLocalCells() * 2,
	             PETSC_DETERMINE, PETSC_DETERMINE, 1, nullptr, 1, nullptr, &A);
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	MatShift(A, 0.5);

	// create MatWrapper
	PETSc::MatWrapper<2> m_operator(A);
	m_operator.apply(x, b);

	for (auto pinfo : d_fine.getPatchInfoVector()) {
		INFO("Patch: " << pinfo.id);
		INFO("x:     " << pinfo.starts[0]);
		INFO("y:     " << pinfo.starts[1]);
		INFO("nx:    " << pinfo.ns[0]);
		INFO("ny:    " << pinfo.ns[1]);
		INFO("dx:    " << pinfo.spacings[0]);
		INFO("dy:    " << pinfo.spacings[1]);
		ComponentView<double, 2> x_ld = x.getComponentView(0, pinfo.local_index);
		ComponentView<double, 2> b_ld = b.getComponentView(0, pinfo.local_index);
		Loop::Nested<2>(x_ld.getStart(), x_ld.getEnd(), [&](const array<int, 2> &coord) {
			INFO("xi:    " << coord[0]);
			INFO("yi:    " << coord[1]);
			CHECK(0.5 * x_ld[coord] == Catch::Approx(b_ld[coord]));
		});
		ComponentView<double, 2> x_ld2 = x.getComponentView(1, pinfo.local_index);
		ComponentView<double, 2> b_ld2 = b.getComponentView(1, pinfo.local_index);
		Loop::Nested<2>(x_ld.getStart(), x_ld.getEnd(), [&](const array<int, 2> &coord) {
			INFO("xi:    " << coord[0]);
			INFO("yi:    " << coord[1]);
			CHECK(0.5 * x_ld2[coord] == Catch::Approx(b_ld2[coord]));
		});
	}
	MatDestroy(&A);
}