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
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/LinearRestrictor.h>
#include <ThunderEgg/PETSc/MatWrapper.h>
#include <ThunderEgg/PETSc/VecWrapper.h>
#include <ThunderEgg/Poisson/MatrixHelper.h>
#include <ThunderEgg/Poisson/StarPatchOperator.h>
#include <ThunderEgg/TriLinearGhostFiller.h>
using namespace std;
using namespace ThunderEgg;
#define MESHES                                                                                     \
	"mesh_inputs/3d_uniform_2x2x2_mpi1.json", "mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json",        \
	"mesh_inputs/3d_mid_refine_4x4x4_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";
TEST_CASE("Poisson::MatrixHelper gives equivalent operator to Poisson::StarPatchOperator",
          "[Poisson::MatrixHelper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	auto                  nx        = GENERATE(8, 10);
	auto                  ny        = GENERATE(8, 10);
	auto                  nz        = GENERATE(8, 10);
	int                   num_ghost = 1;
	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d_fine = domain_reader.getFinerDomain();

	auto gfun = [](const std::array<double, 3> &coord) {
		double x = coord[0];
		double y = coord[1];
		double z = coord[2];
		return sin(M_PI * y) * cos(2 * M_PI * x) * cos(M_PI * z);
	};

	auto f_vec          = PETSc::VecWrapper<3>::GetNewVector(d_fine);
	auto f_vec_expected = PETSc::VecWrapper<3>::GetNewVector(d_fine);

	auto g_vec = PETSc::VecWrapper<3>::GetNewVector(d_fine);
	DomainTools<3>::setValues(d_fine, g_vec, gfun);

	auto gf         = make_shared<TriLinearGhostFiller>(d_fine);
	auto p_operator = make_shared<Poisson::StarPatchOperator<3>>(d_fine, gf);
	p_operator->apply(g_vec, f_vec_expected);

	// generate matrix with matrix_helper
	Poisson::MatrixHelper mh(d_fine);
	Mat                   A          = mh.formCRSMatrix();
	auto                  m_operator = make_shared<PETSc::MatWrapper<3>>(A);
	m_operator->apply(g_vec, f_vec);

	REQUIRE(f_vec->infNorm() > 0);

	for (auto pinfo : d_fine->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		INFO("nz:    " << pinfo->ns[1]);
		INFO("dx:    " << pinfo->spacings[0]);
		INFO("dy:    " << pinfo->spacings[1]);
		INFO("dz:    " << pinfo->spacings[1]);
		LocalData<3> f_vec_ld          = f_vec->getLocalData(0, pinfo->local_index);
		LocalData<3> f_vec_expected_ld = f_vec_expected->getLocalData(0, pinfo->local_index);
		nested_loop<3>(f_vec_ld.getStart(), f_vec_ld.getEnd(), [&](const array<int, 3> &coord) {
			INFO("xi:    " << coord[0]);
			INFO("yi:    " << coord[1]);
			INFO("zi:    " << coord[2]);
			CHECK(f_vec_ld[coord] == Approx(f_vec_expected_ld[coord]));
		});
	}
	MatDestroy(&A);
}
TEST_CASE(
"Poisson::MatrixHelper gives equivalent operator to Poisson::StarPatchOperator with Neumann BC",
"[Poisson::MatrixHelper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	auto                  nx        = GENERATE(8, 10);
	auto                  ny        = GENERATE(8, 10);
	auto                  nz        = GENERATE(8, 10);
	int                   num_ghost = 1;
	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost, true);
	shared_ptr<Domain<3>> d_fine = domain_reader.getFinerDomain();

	auto gfun = [](const std::array<double, 3> &coord) {
		double x = coord[0];
		double y = coord[1];
		double z = coord[2];
		return sin(M_PI * y) * cos(2 * M_PI * x) * cos(M_PI * z);
	};

	auto f_vec          = PETSc::VecWrapper<3>::GetNewVector(d_fine);
	auto f_vec_expected = PETSc::VecWrapper<3>::GetNewVector(d_fine);

	auto g_vec = PETSc::VecWrapper<3>::GetNewVector(d_fine);
	DomainTools<3>::setValues(d_fine, g_vec, gfun);

	auto gf         = make_shared<TriLinearGhostFiller>(d_fine);
	auto p_operator = make_shared<Poisson::StarPatchOperator<3>>(d_fine, gf, true);
	p_operator->apply(g_vec, f_vec_expected);

	// generate matrix with matrix_helper
	Poisson::MatrixHelper mh(d_fine);
	Mat                   A          = mh.formCRSMatrix();
	auto                  m_operator = make_shared<PETSc::MatWrapper<3>>(A);
	m_operator->apply(g_vec, f_vec);

	REQUIRE(f_vec->infNorm() > 0);

	for (auto pinfo : d_fine->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("z:     " << pinfo->starts[2]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		INFO("nz:    " << pinfo->ns[2]);
		INFO("dx:    " << pinfo->spacings[0]);
		INFO("dy:    " << pinfo->spacings[1]);
		INFO("dz:    " << pinfo->spacings[2]);
		LocalData<3> f_vec_ld          = f_vec->getLocalData(0, pinfo->local_index);
		LocalData<3> f_vec_expected_ld = f_vec_expected->getLocalData(0, pinfo->local_index);
		nested_loop<3>(f_vec_ld.getStart(), f_vec_ld.getEnd(), [&](const array<int, 3> &coord) {
			INFO("xi:    " << coord[0]);
			INFO("yi:    " << coord[1]);
			INFO("zi:    " << coord[2]);
			CHECK(f_vec_ld[coord] == Approx(f_vec_expected_ld[coord]));
		});
	}
	MatDestroy(&A);
}
TEST_CASE("Poisson::MatrixHelper constructor throws error with odd number of cells",
          "[Poisson::MatrixHelper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto axis = GENERATE(0, 1, 2);
	INFO("axis: " << axis);
	int n_even    = 10;
	int n_odd     = 11;
	int num_ghost = 1;

	array<int, 3> ns;
	ns.fill(n_even);
	ns[axis] = n_odd;
	DomainReader<3>       domain_reader(mesh_file, ns, num_ghost);
	shared_ptr<Domain<3>> d = domain_reader.getFinerDomain();

	CHECK_THROWS_AS(Poisson::MatrixHelper(d), RuntimeError);
}