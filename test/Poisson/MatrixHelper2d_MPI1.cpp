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
#include <ThunderEgg/BiQuadraticGhostFiller.h>
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/LinearRestrictor.h>
#include <ThunderEgg/PETSc/MatWrapper.h>
#include <ThunderEgg/Poisson/MatrixHelper2d.h>
#include <ThunderEgg/Poisson/StarPatchOperator.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

#define MESHES \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";

TEST_CASE("Poisson::MatrixHelper2d gives equivalent operator to Poisson::StarPatchOperator",
          "[Poisson::MatrixHelper2d]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int             n         = 32;
	int             num_ghost = 1;
	bitset<4>       neumann;
	DomainReader<2> domain_reader(mesh_file, {n, n}, num_ghost);
	Domain<2>       d_fine = domain_reader.getFinerDomain();

	auto gfun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return sinl(M_PI * y) * cosl(2 * M_PI * x);
	};

	Vector<2> f_vec(d_fine, 1);
	Vector<2> f_vec_expected(d_fine, 1);

	Vector<2> g_vec(d_fine, 1);
	DomainTools::SetValues<2>(d_fine, g_vec, gfun);

	BiQuadraticGhostFiller        gf(d_fine, GhostFillingType::Faces);
	Poisson::StarPatchOperator<2> p_operator(d_fine, gf);
	p_operator.apply(g_vec, f_vec_expected);

	// generate matrix with matrix_helper
	Poisson::MatrixHelper2d mh(d_fine, neumann);
	Mat                     A = mh.formCRSMatrix();
	PETSc::MatWrapper<2>    m_operator(A);
	m_operator.apply(g_vec, f_vec);

	REQUIRE(f_vec.infNorm() > 0);

	for (auto pinfo : d_fine.getPatchInfoVector()) {
		INFO("Patch: " << pinfo.id);
		INFO("x:     " << pinfo.starts[0]);
		INFO("y:     " << pinfo.starts[1]);
		INFO("nx:    " << pinfo.ns[0]);
		INFO("ny:    " << pinfo.ns[1]);
		INFO("dx:    " << pinfo.spacings[0]);
		INFO("dy:    " << pinfo.spacings[1]);
		ComponentView<double, 2> f_vec_ld          = f_vec.getComponentView(0, pinfo.local_index);
		ComponentView<double, 2> f_vec_expected_ld = f_vec_expected.getComponentView(0, pinfo.local_index);
		Loop::Nested<2>(f_vec_ld.getStart(), f_vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			INFO("xi:    " << coord[0]);
			INFO("yi:    " << coord[1]);
			CHECK(f_vec_ld[coord] == Catch::Approx(f_vec_expected_ld[coord]));
		});
	}
	MatDestroy(&A);
}
TEST_CASE(
"Poisson::MatrixHelper2d gives equivalent operator to Poisson::StarPatchOperator with Neumann BC",
"[Poisson::MatrixHelper2d]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int             n         = 32;
	int             num_ghost = 1;
	bitset<4>       neumann   = 0xF;
	DomainReader<2> domain_reader(mesh_file, {n, n}, num_ghost);
	Domain<2>       d_fine = domain_reader.getFinerDomain();

	auto gfun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return sinl(M_PI * y) * cosl(2 * M_PI * x);
	};

	Vector<2> f_vec(d_fine, 1);
	Vector<2> f_vec_expected(d_fine, 1);

	Vector<2> g_vec(d_fine, 1);
	DomainTools::SetValues<2>(d_fine, g_vec, gfun);

	BiQuadraticGhostFiller        gf(d_fine, GhostFillingType::Faces);
	Poisson::StarPatchOperator<2> p_operator(d_fine, gf, true);
	p_operator.apply(g_vec, f_vec_expected);

	// generate matrix with matrix_helper
	Poisson::MatrixHelper2d mh(d_fine, neumann);
	Mat                     A = mh.formCRSMatrix();
	PETSc::MatWrapper<2>    m_operator(A);
	m_operator.apply(g_vec, f_vec);

	REQUIRE(f_vec.infNorm() > 0);

	for (auto pinfo : d_fine.getPatchInfoVector()) {
		INFO("Patch: " << pinfo.id);
		INFO("x:     " << pinfo.starts[0]);
		INFO("y:     " << pinfo.starts[1]);
		INFO("nx:    " << pinfo.ns[0]);
		INFO("ny:    " << pinfo.ns[1]);
		INFO("dx:    " << pinfo.spacings[0]);
		INFO("dy:    " << pinfo.spacings[1]);
		ComponentView<double, 2> f_vec_ld          = f_vec.getComponentView(0, pinfo.local_index);
		ComponentView<double, 2> f_vec_expected_ld = f_vec_expected.getComponentView(0, pinfo.local_index);
		Loop::Nested<2>(f_vec_ld.getStart(), f_vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			INFO("xi:    " << coord[0]);
			INFO("yi:    " << coord[1]);
			CHECK(f_vec_ld[coord] == Catch::Approx(f_vec_expected_ld[coord]));
		});
	}
	MatDestroy(&A);
}