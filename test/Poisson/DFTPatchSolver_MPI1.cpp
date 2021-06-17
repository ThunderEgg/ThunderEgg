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
#include <ThunderEgg/GMG/LinearRestrictor.h>
#include <ThunderEgg/Poisson/DFTPatchSolver.h>
#include <ThunderEgg/Poisson/StarPatchOperator.h>
#include <ThunderEgg/ValVector.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

#define MESHES \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";

TEST_CASE("Test Poisson::DFTPatchSolver gets 2nd order convergence", "[Poisson::StarPatchOperator]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	auto nx = GENERATE(10, 13);
	auto ny = GENERATE(10, 13);
	INFO("NX        " << nx);
	INFO("NY        " << ny);
	int       num_ghost = 1;
	bitset<4> neumann;
	double    errors[2];
	for (int i = 1; i <= 2; i++) {
		INFO("MULT      " << i);
		DomainReader<2>       domain_reader(mesh_file, {i * nx, i * ny}, num_ghost);
		shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

		auto ffun = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return -5 * M_PI * M_PI * sinl(M_PI * y) * cosl(2 * M_PI * x);
		};
		auto gfun = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return sinl(M_PI * y) * cosl(2 * M_PI * x);
		};

		auto g_vec = ValVector<2>::GetNewVector(d_fine, 1);
		DomainTools::SetValuesWithGhost<2>(d_fine, g_vec, gfun);
		auto g_vec_expected = ValVector<2>::GetNewVector(d_fine, 1);
		DomainTools::SetValues<2>(d_fine, g_vec_expected, gfun);

		auto f_vec = ValVector<2>::GetNewVector(d_fine, 1);
		DomainTools::SetValues<2>(d_fine, f_vec, ffun);

		auto gf         = make_shared<BiLinearGhostFiller>(d_fine, GhostFillingType::Faces);
		auto p_operator = make_shared<Poisson::StarPatchOperator<2>>(d_fine, gf);
		auto p_solver   = make_shared<Poisson::DFTPatchSolver<2>>(p_operator, neumann);
		p_operator->addDrichletBCToRHS(f_vec, gfun);

		p_solver->smooth(*f_vec, *g_vec);

		auto error_vec = ValVector<2>::GetNewVector(d_fine, 1);
		error_vec->addScaled(1.0, *g_vec, -1.0, *g_vec_expected);
		errors[i - 1] = error_vec->twoNorm() / g_vec_expected->twoNorm();
	}
	INFO("Errors: " << errors[0] << ", " << errors[1]);
	CHECK(log(errors[0] / errors[1]) / log(2) > 1.8);
}
TEST_CASE("Test Poisson::DFTPatchSolver gets 2nd order convergence with neumann boundary",
          "[Poisson::StarPatchOperator]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	auto nx = GENERATE(10, 13);
	auto ny = GENERATE(10, 13);
	INFO("NX        " << nx);
	INFO("NY        " << ny);
	int       num_ghost = 1;
	bitset<4> neumann   = 0xF;
	double    errors[2];
	for (int i = 1; i <= 2; i++) {
		INFO("MULT      " << i);
		DomainReader<2>       domain_reader(mesh_file, {i * nx, i * ny}, num_ghost);
		shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

		auto ffun = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return -5 * M_PI * M_PI * sinl(M_PI * y) * cosl(2 * M_PI * x);
		};
		auto gfun = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return sinl(M_PI * y) * cosl(2 * M_PI * x);
		};
		auto gfun_x = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return -2 * M_PI * sin(M_PI * y) * sin(2 * M_PI * x);
		};
		auto gfun_y = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return M_PI * cos(M_PI * y) * cos(2 * M_PI * x);
		};

		auto g_vec = ValVector<2>::GetNewVector(d_fine, 1);
		DomainTools::SetValuesWithGhost<2>(d_fine, g_vec, gfun);
		auto g_vec_expected = ValVector<2>::GetNewVector(d_fine, 1);
		DomainTools::SetValues<2>(d_fine, g_vec_expected, gfun);

		auto f_vec = ValVector<2>::GetNewVector(d_fine, 1);
		DomainTools::SetValues<2>(d_fine, f_vec, ffun);

		auto gf         = make_shared<BiLinearGhostFiller>(d_fine, GhostFillingType::Faces);
		auto p_operator = make_shared<Poisson::StarPatchOperator<2>>(d_fine, gf, true);
		auto p_solver   = make_shared<Poisson::DFTPatchSolver<2>>(p_operator, neumann);
		p_operator->addNeumannBCToRHS(f_vec, gfun, {gfun_x, gfun_y});

		p_solver->smooth(*f_vec, *g_vec);

		auto error_vec = ValVector<2>::GetNewVector(d_fine, 1);
		g_vec->shift(-d_fine->integrate(g_vec) / d_fine->volume());
		g_vec_expected->shift(-d_fine->integrate(g_vec_expected) / d_fine->volume());
		error_vec->addScaled(1.0, *g_vec, -1.0, *g_vec_expected);
		errors[i - 1] = error_vec->twoNorm() / g_vec_expected->twoNorm();
	}
	INFO("Errors: " << errors[0] << ", " << errors[1]);
	CHECK(log(errors[0] / errors[1]) / log(2) > 1.8);
}
TEST_CASE(
"Test Poisson::DFTPatchSolver gets 2nd order convergence with neumann boundary single patch",
"[Poisson::StarPatchOperator]")
{
	auto mesh_file = "mesh_inputs/2d_uniform_2x2_mpi1.json";
	INFO("MESH FILE " << mesh_file);
	auto nx = GENERATE(10, 13);
	auto ny = GENERATE(10, 13);
	INFO("NX        " << nx);
	INFO("NY        " << ny);
	int       num_ghost = 1;
	bitset<4> neumann   = 0xF;
	double    errors[2];
	for (int i = 1; i <= 2; i++) {
		INFO("MULT      " << i);
		DomainReader<2>       domain_reader(mesh_file, {i * nx, i * ny}, num_ghost);
		shared_ptr<Domain<2>> d_fine = domain_reader.getCoarserDomain();

		auto ffun = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return -5 * M_PI * M_PI * sinl(M_PI * y) * cosl(2 * M_PI * x);
		};
		auto gfun = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return sinl(M_PI * y) * cosl(2 * M_PI * x);
		};
		auto gfun_x = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return -2 * M_PI * sin(M_PI * y) * sin(2 * M_PI * x);
		};
		auto gfun_y = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return M_PI * cos(M_PI * y) * cos(2 * M_PI * x);
		};

		auto g_vec = ValVector<2>::GetNewVector(d_fine, 1);
		DomainTools::SetValuesWithGhost<2>(d_fine, g_vec, gfun);
		auto g_vec_expected = ValVector<2>::GetNewVector(d_fine, 1);
		DomainTools::SetValues<2>(d_fine, g_vec_expected, gfun);

		auto f_vec = ValVector<2>::GetNewVector(d_fine, 1);
		DomainTools::SetValues<2>(d_fine, f_vec, ffun);

		auto gf         = make_shared<BiLinearGhostFiller>(d_fine, GhostFillingType::Faces);
		auto p_operator = make_shared<Poisson::StarPatchOperator<2>>(d_fine, gf, true);
		auto p_solver   = make_shared<Poisson::DFTPatchSolver<2>>(p_operator, neumann);
		p_operator->addNeumannBCToRHS(f_vec, gfun, {gfun_x, gfun_y});

		p_solver->smooth(*f_vec, *g_vec);

		auto error_vec = ValVector<2>::GetNewVector(d_fine, 1);
		g_vec->shift(-d_fine->integrate(g_vec) / d_fine->volume());
		g_vec_expected->shift(-d_fine->integrate(g_vec_expected) / d_fine->volume());
		error_vec->addScaled(1.0, *g_vec, -1.0, *g_vec_expected);
		errors[i - 1] = error_vec->twoNorm() / g_vec_expected->twoNorm();
	}
	INFO("Errors: " << errors[0] << ", " << errors[1]);
	CHECK(log(errors[0] / errors[1]) / log(2) > 1.8);
}