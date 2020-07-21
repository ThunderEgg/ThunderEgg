/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
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
#include <Thunderegg/BiLinearGhostFiller.h>
#include <Thunderegg/DomainTools.h>
#include <Thunderegg/GMG/LinearRestrictor.h>
#include <Thunderegg/Poisson/DftPatchSolver.h>
#include <Thunderegg/Poisson/StarPatchOperator.h>
#include <Thunderegg/ValVector.h>
using namespace std;
using namespace Thunderegg;
#define MESHES                                                                                     \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";
TEST_CASE("Test Poisson::DftPatchSolver gets 2nd order convergence", "[Poisson::StarPatchOperator]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	auto nx = GENERATE(10, 13);
	auto ny = GENERATE(10, 13);
	INFO("NX        " << nx);
	INFO("NY        " << ny);
	int    num_ghost = 1;
	double errors[2];
	for (int i = 1; i <= 2; i++) {
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

		auto g_vec = ValVector<2>::GetNewVector(d_fine);
		DomainTools<2>::setValuesWithGhost(d_fine, g_vec, gfun);
		auto g_vec_expected = ValVector<2>::GetNewVector(d_fine);
		DomainTools<2>::setValues(d_fine, g_vec_expected, gfun);

		auto f_vec = ValVector<2>::GetNewVector(d_fine);
		DomainTools<2>::setValues(d_fine, f_vec, ffun);
		Poisson::StarPatchOperator<2>::addDrichletBCToRHS(d_fine, f_vec, gfun);

		auto gf         = make_shared<BiLinearGhostFiller>(d_fine);
		auto p_operator = make_shared<Poisson::StarPatchOperator<2>>(d_fine, gf);
		auto p_solver   = make_shared<Poisson::DftPatchSolver<2>>(p_operator);
		p_solver->smooth(f_vec, g_vec);

		auto error_vec = ValVector<2>::GetNewVector(d_fine);
		error_vec->addScaled(1.0, g_vec, -1.0, g_vec_expected);
		errors[i] = error_vec->twoNorm() / g_vec_expected->twoNorm();
	}
	INFO("Errors: " << errors[0] << ", " << errors[1]);
	CHECK(log(errors[0] / errors[1]) / log(2) > 1.8);
}