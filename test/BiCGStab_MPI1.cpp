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

#include "catch.hpp"
#include "utils/DomainReader.h"
#include <ThunderEgg/BiCGStab.h>
#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/Poisson/StarPatchOperator.h>
#include <ThunderEgg/ValVectorGenerator.h>
using namespace std;
using namespace ThunderEgg;

TEST_CASE("BiCGStab solves poisson problem withing given tolerance", "[BiCGStab]")
{
	string mesh_file = "mesh_inputs/2d_uniform_2x2_mpi1.json";
	INFO("MESH FILE " << mesh_file);
	DomainReader<2>       domain_reader(mesh_file, {32, 32}, 1);
	shared_ptr<Domain<2>> domain = domain_reader.getCoarserDomain();

	auto ffun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return -5 * M_PI * M_PI * sin(M_PI * y) * cos(2 * M_PI * x);
	};
	auto gfun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return sin(M_PI * y) * cos(2 * M_PI * x);
	};

	auto f_vec = ValVector<2>::GetNewVector(domain, 1);
	DomainTools::SetValues<2>(domain, f_vec, ffun);
	auto residual = ValVector<2>::GetNewVector(domain, 1);

	auto g_vec = ValVector<2>::GetNewVector(domain, 1);

	auto gf = make_shared<BiLinearGhostFiller>(domain);

	auto p_operator = make_shared<Poisson::StarPatchOperator<2>>(domain, gf);
	p_operator->addDrichletBCToRHS(f_vec, gfun);

	double tolerance = GENERATE(1e-9, 1e-7, 1e-5);

	BiCGStab<2>::solve(make_shared<ValVectorGenerator<2>>(domain, 1), p_operator, g_vec, f_vec,
	                   nullptr, 1000, tolerance);

	p_operator->apply(g_vec, residual);
	residual->addScaled(-1, f_vec);
	CHECK(residual->dot(residual) / f_vec->dot(f_vec) <= tolerance);
}
TEST_CASE("giving a good initial guess reduces the iterations", "[BiCGStab]")
{
	string mesh_file = "mesh_inputs/2d_uniform_2x2_mpi1.json";
	INFO("MESH FILE " << mesh_file);
	DomainReader<2>       domain_reader(mesh_file, {32, 32}, 1);
	shared_ptr<Domain<2>> domain = domain_reader.getCoarserDomain();

	auto ffun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return -5 * M_PI * M_PI * sin(M_PI * y) * cos(2 * M_PI * x);
	};
	auto gfun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return sin(M_PI * y) * cos(2 * M_PI * x);
	};

	auto f_vec = ValVector<2>::GetNewVector(domain, 1);
	DomainTools::SetValues<2>(domain, f_vec, ffun);
	auto residual = ValVector<2>::GetNewVector(domain, 1);

	auto g_vec = ValVector<2>::GetNewVector(domain, 1);

	auto gf = make_shared<BiLinearGhostFiller>(domain);

	auto p_operator = make_shared<Poisson::StarPatchOperator<2>>(domain, gf);
	p_operator->addDrichletBCToRHS(f_vec, gfun);

	double tolerance = 1e-5;

	BiCGStab<2>::solve(make_shared<ValVectorGenerator<2>>(domain, 1), p_operator, g_vec, f_vec,
	                   nullptr, 1000, tolerance);

	int iterations_with_solved_guess
	= BiCGStab<2>::solve(make_shared<ValVectorGenerator<2>>(domain, 1), p_operator, g_vec, f_vec,
	                     nullptr, 1000, tolerance);

	CHECK(iterations_with_solved_guess == 0);
}