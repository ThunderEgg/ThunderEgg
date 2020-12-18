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
#include <sstream>
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
TEST_CASE("outputs iteration count and residual to output", "[BiCGStab]")
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

	double tolerance = 1e-7;

	std::stringstream ss;

	BiCGStab<2>::solve(make_shared<ValVectorGenerator<2>>(domain, 1), p_operator, g_vec, f_vec,
	                   nullptr, 1000, tolerance, nullptr, true, ss);

	INFO(ss.str());
	int    prev_iteration;
	double resid;
	ss >> prev_iteration >> resid;
	while (prev_iteration < 18) {
		int iteration;
		ss >> iteration >> resid;
		CHECK(iteration == prev_iteration + 1);
		prev_iteration = iteration;
	}
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
namespace
{
class MockVector : public Vector<2>
{
	public:
	mutable int norm_calls = 0;
	double      dot_value;
	MockVector(double dot_value) : Vector<2>(MPI_COMM_WORLD, 1, 0, 10), dot_value(dot_value) {}
	LocalData<2>       getLocalData(int, int) override {}
	const LocalData<2> getLocalData(int, int) const override {}
	double             dot(std::shared_ptr<const Vector<2>>) const override
	{
		return dot_value;
	}
	double twoNorm() const override
	{
		norm_calls++;
		if (norm_calls == 1) {
			return 1;
		} else {
			return norm_calls * 1e6;
		}
	}
};
class MockVectorGenerator : public VectorGenerator<2>
{
	public:
	std::shared_ptr<MockVector> vec;
	MockVectorGenerator(std::shared_ptr<MockVector> vec) : vec(vec) {}
	std::shared_ptr<Vector<2>> getNewVector() const override
	{
		return vec;
	}
};
class MockOperator : public Operator<2>
{
	public:
	void apply(std::shared_ptr<const Vector<2>>, std::shared_ptr<Vector<2>>) const override {}
};
} // namespace
TEST_CASE("throws breakdown exception when rho is 0", "[BiCGStab]")
{
	auto vec = make_shared<MockVector>(0);
	CHECK_THROWS_AS(BiCGStab<2>::solve(make_shared<MockVectorGenerator>(vec),
	                                   make_shared<MockOperator>(), vec, vec),
	                BreakdownError);
}
TEST_CASE("throws divergence exception when residual keeps increasing", "[BiCGStab]")
{
	auto vec = make_shared<MockVector>(1);
	CHECK_THROWS_AS(BiCGStab<2>::solve(make_shared<MockVectorGenerator>(vec),
	                                   make_shared<MockOperator>(), vec, vec),
	                DivergenceError);
}