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
#include <ThunderEgg/Iterative/BiCGStab.h>
#include <ThunderEgg/Poisson/StarPatchOperator.h>
#include <ThunderEgg/ValVectorGenerator.h>

#include <sstream>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::Iterative;

TEST_CASE("BiCGStab default max iterations", "[BiCGStab]")
{
	BiCGStab<2> bcgs;
	CHECK(bcgs.getMaxIterations() == 1000);
}
TEST_CASE("BiCGStab set max iterations", "[BiCGStab]")
{
	BiCGStab<2> bcgs;
	int         iterations = GENERATE(1, 2, 3);
	bcgs.setMaxIterations(iterations);
	CHECK(bcgs.getMaxIterations() == iterations);
}
TEST_CASE("BiCGStab default tolerance", "[BiCGStab]")
{
	BiCGStab<2> bcgs;
	CHECK(bcgs.getTolerance() == 1e-12);
}
TEST_CASE("BiCGStab set tolerance", "[BiCGStab]")
{
	BiCGStab<2> bcgs;
	double      tolerance = GENERATE(1.2, 2.3, 3.4);
	bcgs.setTolerance(tolerance);
	CHECK(bcgs.getTolerance() == tolerance);
}
TEST_CASE("BiCGStab default timer", "[BiCGStab]")
{
	BiCGStab<2> bcgs;
	CHECK(bcgs.getTimer() == nullptr);
}
TEST_CASE("BiCGStab set timer", "[BiCGStab]")
{
	BiCGStab<2> bcgs;
	auto        timer = make_shared<Timer>(MPI_COMM_WORLD);
	bcgs.setTimer(timer);
	CHECK(bcgs.getTimer() == timer);
}
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

	auto gf = make_shared<BiLinearGhostFiller>(domain, GhostFillingType::Faces);

	auto p_operator = make_shared<Poisson::StarPatchOperator<2>>(domain, gf);
	p_operator->addDrichletBCToRHS(f_vec, gfun);

	double tolerance = GENERATE(1e-9, 1e-7, 1e-5);

	BiCGStab<2> solver;
	solver.setMaxIterations(1000);
	solver.setTolerance(tolerance);
	solver.solve(make_shared<ValVectorGenerator<2>>(domain, 1), p_operator, g_vec, f_vec);

	p_operator->apply(g_vec, residual);
	residual->addScaled(-1, f_vec);
	CHECK(residual->dot(residual) / f_vec->dot(f_vec) <= tolerance);
}
TEST_CASE("BiCGStab handles zero rhs vector", "[BiCGStab]")
{
	string mesh_file = "mesh_inputs/2d_uniform_2x2_mpi1.json";
	INFO("MESH FILE " << mesh_file);
	DomainReader<2>       domain_reader(mesh_file, {32, 32}, 1);
	shared_ptr<Domain<2>> domain = domain_reader.getCoarserDomain();

	auto f_vec = ValVector<2>::GetNewVector(domain, 1);

	auto g_vec = ValVector<2>::GetNewVector(domain, 1);

	auto gf = make_shared<BiLinearGhostFiller>(domain, GhostFillingType::Faces);

	auto p_operator = make_shared<Poisson::StarPatchOperator<2>>(domain, gf);

	double tolerance = GENERATE(1e-9, 1e-7, 1e-5);

	BiCGStab<2> solver;
	solver.setMaxIterations(1000);
	solver.setTolerance(tolerance);
	solver.solve(make_shared<ValVectorGenerator<2>>(domain, 1), p_operator, g_vec, f_vec);

	CHECK(g_vec->infNorm() == 0);
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

	auto gf = make_shared<BiLinearGhostFiller>(domain, GhostFillingType::Faces);

	auto p_operator = make_shared<Poisson::StarPatchOperator<2>>(domain, gf);
	p_operator->addDrichletBCToRHS(f_vec, gfun);

	double tolerance = 1e-7;

	std::stringstream ss;

	BiCGStab<2> solver;
	solver.setMaxIterations(1000);
	solver.setTolerance(tolerance);
	solver.solve(make_shared<ValVectorGenerator<2>>(domain, 1), p_operator, g_vec, f_vec, nullptr,
	             true, ss);

	INFO(ss.str());
	int    prev_iteration;
	double resid;
	ss >> prev_iteration >> resid;
	while (prev_iteration < 5) {
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

	auto gf = make_shared<BiLinearGhostFiller>(domain, GhostFillingType::Faces);

	auto p_operator = make_shared<Poisson::StarPatchOperator<2>>(domain, gf);
	p_operator->addDrichletBCToRHS(f_vec, gfun);

	double tolerance = 1e-5;

	BiCGStab<2> solver;
	solver.setMaxIterations(1000);
	solver.setTolerance(tolerance);
	solver.solve(make_shared<ValVectorGenerator<2>>(domain, 1), p_operator, g_vec, f_vec);

	int iterations_with_solved_guess
	= solver.solve(make_shared<ValVectorGenerator<2>>(domain, 1), p_operator, g_vec, f_vec);

	CHECK(iterations_with_solved_guess == 0);
}
namespace
{
class MockVector : public Vector<2>
{
	public:
	mutable int norm_calls = 0;
	double      dot_value;
	MockVector(double dot_value)
	: Vector<2>(MPI_COMM_WORLD, 1, 0, 10), dot_value(dot_value) {}
	LocalData<2> getLocalData(int, int) override
	{
		return LocalData<2>();
	}
	const LocalData<2> getLocalData(int, int) const override
	{
		return LocalData<2>();
	}
	double dot(std::shared_ptr<const Vector<2>>) const override
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
	MockVectorGenerator(std::shared_ptr<MockVector> vec)
	: vec(vec) {}
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
class I2Operator : public Operator<2>
{
	public:
	void apply(std::shared_ptr<const Vector<2>> x, std::shared_ptr<Vector<2>> y) const override
	{
		y->copy(x);
		y->scale(2);
	}
};
} // namespace
TEST_CASE("BiCGStab solves poisson 2I problem", "[BiCGStab]")
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

	auto f_vec = ValVector<2>::GetNewVector(domain, 1);
	DomainTools::SetValues<2>(domain, f_vec, ffun);
	auto residual = ValVector<2>::GetNewVector(domain, 1);

	auto g_vec = ValVector<2>::GetNewVector(domain, 1);

	auto gf = make_shared<BiLinearGhostFiller>(domain, GhostFillingType::Faces);

	auto op = make_shared<I2Operator>();

	double tolerance = GENERATE(1e-9, 1e-7, 1e-5);

	BiCGStab<2> solver;
	solver.setMaxIterations(1000);
	solver.setTolerance(tolerance);
	solver.solve(make_shared<ValVectorGenerator<2>>(domain, 1), op, g_vec, f_vec);

	op->apply(g_vec, residual);
	residual->addScaled(-1, f_vec);
	CHECK(residual->dot(residual) / f_vec->dot(f_vec) <= tolerance);
}
TEST_CASE("throws breakdown exception when rho is 0", "[BiCGStab]")
{
	auto        vec = make_shared<MockVector>(0);
	BiCGStab<2> solver;
	CHECK_THROWS_AS(
	solver.solve(make_shared<MockVectorGenerator>(vec), make_shared<MockOperator>(), vec, vec),
	BreakdownError);
}