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
#include <Thunderegg/Poisson/StarPatchOperator.h>
#include <Thunderegg/ValVector.h>
using namespace std;
using namespace Thunderegg;
#define MESHES                                                                                     \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";
TEST_CASE("Test Poisson::StarPatchOperator add ghost to RHS", "[Poisson::StarPatchOperator]")
{
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
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
	auto hfun = [](const std::array<double, 2> &coord) { return 1; };

	auto f_vec = ValVector<2>::GetNewVector(d_fine);
	DomainTools<2>::setValuesWithGhost(d_fine, f_vec, ffun);

	auto g_vec = ValVector<2>::GetNewVector(d_fine);
	DomainTools<2>::setValuesWithGhost(d_fine, g_vec, gfun);

	auto g_zero = ValVector<2>::GetNewVector(d_fine);
	DomainTools<2>::setValuesWithGhost(d_fine, g_zero, gfun);
	g_zero->set(0);

	auto h_vec = ValVector<2>::GetNewVector(d_fine);
	DomainTools<2>::setValuesWithGhost(d_fine, h_vec, hfun);

	shared_ptr<BiLinearGhostFiller>           gf(new BiLinearGhostFiller(d_fine));
	shared_ptr<Poisson::StarPatchOperator<2>> p_operator(
	new Poisson::StarPatchOperator<2>(d_fine, gf));
	Poisson::StarPatchOperator<2>::addDrichletBCToRHS(d_fine, f_vec, gfun);

	auto f_expected = ValVector<2>::GetNewVector(d_fine);
	f_expected->copy(f_vec);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		auto u = g_vec->getLocalData(pinfo->local_index);
		auto f = f_expected->getLocalData(pinfo->local_index);
		for (Side<2> s : Side<2>::getValues()) {
			if (pinfo->hasNbr(s)) {
				double h2      = std::pow(pinfo->spacings[s.getAxisIndex()], 2);
				auto   f_slice = f.getSliceOnSide(s);
				auto   u_inner = u.getSliceOnSide(s);
				auto   u_ghost = u.getSliceOnSide(s, -1);
				nested_loop<1>(f_slice.getStart(), f_slice.getEnd(), [&](std::array<int, 1> coord) {
					f_slice[coord] += -(u_inner[coord] + u_ghost[coord]) / (h2);
				});
			}
		}
	}

	for (auto pinfo : d_fine->getPatchInfoVector()) {
		p_operator->addGhostToRHS(pinfo, g_vec->getLocalData(pinfo->local_index),
		                          f_vec->getLocalData(pinfo->local_index));
	}

	for (auto pinfo : d_fine->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		LocalData<2> vec_ld      = f_vec->getLocalData(pinfo->local_index);
		LocalData<2> expected_ld = f_expected->getLocalData(pinfo->local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			INFO("xi:    " << coord[0]);
			INFO("yi:    " << coord[1]);
			REQUIRE(vec_ld[coord] == Approx(expected_ld[coord]));
		});
	}
}
TEST_CASE("Test Poisson::StarPatchOperator apply on linear lhs constant coeff",
          "[Poisson::StarPatchOperator]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto ffun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return 1 + 5 * x + 13 * y;
	};
	auto gfun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return 0;
	};
	auto hfun = [](const std::array<double, 2> &coord) { return 1; };

	auto f_vec = ValVector<2>::GetNewVector(d_fine);

	auto g_vec = ValVector<2>::GetNewVector(d_fine);
	DomainTools<2>::setValuesWithGhost(d_fine, g_vec, gfun);

	auto h_vec = ValVector<2>::GetNewVector(d_fine);
	DomainTools<2>::setValuesWithGhost(d_fine, h_vec, hfun);

	shared_ptr<BiLinearGhostFiller>           gf(new BiLinearGhostFiller(d_fine));
	shared_ptr<Poisson::StarPatchOperator<2>> p_operator(
	new Poisson::StarPatchOperator<2>(d_fine, gf));
	p_operator->apply(f_vec, g_vec);

	for (auto pinfo : d_fine->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		LocalData<2> vec_ld = g_vec->getLocalData(pinfo->local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			INFO("xi:    " << coord[0]);
			INFO("yi:    " << coord[1]);
			CHECK(vec_ld[coord] + 1 == Approx(1));
		});
	}
}
TEST_CASE("Test Poisson::StarPatchOperator gets 2nd order convergence",
          "[Poisson::StarPatchOperator]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int    ns[2]     = {32, 64};
	int    num_ghost = 1;
	double errors[2];
	for (int i = 0; i < 2; i++) {
		DomainReader<2>       domain_reader(mesh_file, {ns[i], ns[i]}, num_ghost);
		shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

		/*
		        auto ffun = [](const std::array<double, 2> &coord) {
		            double x = coord[0];
		            double y = coord[1];
		            return exp(x * y)
		                   * (pow(x, 5) * (-1 + y) * pow(y, 2) + y * (-2 + (5 - 3 * y) * y)
		                      - pow(x, 4) * y * (4 + (-7 + y) * y)
		                      - 2 * x * (1 + 2 * (-2 + y) * (-1 + y) * pow(y, 2))
		                      - pow(x, 2) * (-5 + 8 * y + (-6 + y) * (-1 + y) * pow(y, 3))
		                      + pow(x, 3) * (-3 + y * (12 - 6 * y - pow(y, 3) + pow(y, 4))));
		        };
		        auto gfun = [](const std::array<double, 2> &coord) {
		            double x = coord[0];
		            double y = coord[1];
		            return (1 - x) * x * (1 - y) * y * exp(x * y);
		        };
		        auto hfun = [](const std::array<double, 2> &coord) {
		            double x = coord[0];
		            double y = coord[1];
		            return 1 + x * y;
		        };
		        */
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
		auto hfun = [](const std::array<double, 2> &coord) { return 1; };

		auto f_vec          = ValVector<2>::GetNewVector(d_fine);
		auto f_vec_expected = ValVector<2>::GetNewVector(d_fine);
		DomainTools<2>::setValues(d_fine, f_vec_expected, ffun);
		Poisson::StarPatchOperator<2>::addDrichletBCToRHS(d_fine, f_vec_expected, gfun);

		auto g_vec = ValVector<2>::GetNewVector(d_fine);
		DomainTools<2>::setValues(d_fine, g_vec, gfun);

		auto h_vec = ValVector<2>::GetNewVector(d_fine);
		DomainTools<2>::setValuesWithGhost(d_fine, h_vec, hfun);

		auto gf         = make_shared<BiLinearGhostFiller>(d_fine);
		auto p_operator = make_shared<Poisson::StarPatchOperator<2>>(d_fine, gf);
		p_operator->apply(g_vec, f_vec);

		auto error_vec = ValVector<2>::GetNewVector(d_fine);
		error_vec->addScaled(1.0, f_vec, -1.0, f_vec_expected);
		errors[i] = error_vec->twoNorm() / f_vec_expected->twoNorm();
	}
	INFO("Errors: " << errors[0] << ", " << errors[1]);
	CHECK(log(errors[0] / errors[1]) / log(2) > 1.8);
}
TEST_CASE("Test Poisson::StarPatchOperator constructor throws exception with no ghost cells",
          "[Poisson::StarPatchOperator]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int n         = 32;
	int num_ghost = 0;

	DomainReader<2>       domain_reader(mesh_file, {n, n}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto gf = make_shared<BiLinearGhostFiller>(d_fine);
	CHECK_THROWS_AS(make_shared<Poisson::StarPatchOperator<2>>(d_fine, gf),
	                Thunderegg::Poisson::StarPatchOperatorException);
}
TEST_CASE("Test Poisson::StarPatchOperator::Generator", "[Poisson::StarPatchOperator]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int n         = 32;
	int num_ghost = 1;

	DomainReader<2>       domain_reader(mesh_file, {n, n}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();

	auto gf            = make_shared<BiLinearGhostFiller>(d_fine);
	auto gf_gen        = BiLinearGhostFiller::Generator(gf);
	auto op            = make_shared<Poisson::StarPatchOperator<2>>(d_fine, gf);
	auto op_gen        = Poisson::StarPatchOperator<2>::Generator(op, gf_gen);
	auto finer_level   = make_shared<GMG::Level<2>>(d_fine, nullptr);
	auto coarser_level = make_shared<GMG::Level<2>>(d_coarse, nullptr);
	finer_level->setCoarser(coarser_level);
	finer_level->setOperator(op);
	coarser_level->setFiner(finer_level);
	auto coarser_op  = op_gen(coarser_level);
	auto coarser_op2 = op_gen(coarser_level);
}