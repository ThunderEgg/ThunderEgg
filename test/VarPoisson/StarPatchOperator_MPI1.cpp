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
#include <Thunderegg/ValVector.h>
#include <Thunderegg/VarPoisson/StarPatchOperator.h>
using namespace std;
using namespace Thunderegg;
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";
TEST_CASE("Test StarPatchOperator add ghost to RHS", "[GMG::StarPatchOperator]")
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
	/*
	for (auto pinfo : d_fine->getPatchInfoVector()) {
	    auto ld = g_zero->getLocalData(0);
	    for (Side<2> s : Side<2>::getValues()) {
	        if (!pinfo->hasNbr(s)) {
	            auto ghosts = ld.getGhostSliceOnSide(s, 1);
	            nested_loop<1>(ghosts.getStart(), ghosts.getEnd(),
	                           [&](const array<int, 1> &coord) { ghosts[coord] = 0; });
	        }
	    }
	}
	*/

	auto h_vec = ValVector<2>::GetNewVector(d_fine);
	DomainTools<2>::setValuesWithGhost(d_fine, h_vec, hfun);

	shared_ptr<BiLinearGhostFiller>              gf(new BiLinearGhostFiller(d_fine));
	shared_ptr<VarPoisson::StarPatchOperator<2>> p_operator(
	new VarPoisson::StarPatchOperator<2>(h_vec, d_fine, gf));
	VarPoisson::StarPatchOperator<2>::addDrichletBCToRHS(d_fine, f_vec, gfun, hfun);

	auto f_expected = ValVector<2>::GetNewVector(d_fine);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		p_operator->applySinglePatch(pinfo, g_zero->getLocalData(pinfo->local_index),
		                             f_expected->getLocalData(pinfo->local_index));
	}
	f_expected->scaleThenAdd(-1.0, f_vec);

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