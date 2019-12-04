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

#include "catch.hpp"
#include <Thunderegg/BiLinearGhostFiller.h>
#include <Thunderegg/DomainTools.h>
#include <Thunderegg/Experimental/DomGen.h>
#include <Thunderegg/Experimental/OctTree.h>
#include <Thunderegg/GMG/LinearRestrictor.h>
#include <Thunderegg/ValVector.h>
using namespace std;
using namespace Thunderegg;
TEST_CASE("Linear Test LinearRestrictor", "[GMG::LinearRestrictor]")
{
	auto nx        = GENERATE(2, 10);
	auto ny        = GENERATE(2, 10);
	int  num_ghost = 1;

	Experimental::Tree<2>   t("middle_refine.bin");
	Experimental::DomGen<2> dg(t, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>>   d        = dg.getFinestDomain();
	shared_ptr<Domain<2>>   d_coarse = dg.getCoarserDomain();

	auto fine_vec = PetscVector<2>::GetNewVector(d);
	auto vec      = PetscVector<2>::GetNewVector(d_coarse);
	auto expected = PetscVector<2>::GetNewVector(d_coarse);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 1 + ((x * 0.3) + y);
	};

	DomainTools<2>::setValuesWithGhost(d, fine_vec, f);
	DomainTools<2>::setValuesWithGhost(d_coarse, expected, f);

	auto ilc        = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d);
	auto restrictor = std::make_shared<GMG::LinearRestrictor<2>>(d_coarse, d, ilc);

	restrictor->restrict(vec, fine_vec);

	for (auto pinfo : d_coarse->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		LocalData<2> vec_ld      = vec->getLocalData(pinfo->local_index);
		LocalData<2> expected_ld = expected->getLocalData(pinfo->local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			REQUIRE(vec_ld[coord] == Approx(expected_ld[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			if (pinfo->hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo->getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 1> &coord) {
					               INFO("coord:  " << coord[0]);
					               CHECK(vec_ghost[coord] == Approx(expected_ghost[coord]));
				               });
			}
		}
	}
}