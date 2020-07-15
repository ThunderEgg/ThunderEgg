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
#include <Thunderegg/PetscMatOp.h>
#include <Thunderegg/PetscVector.h>
#include <Thunderegg/Poisson/MatrixHelper2d.h>
#include <Thunderegg/Poisson/StarPatchOperator.h>
using namespace std;
using namespace Thunderegg;
#define MESHES                                                                                     \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";
TEST_CASE("Poisson::MatrixHelper2d gives equivalent operator to Poisson::StarPatchOperator",
          "[Poisson::MatrixHelper2d]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int                   n         = 32;
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {n, n}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto gfun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return sinl(M_PI * y) * cosl(2 * M_PI * x);
	};

	auto f_vec          = PetscVector<2>::GetNewVector(d_fine);
	auto f_vec_expected = PetscVector<2>::GetNewVector(d_fine);

	auto g_vec = PetscVector<2>::GetNewVector(d_fine);
	DomainTools<2>::setValues(d_fine, g_vec, gfun);

	auto gf         = make_shared<BiLinearGhostFiller>(d_fine);
	auto p_operator = make_shared<Poisson::StarPatchOperator<2>>(d_fine, gf);
	p_operator->apply(g_vec, f_vec_expected);

	// generate matrix with matrix_helper
	Poisson::MatrixHelper2d mh(d_fine);
	auto                    m_operator = make_shared<PetscMatOp<2>>(mh.formCRSMatrix());
	m_operator->apply(g_vec, f_vec);

	for (auto pinfo : d_fine->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		LocalData<2> f_vec_ld          = f_vec->getLocalData(pinfo->local_index);
		LocalData<2> f_vec_expected_ld = f_vec_expected->getLocalData(pinfo->local_index);
		nested_loop<2>(f_vec_ld.getStart(), f_vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			CHECK(f_vec_ld[coord] == Approx(f_vec_expected_ld[coord]));
		});
	}
}