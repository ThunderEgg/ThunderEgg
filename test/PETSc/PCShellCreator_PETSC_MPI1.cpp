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
#include "catch.hpp"
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/PETSc/PCShellCreator.h>
#include <ThunderEgg/PETSc/VecWrapper.h>
#include <ThunderEgg/ValVectorGenerator.h>
using namespace std;
using namespace ThunderEgg;
#define MESHES                                                                                     \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";

// mock operator
class HalfIdentity : public Operator<2>
{
	public:
	void apply(std::shared_ptr<const Vector<2>> x, std::shared_ptr<Vector<2>> b) const override
	{
		b->copy(x);
		b->scale(0.5);
	}
};
TEST_CASE("PETSc::PCShellCreator works with 0.5I", "[PETSc::PCShellCreator]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int                   n         = 32;
	int                   num_ghost = 0;
	DomainReader<2>       domain_reader(mesh_file, {n, n}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();
	auto                  vg     = make_shared<ValVectorGenerator<2>>(d_fine);

	auto gfun = [](const std::array<double, 2> &coord) {
		double x = coord[0];
		double y = coord[1];
		return sinl(M_PI * y) * cosl(2 * M_PI * x);
	};

	auto x = PETSc::VecWrapper<2>::GetNewVector(d_fine);
	DomainTools<2>::setValues(d_fine, x, gfun);
	auto b = PETSc::VecWrapper<2>::GetNewVector(d_fine);

	// create an Identity matrix
	auto TE_A = make_shared<HalfIdentity>();
	PC   P    = PETSc::PCShellCreator<2>::GetNewPCShell(TE_A, TE_A, vg);

	PCApply(P, x->getVec(), b->getVec());

	for (auto pinfo : d_fine->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		INFO("dx:    " << pinfo->spacings[0]);
		INFO("dy:    " << pinfo->spacings[1]);
		LocalData<2> x_ld = x->getLocalData(pinfo->local_index);
		LocalData<2> b_ld = b->getLocalData(pinfo->local_index);
		nested_loop<2>(x_ld.getStart(), x_ld.getEnd(), [&](const array<int, 2> &coord) {
			INFO("xi:    " << coord[0]);
			INFO("yi:    " << coord[1]);
			CHECK(0.5 * x_ld[coord] == Approx(b_ld[coord]));
		});
	}
	PCDestroy(&P);
}