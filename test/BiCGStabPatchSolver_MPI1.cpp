#include "BiCGStabPatchSolver_MOCKS.h"
#include "catch.hpp"
#include "utils/DomainReader.h"
#include <ThunderEgg/BiCGStabPatchSolver.h>
#include <ThunderEgg/ValVector.h>
#include <list>
#include <sstream>
using namespace std;
using namespace ThunderEgg;

constexpr auto single_mesh_file  = "mesh_inputs/2d_uniform_2x2_mpi1.json";
constexpr auto refined_mesh_file = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json";
constexpr auto cross_mesh_file   = "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json";

TEST_CASE("BiCGStabPatchSolver smooth", "[BiCGStabPatchSolver]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 5);
	auto                  ny        = GENERATE(2, 5);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto num_components = GENERATE(1, 2, 3);
	INFO("num_components: " << num_components);
	auto u_expected = ValVector<2>::GetNewVector(d_fine, num_components);
	for (int i = 0; i < u_expected->getNumLocalPatches(); i++) {
		for (int c = 0; c < u_expected->getNumComponents(); c++) {
			auto ld = u_expected->getLocalData(c, i);
			nested_loop<2>(ld.getStart(), ld.getEnd(),
			               [&](const std::array<int, 2> &coord) { ld[coord] = i + c + 1; });
		}
	}
	auto u = ValVector<2>::GetNewVector(d_fine, num_components);
	auto f = ValVector<2>::GetNewVector(d_fine, num_components);

	auto mgf = make_shared<MockGhostFiller<2>>();
	// the patch operator is just a 0.5I operator
	auto mpo = make_shared<MockPatchOperator<2>>(d_fine, mgf);

	mpo->apply(u_expected, f);

	BiCGStabPatchSolver<2> bcgs_solver(mpo);

	bcgs_solver.smooth(f, u);

	for (int i = 0; i < u->getNumLocalPatches(); i++) {
		INFO("PATCH_INDEX: " << i);
		for (int c = 0; c < u->getNumComponents(); c++) {
			INFO("c: " << c);
			auto ld          = u->getLocalData(c, i);
			auto ld_expected = u_expected->getLocalData(c, i);
			nested_loop<2>(ld.getStart(), ld.getEnd(), [&](const std::array<int, 2> &coord) {
				CHECK(ld[coord] == Approx(ld_expected[coord]).epsilon(1e-8));
			});
		}
	}
	CHECK(mgf->numCalls() == 2);
	CHECK(mpo->rhsWasModified());
	CHECK(mpo->interiorDirichlet());
}