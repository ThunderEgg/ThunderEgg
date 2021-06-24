#include "../utils/DomainReader.h"
#include "PatchSolver_MOCKS.h"
#include <ThunderEgg/Iterative/PatchSolver.h>

#include <list>
#include <sstream>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

constexpr auto single_mesh_file  = "mesh_inputs/2d_uniform_2x2_mpi1.json";
constexpr auto refined_mesh_file = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json";
constexpr auto cross_mesh_file   = "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json";

TEST_CASE("Iterative::PatchSolver passes vectors of a single patch length",
          "[Iterative::PatchSolver]")
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
	auto u = make_shared<Vector<2>>(*d_fine, num_components);
	auto f = make_shared<Vector<2>>(*d_fine, num_components);

	auto mgf = make_shared<MockGhostFiller<2>>();
	// the patch operator is just a 0.5I operator
	auto mpo = make_shared<MockPatchOperator<2>>(d_fine, mgf);
	auto ms  = make_shared<MockSolver<2>>(
    [](std::shared_ptr<VectorGenerator<2>> vg, std::shared_ptr<const Operator<2>> A,
       std::shared_ptr<Vector<2>> x, std::shared_ptr<const Vector<2>> b,
       std::shared_ptr<const Operator<2>> Mr) {
        CHECK(x->getNumLocalPatches() == 1);
        return 1;
    });

	Iterative::PatchSolver<2> bcgs_solver(ms, mpo);

	bcgs_solver.smooth(*f, *u);
}
TEST_CASE("Iterative::PatchSolver passes modified operator", "[Iterative::PatchSolver]")
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
	auto u = make_shared<Vector<2>>(*d_fine, num_components);
	auto f = make_shared<Vector<2>>(*d_fine, num_components);

	bool called = false;
	auto mgf    = make_shared<MockGhostFiller<2>>();
	// the patch operator is just a 0.5I operator
	auto mpo = make_shared<MockPatchOperator<2>>(d_fine, mgf);
	auto ms  = make_shared<MockSolver<2>>(
    [&](std::shared_ptr<VectorGenerator<2>> vg, std::shared_ptr<const Operator<2>> A,
        std::shared_ptr<Vector<2>> x, std::shared_ptr<const Vector<2>> b,
        std::shared_ptr<const Operator<2>> Mr) {
        if (!called) {
            called = true;
            CHECK(A != mpo);
            A->apply(*b, *x);
        }
        return 1;
    });

	Iterative::PatchSolver<2> bcgs_solver(ms, mpo);

	bcgs_solver.smooth(*f, *u);
	CHECK(mpo->getNumApplyCalls() == 1);
	CHECK(mpo->rhsWasModified());
	CHECK(mpo->boundaryConditionsEnforced());
	CHECK(mpo->internalBoundaryConditionsEnforced());
}
TEST_CASE("Iterative::PatchSolver propagates BreakdownError", "[Iterative::PatchSolver]")
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
	auto u = make_shared<Vector<2>>(*d_fine, num_components);
	auto f = make_shared<Vector<2>>(*d_fine, num_components);
	f->set(1);

	auto mgf = make_shared<MockGhostFiller<2>>();
	// the patch operator is just a 0.5I operator
	auto mpo = make_shared<NonLinMockPatchOperator<2>>(d_fine, mgf);
	auto ms  = make_shared<MockSolver<2>>(
    [](std::shared_ptr<VectorGenerator<2>> vg, std::shared_ptr<const Operator<2>> A,
       std::shared_ptr<Vector<2>> x, std::shared_ptr<const Vector<2>> b,
       std::shared_ptr<const Operator<2>> Mr) {
        throw BreakdownError("Blah");
        return 1;
    });

	Iterative::PatchSolver<2> bcgs_solver(ms, mpo);

	CHECK_THROWS_AS(bcgs_solver.smooth(*f, *u), BreakdownError);
}
TEST_CASE("Iterative::PatchSolver does not propagate BreakdownError", "[Iterative::PatchSolver]")
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
	auto u = make_shared<Vector<2>>(*d_fine, num_components);
	auto f = make_shared<Vector<2>>(*d_fine, num_components);
	f->set(1);

	auto mgf = make_shared<MockGhostFiller<2>>();
	// the patch operator is just a 0.5I operator
	auto mpo = make_shared<NonLinMockPatchOperator<2>>(d_fine, mgf);
	auto ms  = make_shared<MockSolver<2>>(
    [](std::shared_ptr<VectorGenerator<2>> vg, std::shared_ptr<const Operator<2>> A,
       std::shared_ptr<Vector<2>> x, std::shared_ptr<const Vector<2>> b,
       std::shared_ptr<const Operator<2>> Mr) {
        throw BreakdownError("Blah");
        return 1;
    });

	Iterative::PatchSolver<2> bcgs_solver(ms, mpo, true);

	CHECK_NOTHROW(bcgs_solver.smooth(*f, *u));
}