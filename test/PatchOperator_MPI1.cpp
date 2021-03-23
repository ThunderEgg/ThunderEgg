#include "PatchOperator_MOCKS.h"
#include "utils/DomainReader.h"
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/MPIGhostFiller.h>
#include <ThunderEgg/ValVector.h>

#include <list>

#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

constexpr auto single_mesh_file  = "mesh_inputs/2d_uniform_2x2_mpi1.json";
constexpr auto refined_mesh_file = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json";
constexpr auto cross_mesh_file   = "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json";

TEST_CASE("Check PatchOperator calls for various domains", "[PatchOperator]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 5);
	auto                  ny        = GENERATE(2, 5);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto u_num_components = GENERATE(1, 2, 3);
	auto u                = ValVector<2>::GetNewVector(d_fine, u_num_components);
	auto f_num_components = GENERATE(1, 2, 3);
	auto f                = ValVector<2>::GetNewVector(d_fine, f_num_components);

	auto                 mgf = make_shared<MockGhostFiller<2>>();
	MockPatchOperator<2> mpo(d_fine, mgf, u, f);

	mpo.apply(u, f);

	CHECK(mgf->wasCalled());
	CHECK(mpo.allPatchesCalled());
}
TEST_CASE("PatchOperator check getDomain", "[PatchOperator]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto                  num_components = GENERATE(1, 2, 3);
	auto                  nx             = GENERATE(2, 5);
	auto                  ny             = GENERATE(2, 5);
	int                   num_ghost      = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto u = ValVector<2>::GetNewVector(d_fine, num_components);
	auto f = ValVector<2>::GetNewVector(d_fine, num_components);

	auto                 mgf = make_shared<MockGhostFiller<2>>();
	MockPatchOperator<2> mpo(d_fine, mgf, u, f);

	CHECK(mpo.getDomain() == d_fine);
}
TEST_CASE("PatchOperator check getGhostFiller", "[PatchOperator]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto                  num_components = GENERATE(1, 2, 3);
	auto                  nx             = GENERATE(2, 5);
	auto                  ny             = GENERATE(2, 5);
	int                   num_ghost      = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto u = ValVector<2>::GetNewVector(d_fine, num_components);
	auto f = ValVector<2>::GetNewVector(d_fine, num_components);

	auto                 mgf = make_shared<MockGhostFiller<2>>();
	MockPatchOperator<2> mpo(d_fine, mgf, u, f);

	CHECK(mpo.getGhostFiller() == mgf);
}