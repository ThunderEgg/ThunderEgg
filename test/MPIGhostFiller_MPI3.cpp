#include "MPIGhostFiller_MOCKS.h"
#include "catch.hpp"
#include "utils/DomainReader.h"
#include <Thunderegg/DomainTools.h>
#include <Thunderegg/MPIGhostFiller.h>
#include <Thunderegg/ValVector.h>
#include <list>
using namespace std;
using namespace Thunderegg;

constexpr auto uniform = "mesh_inputs/2d_uniform_2x2_mpi3.json";
constexpr auto refined = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi3.json";

TEST_CASE("Calls for various domains 1-side cases MPI3", "[MPIGhostFiller]")
{
	auto mesh_file = GENERATE(as<std::string>{}, uniform, refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2);
	auto                  ny        = GENERATE(2);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<2>::GetNewVector(d_fine);

	CallMockMPIGhostFiller<2> mgf(d_fine, 1);

	mgf.fillGhost(vec);

	CHECK(mgf.called == true);

	mgf.checkCalls();
}
TEST_CASE("Exchange for various domains 1-side cases MPI3", "[MPIGhostFiller]")
{
	auto mesh_file = GENERATE(as<std::string>{}, uniform, refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2);
	auto                  ny        = GENERATE(2);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<2>::GetNewVector(d_fine);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		auto data = vec->getLocalData(pinfo->local_index);
		nested_loop<2>(data.getStart(), data.getEnd(),
		               [&](const std::array<int, 2> &coord) { data[coord] = pinfo->id; });
	}

	ExchangeMockMPIGhostFiller<2> mgf(d_fine, 1);

	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}
TEST_CASE("Two Exchanges for various domains 1-side cases MPI3", "[MPIGhostFiller]")
{
	auto mesh_file = GENERATE(as<std::string>{}, uniform, refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2);
	auto                  ny        = GENERATE(2);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<2>::GetNewVector(d_fine);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		auto data = vec->getLocalData(pinfo->local_index);
		nested_loop<2>(data.getStart(), data.getEnd(),
		               [&](const std::array<int, 2> &coord) { data[coord] = pinfo->id; });
	}

	ExchangeMockMPIGhostFiller<2> mgf(d_fine, 1);

	mgf.fillGhost(vec);
	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}