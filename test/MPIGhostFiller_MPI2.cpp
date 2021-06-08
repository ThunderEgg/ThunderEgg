#include "MPIGhostFiller_MOCKS.h"
#include "utils/DomainReader.h"
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/MPIGhostFiller.h>
#include <ThunderEgg/ValVector.h>

#include <list>

#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

constexpr auto uniform    = "mesh_inputs/2d_uniform_2x2_nw_on_1_mpi2.json";
constexpr auto refined    = "mesh_inputs/2d_uniform_2x2_refined_nw_on_1_mpi2.json";
constexpr auto td_refined = "mesh_inputs/3d_refined_bnw_2x2x2_mpi2.json";

TEST_CASE("Calls for various domains 2d face cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, uniform, refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<2>::GetNewVector(d_fine, num_components);

	CallMockMPIGhostFiller<2> mgf(d_fine, num_components, GhostFillingType::Faces);

	mgf.fillGhost(vec);

	CHECK(mgf.called == true);

	mgf.checkCalls();
}
TEST_CASE("Calls for various domains 2d corner cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, uniform, refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<2>::GetNewVector(d_fine, num_components);

	CallMockMPIGhostFiller<2> mgf(d_fine, num_components, GhostFillingType::Corners);

	mgf.fillGhost(vec);

	CHECK(mgf.called == true);

	mgf.checkCalls();
}
TEST_CASE("Calls for various domains 3d face cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, td_refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	auto                  nz        = GENERATE(2, 3);
	int                   num_ghost = 1;
	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<3>::GetNewVector(d_fine, num_components);

	CallMockMPIGhostFiller<3> mgf(d_fine, num_components, GhostFillingType::Faces);

	mgf.fillGhost(vec);

	CHECK(mgf.called == true);

	mgf.checkCalls();
}
TEST_CASE("Calls for various domains 3d edge cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, td_refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	auto                  nz        = GENERATE(2, 3);
	int                   num_ghost = 1;
	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<3>::GetNewVector(d_fine, num_components);

	CallMockMPIGhostFiller<3> mgf(d_fine, num_components, GhostFillingType::Edges);

	mgf.fillGhost(vec);

	CHECK(mgf.called == true);

	mgf.checkCalls();
}
TEST_CASE("Calls for various domains 3d corners cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, td_refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	auto                  nz        = GENERATE(2, 3);
	int                   num_ghost = 1;
	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<3>::GetNewVector(d_fine, num_components);

	CallMockMPIGhostFiller<3> mgf(d_fine, num_components, GhostFillingType::Corners);

	mgf.fillGhost(vec);

	CHECK(mgf.called == true);

	mgf.checkCalls();
}
TEST_CASE("Exchange for various domains 2d face cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, uniform, refined);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	int  num_ghost = GENERATE(1, 2);
	INFO("nx: " << nx);
	INFO("ny: " << ny);
	INFO("num_ghost: " << ny);
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<2>::GetNewVector(d_fine, num_components);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		for (int c = 0; c < num_components; c++) {
			auto data = vec->getComponentView(c, pinfo.local_index);
			nested_loop<2>(data.getStart(), data.getEnd(),
			               [&](const std::array<int, 2> &coord) { data[coord] = pinfo.id; });
		}
	}

	ExchangeMockMPIGhostFiller<2> mgf(d_fine, GhostFillingType::Faces);

	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}
TEST_CASE("Exchange for various domains 2d corner cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, uniform, refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	int                   num_ghost = GENERATE(1, 2);
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<2>::GetNewVector(d_fine, num_components);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		for (int c = 0; c < num_components; c++) {
			auto data = vec->getComponentView(c, pinfo.local_index);
			nested_loop<2>(data.getStart(), data.getEnd(),
			               [&](const std::array<int, 2> &coord) { data[coord] = pinfo.id; });
		}
	}

	ExchangeMockMPIGhostFiller<2> mgf(d_fine, GhostFillingType::Corners);

	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}
TEST_CASE("Exchange for various domains 3d face cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, td_refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	auto                  nz        = GENERATE(2, 3);
	int                   num_ghost = GENERATE(1, 2);
	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<3>::GetNewVector(d_fine, num_components);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		for (int c = 0; c < num_components; c++) {
			auto data = vec->getComponentView(c, pinfo.local_index);
			nested_loop<3>(data.getStart(), data.getEnd(),
			               [&](const std::array<int, 3> &coord) { data[coord] = pinfo.id; });
		}
	}

	ExchangeMockMPIGhostFiller<3> mgf(d_fine, GhostFillingType::Faces);

	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}
TEST_CASE("Exchange for various domains 3d edge cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, td_refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	auto                  nz        = GENERATE(2, 3);
	int                   num_ghost = GENERATE(1, 2);
	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<3>::GetNewVector(d_fine, num_components);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		for (int c = 0; c < num_components; c++) {
			auto data = vec->getComponentView(c, pinfo.local_index);
			nested_loop<3>(data.getStart(), data.getEnd(),
			               [&](const std::array<int, 3> &coord) { data[coord] = pinfo.id; });
		}
	}

	ExchangeMockMPIGhostFiller<3> mgf(d_fine, GhostFillingType::Edges);

	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}
TEST_CASE("Exchange for various domains 3d corner cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, td_refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	auto                  nz        = GENERATE(2, 3);
	int                   num_ghost = GENERATE(1, 2);
	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<3>::GetNewVector(d_fine, num_components);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		for (int c = 0; c < num_components; c++) {
			auto data = vec->getComponentView(c, pinfo.local_index);
			nested_loop<3>(data.getStart(), data.getEnd(),
			               [&](const std::array<int, 3> &coord) { data[coord] = pinfo.id; });
		}
	}

	ExchangeMockMPIGhostFiller<3> mgf(d_fine, GhostFillingType::Corners);

	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}
TEST_CASE("Two Exchanges for various domains 2d face cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, uniform, refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	int                   num_ghost = GENERATE(1, 2);
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<2>::GetNewVector(d_fine, num_components);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		for (int c = 0; c < num_components; c++) {
			auto data = vec->getComponentView(c, pinfo.local_index);
			nested_loop<2>(data.getStart(), data.getEnd(),
			               [&](const std::array<int, 2> &coord) { data[coord] = pinfo.id; });
		}
	}

	ExchangeMockMPIGhostFiller<2> mgf(d_fine, GhostFillingType::Faces);

	mgf.fillGhost(vec);
	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}
TEST_CASE("Two Exchanges for various domains 2d corner cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, uniform, refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	int                   num_ghost = GENERATE(1, 2);
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<2>::GetNewVector(d_fine, num_components);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		for (int c = 0; c < num_components; c++) {
			auto data = vec->getComponentView(c, pinfo.local_index);
			nested_loop<2>(data.getStart(), data.getEnd(),
			               [&](const std::array<int, 2> &coord) { data[coord] = pinfo.id; });
		}
	}

	ExchangeMockMPIGhostFiller<2> mgf(d_fine, GhostFillingType::Corners);

	mgf.fillGhost(vec);
	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}
TEST_CASE("Two Exchange for various domains 3d face cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, td_refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	auto                  nz        = GENERATE(2, 3);
	int                   num_ghost = GENERATE(1, 2);
	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<3>::GetNewVector(d_fine, num_components);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		for (int c = 0; c < num_components; c++) {
			auto data = vec->getComponentView(c, pinfo.local_index);
			nested_loop<3>(data.getStart(), data.getEnd(),
			               [&](const std::array<int, 3> &coord) { data[coord] = pinfo.id; });
		}
	}

	ExchangeMockMPIGhostFiller<3> mgf(d_fine, GhostFillingType::Faces);

	mgf.fillGhost(vec);
	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}
TEST_CASE("Two Exchange for various domains 3d edge cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, td_refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	auto                  nz        = GENERATE(2, 3);
	int                   num_ghost = GENERATE(1, 2);
	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<3>::GetNewVector(d_fine, num_components);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		for (int c = 0; c < num_components; c++) {
			auto data = vec->getComponentView(c, pinfo.local_index);
			nested_loop<3>(data.getStart(), data.getEnd(),
			               [&](const std::array<int, 3> &coord) { data[coord] = pinfo.id; });
		}
	}

	ExchangeMockMPIGhostFiller<3> mgf(d_fine, GhostFillingType::Edges);

	mgf.fillGhost(vec);
	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}
TEST_CASE("Two Exchange for various domains 3d corner cases MPI2", "[MPIGhostFiller]")
{
	auto num_components = GENERATE(1, 2);
	auto mesh_file      = GENERATE(as<std::string>{}, td_refined);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2, 3);
	auto                  ny        = GENERATE(2, 3);
	auto                  nz        = GENERATE(2, 3);
	int                   num_ghost = GENERATE(1, 2);
	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<3>::GetNewVector(d_fine, num_components);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		for (int c = 0; c < num_components; c++) {
			auto data = vec->getComponentView(c, pinfo.local_index);
			nested_loop<3>(data.getStart(), data.getEnd(),
			               [&](const std::array<int, 3> &coord) { data[coord] = pinfo.id; });
		}
	}

	ExchangeMockMPIGhostFiller<3> mgf(d_fine, GhostFillingType::Corners);

	mgf.fillGhost(vec);
	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}