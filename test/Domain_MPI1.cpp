#include <ThunderEgg/Domain.h>

#include "utils/DomainReader.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Domain constructors work", "[Domain]")
{
	Communicator comm(MPI_COMM_WORLD);

	vector<PatchInfo<2>> pinfos(1);

	auto n         = GENERATE(1, 2, 10, 13);
	auto spacing   = GENERATE(0.01, 1.0, 3.14);
	auto num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].id = 0;
	pinfos[0].ns.fill(n);
	pinfos[0].spacings.fill(spacing);
	pinfos[0].num_ghost_cells = num_ghost;
	Domain<2> d(comm, 0, {n, n}, num_ghost, pinfos.begin(), pinfos.end());

	// check getters
	for (int ni : d.getNs()) {
		CHECK(ni == n);
	}
	CHECK(d.getNumGlobalPatches() == 1);
	CHECK(d.getNumLocalPatches() == 1);
	CHECK(d.getNumGlobalCells() == n * n);
	CHECK(d.getNumLocalCells() == n * n);
	CHECK(d.getNumLocalCellsWithGhost() == (n + 2 * num_ghost) * (n + 2 * num_ghost));
	CHECK(d.getNumCellsInPatch() == n * n);
	CHECK(d.getNumGhostCells() == num_ghost);
	CHECK(d.volume() == Catch::Approx(spacing * spacing * n * n));

	int result;
	int err = MPI_Comm_compare(comm.getMPIComm(), d.getCommunicator().getMPIComm(), &result);
	REQUIRE(err == MPI_SUCCESS);
	CHECK(result == MPI_CONGRUENT);
}
TEST_CASE("Domain setTimer", "[Domain]")
{
	Communicator comm(MPI_COMM_WORLD);

	vector<PatchInfo<2>> pinfos(1);

	int    n         = 10;
	double spacing   = 0.01;
	int    num_ghost = 1;

	pinfos[0].id = 0;
	pinfos[0].ns.fill(n);
	pinfos[0].spacings.fill(spacing);
	pinfos[0].num_ghost_cells = num_ghost;
	Domain<2> d(comm, 0, {n, n}, num_ghost, pinfos.begin(), pinfos.end());

	auto timer = make_shared<Timer>(comm);
	d.setTimer(timer);
	CHECK(d.getTimer() == timer);
	CHECK(d.hasTimer());
}
TEST_CASE("Domain setTimer adds domain to timer", "[Domain]")
{
	Communicator comm(MPI_COMM_WORLD);

	vector<PatchInfo<2>> pinfos(1);

	int    n         = 10;
	double spacing   = 0.01;
	int    num_ghost = 1;

	pinfos[0].id = 0;
	pinfos[0].ns.fill(n);
	pinfos[0].spacings.fill(spacing);
	pinfos[0].num_ghost_cells = num_ghost;
	Domain<2> d(comm, 0, {n, n}, num_ghost, pinfos.begin(), pinfos.end());

	auto timer = make_shared<Timer>(comm);
	d.setTimer(timer);
	// will throw exception if domain not added to timer
	timer->startDomainTiming(0, "A");
	timer->stopDomainTiming(0, "A");
}
TEST_CASE("Domain getTimer default is no timer", "[Domain]")
{
	Communicator comm(MPI_COMM_WORLD);

	vector<PatchInfo<2>> pinfos(1);

	int    n         = 10;
	double spacing   = 0.01;
	int    num_ghost = 1;

	pinfos[0].id = 0;
	pinfos[0].ns.fill(n);
	pinfos[0].spacings.fill(spacing);
	pinfos[0].num_ghost_cells = num_ghost;
	Domain<2> d(comm, 0, {n, n}, num_ghost, pinfos.begin(), pinfos.end());

	CHECK(d.getTimer() == nullptr);
	CHECK_FALSE(d.hasTimer());
}
TEST_CASE("Domain id in constructor", "[Domain]")
{
	Communicator comm(MPI_COMM_WORLD);

	vector<PatchInfo<2>> pinfos(1);

	int    n         = 10;
	double spacing   = 0.01;
	int    num_ghost = 1;

	pinfos[0].id = 0;
	pinfos[0].ns.fill(n);
	pinfos[0].spacings.fill(spacing);
	pinfos[0].num_ghost_cells = num_ghost;
	auto      id              = GENERATE(1, 2, 9);
	Domain<2> d(comm, id, {n, n}, num_ghost, pinfos.begin(), pinfos.end());

	CHECK(d.getId() == id);
}
TEST_CASE("Domain to_json", "[Domain]")
{
	Communicator comm(MPI_COMM_WORLD);

	vector<PatchInfo<2>> pinfos(1);

	int    n         = 10;
	double spacing   = 0.01;
	int    num_ghost = 1;

	pinfos[0].id = 0;
	pinfos[0].ns.fill(n);
	pinfos[0].spacings.fill(spacing);
	pinfos[0].num_ghost_cells = num_ghost;
	Domain<2> d(comm, 0, {n, n}, num_ghost, pinfos.begin(), pinfos.end());

	nlohmann::json j = d;
	REQUIRE(j.is_array());
	REQUIRE(j.size() == 1);
	REQUIRE(j[0]["id"] == 0);
}
TEST_CASE("Domain<2> numLocalPatches",
          "[Schur::InterfaceDomain]")
{
	DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", {10, 10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	CHECK(domain->getNumLocalPatches() == 15);
}
TEST_CASE("Domain<2> numGlobalPatches",
          "[Schur::InterfaceDomain]")
{
	DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", {10, 10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();

	CHECK(domain->getNumGlobalPatches() == 15);
}
TEST_CASE("Domain<2> numLocalCells",
          "[Schur::InterfaceDomain]")
{
	DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", {10, 10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();
	CHECK(domain->getNumLocalCells() == 15 * 10 * 10 * 10);
}
TEST_CASE("Domain<2> numGlobalCells",
          "[Schur::InterfaceDomain]")
{
	DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", {10, 10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();

	CHECK(domain->getNumGlobalCells() == 15 * 10 * 10 * 10);
}
TEST_CASE("Domain<2> getNumGhostCells",
          "[Schur::InterfaceDomain]")
{
	DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", {10, 10, 10}, 2);
	auto            domain = domain_reader.getFinerDomain();
	CHECK(domain->getNumGhostCells() == 2);
}
TEST_CASE("Domain<2> numLocalCellsWithGhost",
          "[Schur::InterfaceDomain]")
{
	DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", {10, 10, 10}, 2);
	auto            domain = domain_reader.getFinerDomain();

	CHECK(domain->getNumLocalCellsWithGhost() == 15 * 14 * 14 * 14);
}
TEST_CASE("Domain<2> getNs",
          "[Schur::InterfaceDomain]")
{
	DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", {10, 11, 12}, 2);
	auto            domain = domain_reader.getFinerDomain();

	CHECK(domain->getNs()[0] == 10);
	CHECK(domain->getNs()[1] == 11);
	CHECK(domain->getNs()[2] == 12);
}
TEST_CASE("Domain<2> getPatchInfoVector size",
          "[Schur::InterfaceDomain]")
{
	DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", {10, 10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();

	CHECK(domain->getPatchInfoVector().size() == 15);
}
TEST_CASE("Domain<3> local indexes match position in pinfo vector",
          "[Domain]")
{
	DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", {10, 10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();

	auto pinfo_vector = domain->getPatchInfoVector();
	for (size_t i = 0; i < pinfo_vector.size(); i++) {
		CHECK(pinfo_vector[i].local_index == (int) i);
	}
}
TEST_CASE("Schur::InterfaceDomain<2> local indexes in neighbor info are consistent",
          "[Schur::InterfaceDomain]")
{
	DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", {10, 10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();

	map<int, int> id_to_local_index_map;

	for (auto pinfo : domain->getPatchInfoVector()) {
		id_to_local_index_map[pinfo.id] = pinfo.local_index;
	}

	auto checkIdAndLocalIndex = [&](int id, int local_index) {
		auto iter = id_to_local_index_map.find(id);
		if (iter != id_to_local_index_map.end()) {
			CHECK(local_index == iter->second);
		} else {
			CHECK(local_index == -1);
		}
	};
	for (auto pinfo : domain->getPatchInfoVector()) {
		for (Side<3> s : Side<3>::getValues()) {
			if (pinfo.hasNbr(s)) {
				NbrType type = pinfo.getNbrType(s);
				if (type == NbrType::Normal) {
					const NormalNbrInfo<2> &info = pinfo.getNormalNbrInfo(s);
					checkIdAndLocalIndex(info.id, info.local_index);
				} else if (type == NbrType::Fine) {
					const FineNbrInfo<2> &info = pinfo.getFineNbrInfo(s);
					for (size_t i = 0; i < info.ids.size(); i++) {
						checkIdAndLocalIndex(info.ids[i], info.local_indexes[i]);
					}
				} else if (type == NbrType::Coarse) {
					const CoarseNbrInfo<2> &info = pinfo.getCoarseNbrInfo(s);
					checkIdAndLocalIndex(info.id, info.local_index);
				}
			}
		}
		for (Edge e : Edge::getValues()) {
			if (pinfo.hasNbr(e)) {
				NbrType type = pinfo.getNbrType(e);
				if (type == NbrType::Normal) {
					const NormalNbrInfo<1> &info = pinfo.getNormalNbrInfo(e);
					checkIdAndLocalIndex(info.id, info.local_index);
				} else if (type == NbrType::Fine) {
					const FineNbrInfo<1> &info = pinfo.getFineNbrInfo(e);
					for (size_t i = 0; i < info.ids.size(); i++) {
						checkIdAndLocalIndex(info.ids[i], info.local_indexes[i]);
					}
				} else if (type == NbrType::Coarse) {
					const CoarseNbrInfo<1> &info = pinfo.getCoarseNbrInfo(e);
					checkIdAndLocalIndex(info.id, info.local_index);
				}
			}
		}
		for (Corner<3> c : Corner<3>::getValues()) {
			if (pinfo.hasNbr(c)) {
				NbrType type = pinfo.getNbrType(c);
				if (type == NbrType::Normal) {
					const NormalNbrInfo<0> &info = pinfo.getNormalNbrInfo(c);
					checkIdAndLocalIndex(info.id, info.local_index);
				} else if (type == NbrType::Fine) {
					const FineNbrInfo<0> &info = pinfo.getFineNbrInfo(c);
					for (size_t i = 0; i < info.ids.size(); i++) {
						checkIdAndLocalIndex(info.ids[i], info.local_indexes[i]);
					}
				} else if (type == NbrType::Coarse) {
					const CoarseNbrInfo<0> &info = pinfo.getCoarseNbrInfo(c);
					checkIdAndLocalIndex(info.id, info.local_index);
				}
			}
		}
	}
}
TEST_CASE("Domain<3> global indexes match position in pinfo vector",
          "[Domain]")
{
	DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", {10, 10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();

	auto pinfo_vector = domain->getPatchInfoVector();
	for (size_t i = 0; i < pinfo_vector.size(); i++) {
		CHECK(pinfo_vector[i].global_index == (int) i);
	}
}
TEST_CASE("Domain<3> global indexes in neighbor info are consistent",
          "[Domain]")
{
	DomainReader<3> domain_reader("mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", {10, 10, 10}, 0);
	auto            domain = domain_reader.getFinerDomain();

	map<int, int> id_to_global_index_map;

	for (auto pinfo : domain->getPatchInfoVector()) {
		id_to_global_index_map[pinfo.id] = pinfo.local_index;
	}

	auto checkIdAndLocalIndex = [&](int id, int global_index) {
		CHECK(global_index == id_to_global_index_map.at(id));
	};
	for (auto pinfo : domain->getPatchInfoVector()) {
		for (Side<3> s : Side<3>::getValues()) {
			if (pinfo.hasNbr(s)) {
				NbrType type = pinfo.getNbrType(s);
				if (type == NbrType::Normal) {
					const NormalNbrInfo<2> &info = pinfo.getNormalNbrInfo(s);
					checkIdAndLocalIndex(info.id, info.global_index);
				} else if (type == NbrType::Fine) {
					const FineNbrInfo<2> &info = pinfo.getFineNbrInfo(s);
					for (size_t i = 0; i < info.ids.size(); i++) {
						checkIdAndLocalIndex(info.ids[i], info.global_indexes[i]);
					}
				} else if (type == NbrType::Coarse) {
					const CoarseNbrInfo<2> &info = pinfo.getCoarseNbrInfo(s);
					checkIdAndLocalIndex(info.id, info.global_index);
				}
			}
		}
		for (Edge e : Edge::getValues()) {
			if (pinfo.hasNbr(e)) {
				NbrType type = pinfo.getNbrType(e);
				if (type == NbrType::Normal) {
					const NormalNbrInfo<1> &info = pinfo.getNormalNbrInfo(e);
					checkIdAndLocalIndex(info.id, info.global_index);
				} else if (type == NbrType::Fine) {
					const FineNbrInfo<1> &info = pinfo.getFineNbrInfo(e);
					for (size_t i = 0; i < info.ids.size(); i++) {
						checkIdAndLocalIndex(info.ids[i], info.global_indexes[i]);
					}
				} else if (type == NbrType::Coarse) {
					const CoarseNbrInfo<1> &info = pinfo.getCoarseNbrInfo(e);
					checkIdAndLocalIndex(info.id, info.global_index);
				}
			}
		}
		for (Corner<3> c : Corner<3>::getValues()) {
			if (pinfo.hasNbr(c)) {
				NbrType type = pinfo.getNbrType(c);
				if (type == NbrType::Normal) {
					const NormalNbrInfo<0> &info = pinfo.getNormalNbrInfo(c);
					checkIdAndLocalIndex(info.id, info.global_index);
				} else if (type == NbrType::Fine) {
					const FineNbrInfo<0> &info = pinfo.getFineNbrInfo(c);
					for (size_t i = 0; i < info.ids.size(); i++) {
						checkIdAndLocalIndex(info.ids[i], info.global_indexes[i]);
					}
				} else if (type == NbrType::Coarse) {
					const CoarseNbrInfo<0> &info = pinfo.getCoarseNbrInfo(c);
					checkIdAndLocalIndex(info.id, info.global_index);
				}
			}
		}
	}
}