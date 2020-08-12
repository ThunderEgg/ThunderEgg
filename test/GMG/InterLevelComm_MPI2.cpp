/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019-2020 ThunderEgg Developers. See AUTHORS.md file at the
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
#include "utils/DomainReader.h"
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/InterLevelComm.h>
#include <ThunderEgg/ValVector.h>
using namespace ThunderEgg;
using namespace std;
const string uniform     = "mesh_inputs/2d_uniform_quad_mpi2.json";
const string mid_uniform = "mesh_inputs/2d_4x4_mid_on_1_mpi2.json";
#define MESHES uniform
#define MESHE_FILES uniform, mid_uniform

TEST_CASE("Check number of local and ghost parents", "[GMG::InterLevelComm]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	auto                  nx        = GENERATE(2);
	auto                  ny        = GENERATE(2);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	INFO("d_fine: " << d_fine->getNumLocalPatches());
	INFO("d_coarse: " << d_coarse->getNumLocalPatches());
	auto ilc = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	INFO("RANK " << rank);
	size_t             num_ghost_parents = 0;
	size_t             num_local_parents = 0;
	map<int, set<int>> my_local_parents_to_children;
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		if (pinfo->parent_rank == rank) {
			num_local_parents++;
		} else {
			num_ghost_parents++;
		}
	}
	CHECK(ilc->getPatchesWithGhostParent().size() == num_ghost_parents);
	CHECK(ilc->getPatchesWithLocalParent().size() == num_local_parents);
}
TEST_CASE("Check that parents have unique local indexes", "[GMG::InterLevelComm]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	auto                  nx        = GENERATE(2);
	auto                  ny        = GENERATE(2);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	INFO("d_fine: " << d_fine->getNumLocalPatches());
	INFO("d_coarse: " << d_coarse->getNumLocalPatches());
	auto ilc = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	INFO("RANK " << rank);

	map<int, set<int>> local_index_id_map;
	for (auto pair : ilc->getPatchesWithLocalParent()) {
		local_index_id_map[pair.first].insert(pair.second->parent_id);
	}
	for (auto pair : local_index_id_map) {
		CHECK(pair.second.size() == 1);
	}
	map<int, set<int>> ghost_local_index_id_map;
	for (auto pair : ilc->getPatchesWithGhostParent()) {
		ghost_local_index_id_map[pair.first].insert(pair.second->parent_id);
	}
	for (auto pair : ghost_local_index_id_map) {
		CHECK(pair.second.size() == 1);
	}
}
TEST_CASE("2-processor getNewGhostVector on uniform quad", "[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto ghost_vec = ilc->getNewGhostVector();

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		CHECK(ghost_vec->getNumLocalPatches() == 0);
	} else {
		CHECK(ghost_vec->getNumLocalPatches() == 1);
	}
}
TEST_CASE("2-processor sendGhostPatches on uniform quad", "[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// info
	INFO("nx: " << nx);
	INFO("ny: " << ny);
	INFO("rank: " << rank);

	auto f = [&](const std::array<double, 2> coord) -> double { return rank + 1; };

	// fill vectors with rank+1
	DomainTools<2>::setValuesWithGhost(d_coarse, coarse_vec, f);
	for (int i = 0; i < ghost_vec->getNumLocalPatches(); i++) {
		auto local_data = ghost_vec->getLocalData(i);
		nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
		               [&](const std::array<int, 2> &coord) { local_data[coord] = rank + 1; });
	}

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec);
	if (rank == 0) {
		// the coarse vec should be filled with 3
		auto local_data = coarse_vec->getLocalData(0);
		nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
		               [&](const std::array<int, 2> &coord) {
			               INFO("xi: " << coord[0]);
			               INFO("yi: " << coord[1]);
			               CHECK(local_data[coord] == 3);
		               });
	} else {
	}
}
TEST_CASE(
"2-processor sendGhostPatches throws exception when start isn't called before finish on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2);
	auto                  ny        = GENERATE(2);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE(
"2-processor sendGhostPatches throws exception when start and finish are called on different ghost vectors on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2);
	auto                  ny        = GENERATE(2);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec   = ilc->getNewGhostVector();
	auto ghost_vec_2 = ilc->getNewGhostVector();

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec_2), RuntimeError);
}
TEST_CASE(
"2-processor sendGhostPatches throws exception when start and finish are called on different vectors on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec   = ValVector<2>::GetNewVector(d_coarse);
	auto coarse_vec_2 = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec_2, ghost_vec), RuntimeError);
}
TEST_CASE(
"2-processor sendGhostPatches throws exception when start and finish are called on different vectors and ghost vectors on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec   = ValVector<2>::GetNewVector(d_coarse);
	auto coarse_vec_2 = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec   = ilc->getNewGhostVector();
	auto ghost_vec_2 = ilc->getNewGhostVector();

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec_2, ghost_vec_2), RuntimeError);
}
TEST_CASE(
"2-processor sendGhostPatches throws exception when start is called twice on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE(
"2-processor sendGhostPatches throws exception when get start is called after send start on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->getGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE(
"2-processor sendGhostPatches throws exception when sned start is called after get start on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE("2-processor getGhostPatches on uniform quad", "[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	auto f = [&](const std::array<double, 2> coord) -> double { return rank + 1; };

	// fill vectors with rank+1
	DomainTools<2>::setValuesWithGhost(d_coarse, coarse_vec, f);
	for (int i = 0; i < ghost_vec->getNumLocalPatches(); i++) {
		auto local_data = ghost_vec->getLocalData(i);
		nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
		               [&](const std::array<int, 2> &coord) { local_data[coord] = rank + 1; });
	}

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	ilc->getGhostPatchesFinish(coarse_vec, ghost_vec);
	if (rank == 0) {
	} else {
		// the coarse vec should be filled with 1
		auto local_data = ghost_vec->getLocalData(0);
		nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
		               [&](const std::array<int, 2> &coord) { CHECK(local_data[coord] == 1); });
	}
}
TEST_CASE(
"2-processor getGhostPatches throws exception when start isn't called before finish on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE(
"2-processor getGhostPatches throws exception when start and finish are called on different ghost vectors on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec   = ilc->getNewGhostVector();
	auto ghost_vec_2 = ilc->getNewGhostVector();

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec, ghost_vec_2), RuntimeError);
}
TEST_CASE(
"2-processor getGhostPatches throws exception when start and finish are called on different vectors on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec   = ValVector<2>::GetNewVector(d_coarse);
	auto coarse_vec_2 = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec_2, ghost_vec), RuntimeError);
}
TEST_CASE(
"2-processor getGhostPatches throws exception when start and finish are called on different vectors and ghost vectors on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec   = ValVector<2>::GetNewVector(d_coarse);
	auto coarse_vec_2 = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec   = ilc->getNewGhostVector();
	auto ghost_vec_2 = ilc->getNewGhostVector();

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec_2, ghost_vec_2), RuntimeError);
}
TEST_CASE("2-processor getGhostPatches throws exception when start is called twice on uniform quad",
          "[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->getGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE("2-processor getGhostPatches throws exception when send finish is called on get start",
          "[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE("2-processor getGhostPatches throws exception when get finish is called on send start",
          "[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE("2-processor getGhostPatches throws exception when send start is called after get start",
          "[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE("2-processor getGhostPatches throws exception when get start is called after send start",
          "[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->getGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE(
"2-processor getGhostPatches throws exception when send finish is called after get start and finish",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	ilc->getGhostPatchesFinish(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE(
"2-processor getGhostPatches throws exception when get finish is called after send start and finish",
"[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE("2-processor getGhostPatches called twice on uniform quad", "[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	auto f = [&](const std::array<double, 2> coord) -> double { return rank + 1; };
	for (int i = 0; i < 2; i++) {
		INFO("Call" << i);
		auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

		auto ghost_vec = ilc->getNewGhostVector();

		// fill vectors with rank+1
		DomainTools<2>::setValuesWithGhost(d_coarse, coarse_vec, f);
		for (int i = 0; i < ghost_vec->getNumLocalPatches(); i++) {
			auto local_data = ghost_vec->getLocalData(i);
			nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
			               [&](const std::array<int, 2> &coord) { local_data[coord] = rank + 1; });
		}

		ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
		ilc->getGhostPatchesFinish(coarse_vec, ghost_vec);
		if (rank == 0) {
		} else {
			// the coarse vec should be filled with 1
			auto local_data = ghost_vec->getLocalData(0);
			nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
			               [&](const std::array<int, 2> &coord) { CHECK(local_data[coord] == 1); });
		}
	}
}
TEST_CASE("2-processor sendGhostPatches called twice on uniform quad", "[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	auto f = [&](const std::array<double, 2> coord) -> double { return rank + 1; };

	for (int i = 0; i < 2; i++) {
		auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);

		auto ghost_vec = ilc->getNewGhostVector();

		// fill vectors with rank+1
		DomainTools<2>::setValuesWithGhost(d_coarse, coarse_vec, f);
		for (int i = 0; i < ghost_vec->getNumLocalPatches(); i++) {
			auto local_data = ghost_vec->getLocalData(i);
			nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
			               [&](const std::array<int, 2> &coord) { local_data[coord] = rank + 1; });
		}

		ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
		ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec);
		if (rank == 0) {
			// the coarse vec should be filled with 3
			auto local_data = coarse_vec->getLocalData(0);
			nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
			               [&](const std::array<int, 2> &coord) { CHECK(local_data[coord] == 3); });
		} else {
		}
	}
}
TEST_CASE("2-processor sendGhostPatches then getGhostPaches called on uniform quad",
          "[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	auto f = [&](const std::array<double, 2> coord) -> double { return rank + 1; };

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);
	auto ghost_vec  = ilc->getNewGhostVector();

	// fill vectors with rank+1
	DomainTools<2>::setValuesWithGhost(d_coarse, coarse_vec, f);
	for (int i = 0; i < ghost_vec->getNumLocalPatches(); i++) {
		auto local_data = ghost_vec->getLocalData(i);
		nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
		               [&](const std::array<int, 2> &coord) { local_data[coord] = rank + 1; });
	}

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec);
	if (rank == 0) {
		// the coarse vec should be filled with 3
		auto local_data = coarse_vec->getLocalData(0);
		nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
		               [&](const std::array<int, 2> &coord) { CHECK(local_data[coord] == 3); });
	} else {
	}

	// fill vectors with rank+1
	DomainTools<2>::setValuesWithGhost(d_coarse, coarse_vec, f);
	for (int i = 0; i < ghost_vec->getNumLocalPatches(); i++) {
		auto local_data = ghost_vec->getLocalData(i);
		nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
		               [&](const std::array<int, 2> &coord) { local_data[coord] = rank + 1; });
	}

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	ilc->getGhostPatchesFinish(coarse_vec, ghost_vec);
	if (rank == 0) {
	} else {
		// the coarse vec should be filled with 1
		auto local_data = ghost_vec->getLocalData(0);
		nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
		               [&](const std::array<int, 2> &coord) { CHECK(local_data[coord] == 1); });
	}
}
TEST_CASE("2-processor getGhostPatches then sendGhostPaches called on uniform quad",
          "[GMG::InterLevelComm]")
{
	auto                  mesh_file = GENERATE(as<std::string>{}, MESHES);
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, d_fine);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	auto f = [&](const std::array<double, 2> coord) -> double { return rank + 1; };

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse);
	auto ghost_vec  = ilc->getNewGhostVector();

	// fill vectors with rank+1
	DomainTools<2>::setValuesWithGhost(d_coarse, coarse_vec, f);
	for (int i = 0; i < ghost_vec->getNumLocalPatches(); i++) {
		auto local_data = ghost_vec->getLocalData(i);
		nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
		               [&](const std::array<int, 2> &coord) { local_data[coord] = rank + 1; });
	}

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	ilc->getGhostPatchesFinish(coarse_vec, ghost_vec);
	if (rank == 0) {
	} else {
		// the coarse vec should be filled with 1
		auto local_data = ghost_vec->getLocalData(0);
		nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
		               [&](const std::array<int, 2> &coord) { CHECK(local_data[coord] == 1); });
	}

	// fill vectors with rank+1
	DomainTools<2>::setValuesWithGhost(d_coarse, coarse_vec, f);
	for (int i = 0; i < ghost_vec->getNumLocalPatches(); i++) {
		auto local_data = ghost_vec->getLocalData(i);
		nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
		               [&](const std::array<int, 2> &coord) { local_data[coord] = rank + 1; });
	}

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec);
	if (rank == 0) {
		// the coarse vec should be filled with 3
		auto local_data = coarse_vec->getLocalData(0);
		nested_loop<2>(local_data.getGhostStart(), local_data.getGhostEnd(),
		               [&](const std::array<int, 2> &coord) { CHECK(local_data[coord] == 3); });
	} else {
	}
}