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
const string mesh_file = "mesh_inputs/2d_uniform_quad_mpi3.json";
TEST_CASE("3-processor InterLevelComm GetPatches on uniform quad", "[GMG::InterLevelComm]")
{
	auto                  nx        = GENERATE(2);
	auto                  ny        = GENERATE(2);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	INFO("d_fine: " << d_fine->getNumLocalPatches());
	INFO("d_coarse: " << d_coarse->getNumLocalPatches());
	auto ilc = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, 1, d_fine);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		CHECK(ilc->getPatchesWithGhostParent().size() == 0);
		CHECK(ilc->getPatchesWithLocalParent().size() == 2);

		map<int, set<int>> parents_to_children;
		for (auto pair : ilc->getPatchesWithLocalParent()) {
			parents_to_children[pair.first].insert(pair.second->id);
		}
		CHECK(parents_to_children.count(0));
		CHECK(parents_to_children[0].count(1));
		CHECK(parents_to_children[0].count(2));
	} else if (rank == 1) {
		CHECK(ilc->getPatchesWithGhostParent().size() == 1);
		CHECK(ilc->getPatchesWithLocalParent().size() == 0);

		map<int, set<int>> parents_to_children;
		for (auto pair : ilc->getPatchesWithGhostParent()) {
			parents_to_children[pair.first].insert(pair.second->id);
		}
		CHECK(parents_to_children.count(0) == 1);
		CHECK(parents_to_children[0].count(3) == 1);
	} else {
		CHECK(ilc->getPatchesWithGhostParent().size() == 1);
		CHECK(ilc->getPatchesWithLocalParent().size() == 0);

		map<int, set<int>> parents_to_children;
		for (auto pair : ilc->getPatchesWithGhostParent()) {
			parents_to_children[pair.first].insert(pair.second->id);
		}
		CHECK(parents_to_children.count(0) == 1);
		CHECK(parents_to_children[0].count(4) == 1);
	}
}
TEST_CASE("3-processor getNewGhostVector on uniform quad", "[GMG::InterLevelComm]")
{
	auto                  num_components = GENERATE(1, 2, 3);
	auto                  nx             = GENERATE(2, 10);
	auto                  ny             = GENERATE(2, 10);
	int                   num_ghost      = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto ilc = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, num_components, d_fine);

	auto ghost_vec = ilc->getNewGhostVector();

	CHECK(ghost_vec->getNumComponents() == num_components);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		CHECK(ghost_vec->getNumLocalPatches() == 0);
	} else {
		CHECK(ghost_vec->getNumLocalPatches() == 1);
	}
}
TEST_CASE("3-processor sendGhostPatches on uniform quad", "[GMG::InterLevelComm]")
{
	auto                  num_components = GENERATE(1, 2, 3);
	auto                  nx             = GENERATE(2, 10);
	auto                  ny             = GENERATE(2, 10);
	int                   num_ghost      = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto ilc = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, num_components, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse, num_components);

	auto ghost_vec = ilc->getNewGhostVector();

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// fill vectors with rank+c+1
	for (int i = 0; i < coarse_vec->getNumLocalPatches(); i++) {
		auto local_datas = coarse_vec->getLocalDatas(i);
		nested_loop<2>(local_datas[0].getGhostStart(), local_datas[0].getGhostEnd(),
		               [&](const std::array<int, 2> &coord) {
			               for (int c = 0; c < num_components; c++) {
				               local_datas[c][coord] = rank + c + 1;
			               }
		               });
	}
	for (int i = 0; i < ghost_vec->getNumLocalPatches(); i++) {
		auto local_datas = ghost_vec->getLocalDatas(i);
		nested_loop<2>(local_datas[0].getGhostStart(), local_datas[0].getGhostEnd(),
		               [&](const std::array<int, 2> &coord) {
			               for (int c = 0; c < num_components; c++) {
				               local_datas[c][coord] = rank + c + 1;
			               }
		               });
	}

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec);
	if (rank == 0) {
		// the coarse vec should be filled with 6+2*c
		auto local_datas = coarse_vec->getLocalDatas(0);
		nested_loop<2>(local_datas[0].getGhostStart(), local_datas[0].getGhostEnd(),
		               [&](const std::array<int, 2> &coord) {
			               for (int c = 0; c < num_components; c++) {
				               INFO("c " << c);
				               INFO("xi: " << coord[0]);
				               INFO("yi: " << coord[1]);
				               CHECK(local_datas[c][coord] == 1 + 2 + 3 + 3 * c);
			               }
		               });
	} else {
	}
}
TEST_CASE(
"3-processor sendGhostPatches throws exception when start isn't called before finish on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, 1, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse, 1);

	auto ghost_vec = ilc->getNewGhostVector();

	CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE(
"3-processor sendGhostPatches throws exception when start and finish are called on different ghost vectors on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, 1, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse, 1);

	auto ghost_vec   = ilc->getNewGhostVector();
	auto ghost_vec_2 = ilc->getNewGhostVector();

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec_2), RuntimeError);
}
TEST_CASE(
"3-processor sendGhostPatches throws exception when start and finish are called on different vectors on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, 1, d_fine);

	auto coarse_vec   = ValVector<2>::GetNewVector(d_coarse, 1);
	auto coarse_vec_2 = ValVector<2>::GetNewVector(d_coarse, 1);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec_2, ghost_vec), RuntimeError);
}
TEST_CASE(
"3-processor sendGhostPatches throws exception when start and finish are called on different vectors and ghost vectors on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, 1, d_fine);

	auto coarse_vec   = ValVector<2>::GetNewVector(d_coarse, 1);
	auto coarse_vec_2 = ValVector<2>::GetNewVector(d_coarse, 1);

	auto ghost_vec   = ilc->getNewGhostVector();
	auto ghost_vec_2 = ilc->getNewGhostVector();

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec_2, ghost_vec_2), RuntimeError);
}
TEST_CASE(
"3-processor sendGhostPatches throws exception when start is called twice on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, 1, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse, 1);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE("3-processor getGhostPatches on uniform quad", "[GMG::InterLevelComm]")
{
	auto                  num_components = GENERATE(1, 2, 3);
	auto                  nx             = GENERATE(2, 10);
	auto                  ny             = GENERATE(2, 10);
	int                   num_ghost      = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto ilc = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, num_components, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse, num_components);

	auto ghost_vec = ilc->getNewGhostVector();

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// fill vectors with rank+c+1
	for (int i = 0; i < coarse_vec->getNumLocalPatches(); i++) {
		auto local_datas = coarse_vec->getLocalDatas(i);
		nested_loop<2>(local_datas[0].getGhostStart(), local_datas[0].getGhostEnd(),
		               [&](const std::array<int, 2> &coord) {
			               for (int c = 0; c < num_components; c++) {
				               local_datas[c][coord] = rank + c + 1;
			               }
		               });
	}
	for (int i = 0; i < ghost_vec->getNumLocalPatches(); i++) {
		auto local_datas = ghost_vec->getLocalDatas(i);
		nested_loop<2>(local_datas[0].getGhostStart(), local_datas[0].getGhostEnd(),
		               [&](const std::array<int, 2> &coord) {
			               for (int c = 0; c < num_components; c++) {
				               local_datas[c][coord] = rank + c + 1;
			               }
		               });
	}

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	ilc->getGhostPatchesFinish(coarse_vec, ghost_vec);
	if (rank == 0) {
	} else {
		// the coarse vec should be filled with 1+c
		auto local_datas = ghost_vec->getLocalDatas(0);
		nested_loop<2>(local_datas[0].getGhostStart(), local_datas[0].getGhostEnd(),
		               [&](const std::array<int, 2> &coord) {
			               for (int c = 0; c < num_components; c++) {
				               CHECK(local_datas[c][coord] == 1 + c);
			               }
		               });
	}
}
TEST_CASE(
"3-processor getGhostPatches throws exception when start isn't called before finish on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, 1, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse, 1);

	auto ghost_vec = ilc->getNewGhostVector();

	CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec, ghost_vec), RuntimeError);
}
TEST_CASE(
"3-processor getGhostPatches throws exception when start and finish are called on different ghost vectors on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, 1, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse, 1);

	auto ghost_vec   = ilc->getNewGhostVector();
	auto ghost_vec_2 = ilc->getNewGhostVector();

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec, ghost_vec_2), RuntimeError);
}
TEST_CASE(
"3-processor getGhostPatches throws exception when start and finish are called on different vectors on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, 1, d_fine);

	auto coarse_vec   = ValVector<2>::GetNewVector(d_coarse, 1);
	auto coarse_vec_2 = ValVector<2>::GetNewVector(d_coarse, 1);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->sendGhostPatchesFinish(coarse_vec_2, ghost_vec), RuntimeError);
}
TEST_CASE(
"3-processor getGhostPatches throws exception when start and finish are called on different vectors and ghost vectors on uniform quad",
"[GMG::InterLevelComm]")
{
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, 1, d_fine);

	auto coarse_vec   = ValVector<2>::GetNewVector(d_coarse, 1);
	auto coarse_vec_2 = ValVector<2>::GetNewVector(d_coarse, 1);

	auto ghost_vec   = ilc->getNewGhostVector();
	auto ghost_vec_2 = ilc->getNewGhostVector();

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->getGhostPatchesFinish(coarse_vec_2, ghost_vec_2), RuntimeError);
}
TEST_CASE("3-processor getGhostPatches throws exception when start is called twice on uniform quad",
          "[GMG::InterLevelComm]")
{
	auto                  nx        = GENERATE(2, 10);
	auto                  ny        = GENERATE(2, 10);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, 1, d_fine);

	auto coarse_vec = ValVector<2>::GetNewVector(d_coarse, 1);

	auto ghost_vec = ilc->getNewGhostVector();

	ilc->getGhostPatchesStart(coarse_vec, ghost_vec);
	CHECK_THROWS_AS(ilc->getGhostPatchesStart(coarse_vec, ghost_vec), RuntimeError);
}