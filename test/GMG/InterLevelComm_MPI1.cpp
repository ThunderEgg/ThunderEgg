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
#include "../utils/DomainReader.h"
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/InterLevelComm.h>
#include <ThunderEgg/ValVector.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace ThunderEgg;
using namespace std;

const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";

TEST_CASE("1-processor InterLevelComm GetPatches on uniform 4x4", "[GMG::InterLevelComm]")
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

	CHECK(ilc->getPatchesWithGhostParent().size() == 0);
	CHECK(ilc->getPatchesWithLocalParent().size() == 16);

	map<int, set<int>> parents_to_children;
	for (auto pair : ilc->getPatchesWithLocalParent()) {
		parents_to_children[pair.first].insert(pair.second.get().id);
	}
	CHECK(parents_to_children.size() == 4);
	CHECK(parents_to_children.count(0) == 1);
	CHECK(parents_to_children.count(1) == 1);
	CHECK(parents_to_children.count(2) == 1);
	CHECK(parents_to_children.count(3) == 1);

	CHECK(parents_to_children[0].count(13) == 1);
	CHECK(parents_to_children[0].count(14) == 1);
	CHECK(parents_to_children[0].count(15) == 1);
	CHECK(parents_to_children[0].count(16) == 1);

	CHECK(parents_to_children[1].count(17) == 1);
	CHECK(parents_to_children[1].count(18) == 1);
	CHECK(parents_to_children[1].count(19) == 1);
	CHECK(parents_to_children[1].count(20) == 1);

	CHECK(parents_to_children[2].count(5) == 1);
	CHECK(parents_to_children[2].count(6) == 1);
	CHECK(parents_to_children[2].count(7) == 1);
	CHECK(parents_to_children[2].count(8) == 1);

	CHECK(parents_to_children[3].count(9) == 1);
	CHECK(parents_to_children[3].count(10) == 1);
	CHECK(parents_to_children[3].count(11) == 1);
	CHECK(parents_to_children[3].count(12) == 1);
}
TEST_CASE("1-processor getNewGhostVector on uniform 4x4", "[GMG::InterLevelComm]")
{
	auto                  num_components = GENERATE(1, 2, 3);
	auto                  nx             = GENERATE(2, 10);
	auto                  ny             = GENERATE(2, 10);
	int                   num_ghost      = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, num_components, d_fine);

	auto ghost_vec = ilc->getNewGhostVector();

	CHECK(ghost_vec->getNumLocalPatches() == 0);
	CHECK(ghost_vec->getNumComponents() == num_components);
}
TEST_CASE("1-processor sendGhostPatches on uniform 4x4", "[GMG::InterLevelComm]")
{
	auto                  num_components = GENERATE(1, 2, 3);
	auto                  nx             = GENERATE(2, 10);
	auto                  ny             = GENERATE(2, 10);
	int                   num_ghost      = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine   = domain_reader.getFinerDomain();
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();
	auto                  ilc      = std::make_shared<GMG::InterLevelComm<2>>(d_coarse, num_components, d_fine);

	auto coarse_vec      = ValVector<2>::GetNewVector(d_coarse, num_components);
	auto coarse_expected = ValVector<2>::GetNewVector(d_coarse, num_components);

	auto ghost_vec = ilc->getNewGhostVector();

	// info
	INFO("nx: " << nx);
	INFO("ny: " << ny);

	auto f
	= [&](const std::array<double, 2> coord) -> double { return 1 + coord[0] + 2 * coord[1]; };

	DomainTools::SetValuesWithGhost<2>(d_coarse, coarse_vec, f);
	int idx = 0;
	for (int i = 0; i < coarse_vec->getNumLocalPatches(); i++) {
		auto vec_lds = coarse_vec->getViews(i);
		nested_loop<2>(vec_lds[0].getStart(), vec_lds[0].getEnd(), [&](const array<int, 2> &coord) {
			for (size_t c = 0; c < vec_lds.size(); c++) {
				vec_lds[c][coord] = idx;
				idx++;
			}
		});
	}
	// since there are not ghost patches, the coarse vec should not be modified
	coarse_expected->copy(coarse_vec);

	ilc->sendGhostPatchesStart(coarse_vec, ghost_vec);
	ilc->sendGhostPatchesFinish(coarse_vec, ghost_vec);
	for (int i = 0; i < coarse_vec->getNumLocalPatches(); i++) {
		auto vec_lds      = coarse_vec->getViews(i);
		auto expected_lds = coarse_expected->getViews(i);
		nested_loop<2>(vec_lds[0].getStart(), vec_lds[0].getEnd(), [&](const array<int, 2> &coord) {
			for (size_t c = 0; c < vec_lds.size(); c++) {
				REQUIRE(vec_lds[c][coord] == Catch::Approx(expected_lds[c][coord]));
			}
		});
	}
}