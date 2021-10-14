/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2021      Scott Aiton
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

#include "P8estDomainGenerator_SHARED.h"
#include <catch2/catch_approx.hpp>
using namespace std;
using namespace ThunderEgg;

#include <catch2/catch_test_macros.hpp>

void Check2x2x2DomainNeighbors(const ThunderEgg::Domain<3> &domain)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<PatchInfo<3>> patches = GetAllPatchesOnRank0(domain);
	if (rank == 0) {
		const PatchInfo<3> *domain_bsw_patch = nullptr;
		const PatchInfo<3> *domain_bse_patch = nullptr;
		const PatchInfo<3> *domain_bnw_patch = nullptr;
		const PatchInfo<3> *domain_bne_patch = nullptr;
		const PatchInfo<3> *domain_tsw_patch = nullptr;
		const PatchInfo<3> *domain_tse_patch = nullptr;
		const PatchInfo<3> *domain_tnw_patch = nullptr;
		const PatchInfo<3> *domain_tne_patch = nullptr;

		for (PatchInfo<3> &patch : patches) {
			double x = patch.starts[0];
			double y = patch.starts[1];
			double z = patch.starts[2];
			if (x == Catch::Approx(0) && y == Catch::Approx(0) && z == Catch::Approx(0)) {
				domain_bsw_patch = &patch;
			}
			if (x == Catch::Approx(0.5) && y == Catch::Approx(0) && z == Catch::Approx(0)) {
				domain_bse_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.5) && z == Catch::Approx(0)) {
				domain_bnw_patch = &patch;
			}
			if (x == Catch::Approx(0.5) && y == Catch::Approx(0.5) && z == Catch::Approx(0)) {
				domain_bne_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0) && z == Catch::Approx(0.5)) {
				domain_tsw_patch = &patch;
			}
			if (x == Catch::Approx(0.5) && y == Catch::Approx(0) && z == Catch::Approx(0.5)) {
				domain_tse_patch = &patch;
			}
			if (x == Catch::Approx(0) && y == Catch::Approx(0.5) && z == Catch::Approx(0.5)) {
				domain_tnw_patch = &patch;
			}
			if (x == Catch::Approx(0.5) && y == Catch::Approx(0.5) && z == Catch::Approx(0.5)) {
				domain_tne_patch = &patch;
			}
		}

		REQUIRE(domain_bsw_patch != nullptr);
		REQUIRE(domain_bse_patch != nullptr);
		REQUIRE(domain_bnw_patch != nullptr);
		REQUIRE(domain_bne_patch != nullptr);
		REQUIRE(domain_tsw_patch != nullptr);
		REQUIRE(domain_tse_patch != nullptr);
		REQUIRE(domain_tnw_patch != nullptr);
		REQUIRE(domain_tne_patch != nullptr);

		//Side hasnbr
		CHECK(domain_bsw_patch->hasNbr(Side<3>::west()) == false);
		CHECK(domain_bsw_patch->hasNbr(Side<3>::east()) == true);
		CHECK(domain_bsw_patch->hasNbr(Side<3>::south()) == false);
		CHECK(domain_bsw_patch->hasNbr(Side<3>::north()) == true);
		CHECK(domain_bsw_patch->hasNbr(Side<3>::bottom()) == false);
		CHECK(domain_bsw_patch->hasNbr(Side<3>::top()) == true);

		CHECK(domain_bse_patch->hasNbr(Side<3>::west()) == true);
		CHECK(domain_bse_patch->hasNbr(Side<3>::east()) == false);
		CHECK(domain_bse_patch->hasNbr(Side<3>::south()) == false);
		CHECK(domain_bse_patch->hasNbr(Side<3>::north()) == true);
		CHECK(domain_bse_patch->hasNbr(Side<3>::bottom()) == false);
		CHECK(domain_bse_patch->hasNbr(Side<3>::top()) == true);

		CHECK(domain_bnw_patch->hasNbr(Side<3>::west()) == false);
		CHECK(domain_bnw_patch->hasNbr(Side<3>::east()) == true);
		CHECK(domain_bnw_patch->hasNbr(Side<3>::south()) == true);
		CHECK(domain_bnw_patch->hasNbr(Side<3>::north()) == false);
		CHECK(domain_bnw_patch->hasNbr(Side<3>::bottom()) == false);
		CHECK(domain_bnw_patch->hasNbr(Side<3>::top()) == true);

		CHECK(domain_bne_patch->hasNbr(Side<3>::west()) == true);
		CHECK(domain_bne_patch->hasNbr(Side<3>::east()) == false);
		CHECK(domain_bne_patch->hasNbr(Side<3>::south()) == true);
		CHECK(domain_bne_patch->hasNbr(Side<3>::north()) == false);
		CHECK(domain_bne_patch->hasNbr(Side<3>::bottom()) == false);
		CHECK(domain_bne_patch->hasNbr(Side<3>::top()) == true);

		CHECK(domain_tsw_patch->hasNbr(Side<3>::west()) == false);
		CHECK(domain_tsw_patch->hasNbr(Side<3>::east()) == true);
		CHECK(domain_tsw_patch->hasNbr(Side<3>::south()) == false);
		CHECK(domain_tsw_patch->hasNbr(Side<3>::north()) == true);
		CHECK(domain_tsw_patch->hasNbr(Side<3>::bottom()) == true);
		CHECK(domain_tsw_patch->hasNbr(Side<3>::top()) == false);

		CHECK(domain_tse_patch->hasNbr(Side<3>::west()) == true);
		CHECK(domain_tse_patch->hasNbr(Side<3>::east()) == false);
		CHECK(domain_tse_patch->hasNbr(Side<3>::south()) == false);
		CHECK(domain_tse_patch->hasNbr(Side<3>::north()) == true);
		CHECK(domain_tse_patch->hasNbr(Side<3>::bottom()) == true);
		CHECK(domain_tse_patch->hasNbr(Side<3>::top()) == false);

		CHECK(domain_tnw_patch->hasNbr(Side<3>::west()) == false);
		CHECK(domain_tnw_patch->hasNbr(Side<3>::east()) == true);
		CHECK(domain_tnw_patch->hasNbr(Side<3>::south()) == true);
		CHECK(domain_tnw_patch->hasNbr(Side<3>::north()) == false);
		CHECK(domain_tnw_patch->hasNbr(Side<3>::bottom()) == true);
		CHECK(domain_tnw_patch->hasNbr(Side<3>::top()) == false);

		CHECK(domain_tne_patch->hasNbr(Side<3>::west()) == true);
		CHECK(domain_tne_patch->hasNbr(Side<3>::east()) == false);
		CHECK(domain_tne_patch->hasNbr(Side<3>::south()) == true);
		CHECK(domain_tne_patch->hasNbr(Side<3>::north()) == false);
		CHECK(domain_tne_patch->hasNbr(Side<3>::bottom()) == true);
		CHECK(domain_tne_patch->hasNbr(Side<3>::top()) == false);

		//SIDE ids
		CHECK(domain_bsw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_bse_patch->id);
		CHECK(domain_bsw_patch->getNormalNbrInfo(Side<3>::north()).id == domain_bnw_patch->id);
		CHECK(domain_bsw_patch->getNormalNbrInfo(Side<3>::top()).id == domain_tsw_patch->id);

		CHECK(domain_bse_patch->getNormalNbrInfo(Side<3>::west()).id == domain_bsw_patch->id);
		CHECK(domain_bse_patch->getNormalNbrInfo(Side<3>::north()).id == domain_bne_patch->id);
		CHECK(domain_bse_patch->getNormalNbrInfo(Side<3>::top()).id == domain_tse_patch->id);

		CHECK(domain_bnw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_bne_patch->id);
		CHECK(domain_bnw_patch->getNormalNbrInfo(Side<3>::south()).id == domain_bsw_patch->id);
		CHECK(domain_bnw_patch->getNormalNbrInfo(Side<3>::top()).id == domain_tnw_patch->id);

		CHECK(domain_bne_patch->getNormalNbrInfo(Side<3>::west()).id == domain_bnw_patch->id);
		CHECK(domain_bne_patch->getNormalNbrInfo(Side<3>::south()).id == domain_bse_patch->id);
		CHECK(domain_bne_patch->getNormalNbrInfo(Side<3>::top()).id == domain_tne_patch->id);

		CHECK(domain_tsw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_tse_patch->id);
		CHECK(domain_tsw_patch->getNormalNbrInfo(Side<3>::north()).id == domain_tnw_patch->id);
		CHECK(domain_tsw_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_bsw_patch->id);

		CHECK(domain_tse_patch->getNormalNbrInfo(Side<3>::west()).id == domain_tsw_patch->id);
		CHECK(domain_tse_patch->getNormalNbrInfo(Side<3>::north()).id == domain_tne_patch->id);
		CHECK(domain_tse_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_bse_patch->id);

		CHECK(domain_tnw_patch->getNormalNbrInfo(Side<3>::east()).id == domain_tne_patch->id);
		CHECK(domain_tnw_patch->getNormalNbrInfo(Side<3>::south()).id == domain_tsw_patch->id);
		CHECK(domain_tnw_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_bnw_patch->id);

		CHECK(domain_tne_patch->getNormalNbrInfo(Side<3>::west()).id == domain_tnw_patch->id);
		CHECK(domain_tne_patch->getNormalNbrInfo(Side<3>::south()).id == domain_tse_patch->id);
		CHECK(domain_tne_patch->getNormalNbrInfo(Side<3>::bottom()).id == domain_bne_patch->id);

		//SIDE ranks
		CHECK(domain_bsw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_bse_patch->rank);
		CHECK(domain_bsw_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_bnw_patch->rank);
		CHECK(domain_bsw_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_tsw_patch->rank);

		CHECK(domain_bse_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_bsw_patch->rank);
		CHECK(domain_bse_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_bne_patch->rank);
		CHECK(domain_bse_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_tse_patch->rank);

		CHECK(domain_bnw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_bne_patch->rank);
		CHECK(domain_bnw_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_bsw_patch->rank);
		CHECK(domain_bnw_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_tnw_patch->rank);

		CHECK(domain_bne_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_bnw_patch->rank);
		CHECK(domain_bne_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_bse_patch->rank);
		CHECK(domain_bne_patch->getNormalNbrInfo(Side<3>::top()).rank == domain_tne_patch->rank);

		CHECK(domain_tsw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_tse_patch->rank);
		CHECK(domain_tsw_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_tnw_patch->rank);
		CHECK(domain_tsw_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_bsw_patch->rank);

		CHECK(domain_tse_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_tsw_patch->rank);
		CHECK(domain_tse_patch->getNormalNbrInfo(Side<3>::north()).rank == domain_tne_patch->rank);
		CHECK(domain_tse_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_bse_patch->rank);

		CHECK(domain_tnw_patch->getNormalNbrInfo(Side<3>::east()).rank == domain_tne_patch->rank);
		CHECK(domain_tnw_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_tsw_patch->rank);
		CHECK(domain_tnw_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_bnw_patch->rank);

		CHECK(domain_tne_patch->getNormalNbrInfo(Side<3>::west()).rank == domain_tnw_patch->rank);
		CHECK(domain_tne_patch->getNormalNbrInfo(Side<3>::south()).rank == domain_tse_patch->rank);
		CHECK(domain_tne_patch->getNormalNbrInfo(Side<3>::bottom()).rank == domain_bne_patch->rank);

		//Edge hasnbr
		CHECK_FALSE(domain_bsw_patch->hasNbr(Edge::bs()));
		CHECK(domain_bsw_patch->hasNbr(Edge::tn()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Edge::bn()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Edge::ts()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Edge::bw()));
		CHECK(domain_bsw_patch->hasNbr(Edge::te()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Edge::be()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Edge::tw()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Edge::sw()));
		CHECK(domain_bsw_patch->hasNbr(Edge::ne()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Edge::se()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Edge::nw()));

		CHECK_FALSE(domain_bse_patch->hasNbr(Edge::bs()));
		CHECK(domain_bse_patch->hasNbr(Edge::tn()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Edge::bn()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Edge::ts()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Edge::bw()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Edge::te()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Edge::be()));
		CHECK(domain_bse_patch->hasNbr(Edge::tw()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Edge::sw()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Edge::ne()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Edge::se()));
		CHECK(domain_bse_patch->hasNbr(Edge::nw()));

		CHECK_FALSE(domain_bnw_patch->hasNbr(Edge::bs()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Edge::tn()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Edge::bn()));
		CHECK(domain_bnw_patch->hasNbr(Edge::ts()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Edge::bw()));
		CHECK(domain_bnw_patch->hasNbr(Edge::te()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Edge::be()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Edge::tw()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Edge::sw()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Edge::ne()));
		CHECK(domain_bnw_patch->hasNbr(Edge::se()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Edge::nw()));

		CHECK_FALSE(domain_bne_patch->hasNbr(Edge::bs()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Edge::tn()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Edge::bn()));
		CHECK(domain_bne_patch->hasNbr(Edge::ts()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Edge::bw()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Edge::te()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Edge::be()));
		CHECK(domain_bne_patch->hasNbr(Edge::tw()));
		CHECK(domain_bne_patch->hasNbr(Edge::sw()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Edge::ne()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Edge::se()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Edge::nw()));

		CHECK_FALSE(domain_tsw_patch->hasNbr(Edge::bs()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Edge::tn()));
		CHECK(domain_tsw_patch->hasNbr(Edge::bn()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Edge::ts()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Edge::bw()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Edge::te()));
		CHECK(domain_tsw_patch->hasNbr(Edge::be()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Edge::tw()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Edge::sw()));
		CHECK(domain_tsw_patch->hasNbr(Edge::ne()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Edge::se()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Edge::nw()));

		CHECK_FALSE(domain_tse_patch->hasNbr(Edge::bs()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Edge::tn()));
		CHECK(domain_tse_patch->hasNbr(Edge::bn()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Edge::ts()));
		CHECK(domain_tse_patch->hasNbr(Edge::bw()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Edge::te()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Edge::be()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Edge::tw()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Edge::sw()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Edge::ne()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Edge::se()));
		CHECK(domain_tse_patch->hasNbr(Edge::nw()));

		CHECK(domain_tnw_patch->hasNbr(Edge::bs()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Edge::tn()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Edge::bn()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Edge::ts()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Edge::bw()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Edge::te()));
		CHECK(domain_tnw_patch->hasNbr(Edge::be()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Edge::tw()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Edge::sw()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Edge::ne()));
		CHECK(domain_tnw_patch->hasNbr(Edge::se()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Edge::nw()));

		CHECK(domain_tne_patch->hasNbr(Edge::bs()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Edge::tn()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Edge::bn()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Edge::ts()));
		CHECK(domain_tne_patch->hasNbr(Edge::bw()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Edge::te()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Edge::be()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Edge::tw()));
		CHECK(domain_tne_patch->hasNbr(Edge::sw()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Edge::ne()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Edge::se()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Edge::nw()));

		//Edge nbr id
		CHECK(domain_bsw_patch->getNormalNbrInfo(Edge::tn()).id == domain_tnw_patch->id);
		CHECK(domain_bsw_patch->getNormalNbrInfo(Edge::te()).id == domain_tse_patch->id);
		CHECK(domain_bsw_patch->getNormalNbrInfo(Edge::ne()).id == domain_bne_patch->id);

		CHECK(domain_bse_patch->getNormalNbrInfo(Edge::tn()).id == domain_tne_patch->id);
		CHECK(domain_bse_patch->getNormalNbrInfo(Edge::tw()).id == domain_tsw_patch->id);
		CHECK(domain_bse_patch->getNormalNbrInfo(Edge::nw()).id == domain_bnw_patch->id);

		CHECK(domain_bnw_patch->getNormalNbrInfo(Edge::ts()).id == domain_tsw_patch->id);
		CHECK(domain_bnw_patch->getNormalNbrInfo(Edge::te()).id == domain_tne_patch->id);
		CHECK(domain_bnw_patch->getNormalNbrInfo(Edge::se()).id == domain_bse_patch->id);

		CHECK(domain_bne_patch->getNormalNbrInfo(Edge::ts()).id == domain_tse_patch->id);
		CHECK(domain_bne_patch->getNormalNbrInfo(Edge::tw()).id == domain_tnw_patch->id);
		CHECK(domain_bne_patch->getNormalNbrInfo(Edge::sw()).id == domain_bsw_patch->id);

		CHECK(domain_tsw_patch->getNormalNbrInfo(Edge::bn()).id == domain_bnw_patch->id);
		CHECK(domain_tsw_patch->getNormalNbrInfo(Edge::be()).id == domain_bse_patch->id);
		CHECK(domain_tsw_patch->getNormalNbrInfo(Edge::ne()).id == domain_tne_patch->id);

		CHECK(domain_tse_patch->getNormalNbrInfo(Edge::bn()).id == domain_bne_patch->id);
		CHECK(domain_tse_patch->getNormalNbrInfo(Edge::bw()).id == domain_bsw_patch->id);
		CHECK(domain_tse_patch->getNormalNbrInfo(Edge::nw()).id == domain_tnw_patch->id);

		CHECK(domain_tnw_patch->getNormalNbrInfo(Edge::bs()).id == domain_bsw_patch->id);
		CHECK(domain_tnw_patch->getNormalNbrInfo(Edge::be()).id == domain_bne_patch->id);
		CHECK(domain_tnw_patch->getNormalNbrInfo(Edge::se()).id == domain_tse_patch->id);

		CHECK(domain_tne_patch->getNormalNbrInfo(Edge::bs()).id == domain_bse_patch->id);
		CHECK(domain_tne_patch->getNormalNbrInfo(Edge::bw()).id == domain_bnw_patch->id);
		CHECK(domain_tne_patch->getNormalNbrInfo(Edge::sw()).id == domain_tsw_patch->id);

		//Edge nbr rank
		CHECK(domain_bsw_patch->getNormalNbrInfo(Edge::tn()).rank == domain_tnw_patch->rank);
		CHECK(domain_bsw_patch->getNormalNbrInfo(Edge::te()).rank == domain_tse_patch->rank);
		CHECK(domain_bsw_patch->getNormalNbrInfo(Edge::ne()).rank == domain_bne_patch->rank);

		CHECK(domain_bse_patch->getNormalNbrInfo(Edge::tn()).rank == domain_tne_patch->rank);
		CHECK(domain_bse_patch->getNormalNbrInfo(Edge::tw()).rank == domain_tsw_patch->rank);
		CHECK(domain_bse_patch->getNormalNbrInfo(Edge::nw()).rank == domain_bnw_patch->rank);

		CHECK(domain_bnw_patch->getNormalNbrInfo(Edge::ts()).rank == domain_tsw_patch->rank);
		CHECK(domain_bnw_patch->getNormalNbrInfo(Edge::te()).rank == domain_tne_patch->rank);
		CHECK(domain_bnw_patch->getNormalNbrInfo(Edge::se()).rank == domain_bse_patch->rank);

		CHECK(domain_bne_patch->getNormalNbrInfo(Edge::ts()).rank == domain_tse_patch->rank);
		CHECK(domain_bne_patch->getNormalNbrInfo(Edge::tw()).rank == domain_tnw_patch->rank);
		CHECK(domain_bne_patch->getNormalNbrInfo(Edge::sw()).rank == domain_bsw_patch->rank);

		CHECK(domain_tsw_patch->getNormalNbrInfo(Edge::bn()).rank == domain_bnw_patch->rank);
		CHECK(domain_tsw_patch->getNormalNbrInfo(Edge::be()).rank == domain_bse_patch->rank);
		CHECK(domain_tsw_patch->getNormalNbrInfo(Edge::ne()).rank == domain_tne_patch->rank);

		CHECK(domain_tse_patch->getNormalNbrInfo(Edge::bn()).rank == domain_bne_patch->rank);
		CHECK(domain_tse_patch->getNormalNbrInfo(Edge::bw()).rank == domain_bsw_patch->rank);
		CHECK(domain_tse_patch->getNormalNbrInfo(Edge::nw()).rank == domain_tnw_patch->rank);

		CHECK(domain_tnw_patch->getNormalNbrInfo(Edge::bs()).rank == domain_bsw_patch->rank);
		CHECK(domain_tnw_patch->getNormalNbrInfo(Edge::be()).rank == domain_bne_patch->rank);
		CHECK(domain_tnw_patch->getNormalNbrInfo(Edge::se()).rank == domain_tse_patch->rank);

		CHECK(domain_tne_patch->getNormalNbrInfo(Edge::bs()).rank == domain_bse_patch->rank);
		CHECK(domain_tne_patch->getNormalNbrInfo(Edge::bw()).rank == domain_bnw_patch->rank);
		CHECK(domain_tne_patch->getNormalNbrInfo(Edge::sw()).rank == domain_tsw_patch->rank);

		// Corner hasNbr
		CHECK_FALSE(domain_bsw_patch->hasNbr(Corner<3>::bsw()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Corner<3>::bse()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Corner<3>::bnw()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Corner<3>::bne()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Corner<3>::tsw()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Corner<3>::tse()));
		CHECK_FALSE(domain_bsw_patch->hasNbr(Corner<3>::tnw()));
		CHECK(domain_bsw_patch->hasNbr(Corner<3>::tne()));

		CHECK_FALSE(domain_bse_patch->hasNbr(Corner<3>::bsw()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Corner<3>::bse()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Corner<3>::bnw()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Corner<3>::bne()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Corner<3>::tsw()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Corner<3>::tse()));
		CHECK(domain_bse_patch->hasNbr(Corner<3>::tnw()));
		CHECK_FALSE(domain_bse_patch->hasNbr(Corner<3>::tne()));

		CHECK_FALSE(domain_bnw_patch->hasNbr(Corner<3>::bsw()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Corner<3>::bse()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Corner<3>::bnw()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Corner<3>::bne()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Corner<3>::tsw()));
		CHECK(domain_bnw_patch->hasNbr(Corner<3>::tse()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Corner<3>::tnw()));
		CHECK_FALSE(domain_bnw_patch->hasNbr(Corner<3>::tne()));

		CHECK_FALSE(domain_bne_patch->hasNbr(Corner<3>::bsw()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Corner<3>::bse()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Corner<3>::bnw()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Corner<3>::bne()));
		CHECK(domain_bne_patch->hasNbr(Corner<3>::tsw()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Corner<3>::tse()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Corner<3>::tnw()));
		CHECK_FALSE(domain_bne_patch->hasNbr(Corner<3>::tne()));

		CHECK_FALSE(domain_tsw_patch->hasNbr(Corner<3>::bsw()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Corner<3>::bse()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Corner<3>::bnw()));
		CHECK(domain_tsw_patch->hasNbr(Corner<3>::bne()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Corner<3>::tsw()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Corner<3>::tse()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Corner<3>::tnw()));
		CHECK_FALSE(domain_tsw_patch->hasNbr(Corner<3>::tne()));

		CHECK_FALSE(domain_tse_patch->hasNbr(Corner<3>::bsw()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Corner<3>::bse()));
		CHECK(domain_tse_patch->hasNbr(Corner<3>::bnw()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Corner<3>::bne()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Corner<3>::tsw()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Corner<3>::tse()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Corner<3>::tnw()));
		CHECK_FALSE(domain_tse_patch->hasNbr(Corner<3>::tne()));

		CHECK_FALSE(domain_tnw_patch->hasNbr(Corner<3>::bsw()));
		CHECK(domain_tnw_patch->hasNbr(Corner<3>::bse()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Corner<3>::bnw()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Corner<3>::bne()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Corner<3>::tsw()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Corner<3>::tse()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Corner<3>::tnw()));
		CHECK_FALSE(domain_tnw_patch->hasNbr(Corner<3>::tne()));

		CHECK(domain_tne_patch->hasNbr(Corner<3>::bsw()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Corner<3>::bse()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Corner<3>::bnw()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Corner<3>::bne()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Corner<3>::tsw()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Corner<3>::tse()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Corner<3>::tnw()));
		CHECK_FALSE(domain_tne_patch->hasNbr(Corner<3>::tne()));

		// Corner nbr id
		CHECK(domain_bsw_patch->getNormalNbrInfo(Corner<3>::tne()).id == domain_tne_patch->id);

		CHECK(domain_bse_patch->getNormalNbrInfo(Corner<3>::tnw()).id == domain_tnw_patch->id);

		CHECK(domain_bnw_patch->getNormalNbrInfo(Corner<3>::tse()).id == domain_tse_patch->id);

		CHECK(domain_bne_patch->getNormalNbrInfo(Corner<3>::tsw()).id == domain_tsw_patch->id);

		CHECK(domain_tsw_patch->getNormalNbrInfo(Corner<3>::bne()).id == domain_bne_patch->id);

		CHECK(domain_tse_patch->getNormalNbrInfo(Corner<3>::bnw()).id == domain_bnw_patch->id);

		CHECK(domain_tnw_patch->getNormalNbrInfo(Corner<3>::bse()).id == domain_bse_patch->id);

		CHECK(domain_tne_patch->getNormalNbrInfo(Corner<3>::bsw()).id == domain_bsw_patch->id);

		// Corner nbr id
		CHECK(domain_bsw_patch->getNormalNbrInfo(Corner<3>::tne()).rank == domain_tne_patch->rank);

		CHECK(domain_bse_patch->getNormalNbrInfo(Corner<3>::tnw()).rank == domain_tnw_patch->rank);

		CHECK(domain_bnw_patch->getNormalNbrInfo(Corner<3>::tse()).rank == domain_tse_patch->rank);

		CHECK(domain_bne_patch->getNormalNbrInfo(Corner<3>::tsw()).rank == domain_tsw_patch->rank);

		CHECK(domain_tsw_patch->getNormalNbrInfo(Corner<3>::bne()).rank == domain_bne_patch->rank);

		CHECK(domain_tse_patch->getNormalNbrInfo(Corner<3>::bnw()).rank == domain_bnw_patch->rank);

		CHECK(domain_tnw_patch->getNormalNbrInfo(Corner<3>::bse()).rank == domain_bse_patch->rank);

		CHECK(domain_tne_patch->getNormalNbrInfo(Corner<3>::bsw()).rank == domain_bsw_patch->rank);
	}
}
