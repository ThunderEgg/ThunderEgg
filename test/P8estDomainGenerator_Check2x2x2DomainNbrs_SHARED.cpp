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

using namespace std;
using namespace ThunderEgg;

#include <doctest.h>

void
Check2x2x2DomainNeighbors(const ThunderEgg::Domain<3>& domain)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  vector<PatchInfo<3>> patches = GetAllPatchesOnRank0(domain);
  if (rank == 0) {
    const PatchInfo<3>* domain_bsw_patch = nullptr;
    const PatchInfo<3>* domain_bse_patch = nullptr;
    const PatchInfo<3>* domain_bnw_patch = nullptr;
    const PatchInfo<3>* domain_bne_patch = nullptr;
    const PatchInfo<3>* domain_tsw_patch = nullptr;
    const PatchInfo<3>* domain_tse_patch = nullptr;
    const PatchInfo<3>* domain_tnw_patch = nullptr;
    const PatchInfo<3>* domain_tne_patch = nullptr;

    for (PatchInfo<3>& patch : patches) {
      double x = patch.starts[0];
      double y = patch.starts[1];
      double z = patch.starts[2];
      if (x == doctest::Approx(0) && y == doctest::Approx(0) && z == doctest::Approx(0)) {
        domain_bsw_patch = &patch;
      }
      if (x == doctest::Approx(0.5) && y == doctest::Approx(0) && z == doctest::Approx(0)) {
        domain_bse_patch = &patch;
      }
      if (x == doctest::Approx(0) && y == doctest::Approx(0.5) && z == doctest::Approx(0)) {
        domain_bnw_patch = &patch;
      }
      if (x == doctest::Approx(0.5) && y == doctest::Approx(0.5) && z == doctest::Approx(0)) {
        domain_bne_patch = &patch;
      }
      if (x == doctest::Approx(0) && y == doctest::Approx(0) && z == doctest::Approx(0.5)) {
        domain_tsw_patch = &patch;
      }
      if (x == doctest::Approx(0.5) && y == doctest::Approx(0) && z == doctest::Approx(0.5)) {
        domain_tse_patch = &patch;
      }
      if (x == doctest::Approx(0) && y == doctest::Approx(0.5) && z == doctest::Approx(0.5)) {
        domain_tnw_patch = &patch;
      }
      if (x == doctest::Approx(0.5) && y == doctest::Approx(0.5) && z == doctest::Approx(0.5)) {
        domain_tne_patch = &patch;
      }
    }

    REQUIRE_NE(domain_bsw_patch, nullptr);
    REQUIRE_NE(domain_bse_patch, nullptr);
    REQUIRE_NE(domain_bnw_patch, nullptr);
    REQUIRE_NE(domain_bne_patch, nullptr);
    REQUIRE_NE(domain_tsw_patch, nullptr);
    REQUIRE_NE(domain_tse_patch, nullptr);
    REQUIRE_NE(domain_tnw_patch, nullptr);
    REQUIRE_NE(domain_tne_patch, nullptr);

    // Side hasnbr
    CHECK_EQ(domain_bsw_patch->hasNbr(Side<3>::west()), false);
    CHECK_EQ(domain_bsw_patch->hasNbr(Side<3>::east()), true);
    CHECK_EQ(domain_bsw_patch->hasNbr(Side<3>::south()), false);
    CHECK_EQ(domain_bsw_patch->hasNbr(Side<3>::north()), true);
    CHECK_EQ(domain_bsw_patch->hasNbr(Side<3>::bottom()), false);
    CHECK_EQ(domain_bsw_patch->hasNbr(Side<3>::top()), true);

    CHECK_EQ(domain_bse_patch->hasNbr(Side<3>::west()), true);
    CHECK_EQ(domain_bse_patch->hasNbr(Side<3>::east()), false);
    CHECK_EQ(domain_bse_patch->hasNbr(Side<3>::south()), false);
    CHECK_EQ(domain_bse_patch->hasNbr(Side<3>::north()), true);
    CHECK_EQ(domain_bse_patch->hasNbr(Side<3>::bottom()), false);
    CHECK_EQ(domain_bse_patch->hasNbr(Side<3>::top()), true);

    CHECK_EQ(domain_bnw_patch->hasNbr(Side<3>::west()), false);
    CHECK_EQ(domain_bnw_patch->hasNbr(Side<3>::east()), true);
    CHECK_EQ(domain_bnw_patch->hasNbr(Side<3>::south()), true);
    CHECK_EQ(domain_bnw_patch->hasNbr(Side<3>::north()), false);
    CHECK_EQ(domain_bnw_patch->hasNbr(Side<3>::bottom()), false);
    CHECK_EQ(domain_bnw_patch->hasNbr(Side<3>::top()), true);

    CHECK_EQ(domain_bne_patch->hasNbr(Side<3>::west()), true);
    CHECK_EQ(domain_bne_patch->hasNbr(Side<3>::east()), false);
    CHECK_EQ(domain_bne_patch->hasNbr(Side<3>::south()), true);
    CHECK_EQ(domain_bne_patch->hasNbr(Side<3>::north()), false);
    CHECK_EQ(domain_bne_patch->hasNbr(Side<3>::bottom()), false);
    CHECK_EQ(domain_bne_patch->hasNbr(Side<3>::top()), true);

    CHECK_EQ(domain_tsw_patch->hasNbr(Side<3>::west()), false);
    CHECK_EQ(domain_tsw_patch->hasNbr(Side<3>::east()), true);
    CHECK_EQ(domain_tsw_patch->hasNbr(Side<3>::south()), false);
    CHECK_EQ(domain_tsw_patch->hasNbr(Side<3>::north()), true);
    CHECK_EQ(domain_tsw_patch->hasNbr(Side<3>::bottom()), true);
    CHECK_EQ(domain_tsw_patch->hasNbr(Side<3>::top()), false);

    CHECK_EQ(domain_tse_patch->hasNbr(Side<3>::west()), true);
    CHECK_EQ(domain_tse_patch->hasNbr(Side<3>::east()), false);
    CHECK_EQ(domain_tse_patch->hasNbr(Side<3>::south()), false);
    CHECK_EQ(domain_tse_patch->hasNbr(Side<3>::north()), true);
    CHECK_EQ(domain_tse_patch->hasNbr(Side<3>::bottom()), true);
    CHECK_EQ(domain_tse_patch->hasNbr(Side<3>::top()), false);

    CHECK_EQ(domain_tnw_patch->hasNbr(Side<3>::west()), false);
    CHECK_EQ(domain_tnw_patch->hasNbr(Side<3>::east()), true);
    CHECK_EQ(domain_tnw_patch->hasNbr(Side<3>::south()), true);
    CHECK_EQ(domain_tnw_patch->hasNbr(Side<3>::north()), false);
    CHECK_EQ(domain_tnw_patch->hasNbr(Side<3>::bottom()), true);
    CHECK_EQ(domain_tnw_patch->hasNbr(Side<3>::top()), false);

    CHECK_EQ(domain_tne_patch->hasNbr(Side<3>::west()), true);
    CHECK_EQ(domain_tne_patch->hasNbr(Side<3>::east()), false);
    CHECK_EQ(domain_tne_patch->hasNbr(Side<3>::south()), true);
    CHECK_EQ(domain_tne_patch->hasNbr(Side<3>::north()), false);
    CHECK_EQ(domain_tne_patch->hasNbr(Side<3>::bottom()), true);
    CHECK_EQ(domain_tne_patch->hasNbr(Side<3>::top()), false);

    // SIDE ids
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Side<3>::east()).id, domain_bse_patch->id);
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Side<3>::north()).id, domain_bnw_patch->id);
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Side<3>::top()).id, domain_tsw_patch->id);

    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Side<3>::west()).id, domain_bsw_patch->id);
    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Side<3>::north()).id, domain_bne_patch->id);
    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Side<3>::top()).id, domain_tse_patch->id);

    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Side<3>::east()).id, domain_bne_patch->id);
    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Side<3>::south()).id, domain_bsw_patch->id);
    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Side<3>::top()).id, domain_tnw_patch->id);

    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Side<3>::west()).id, domain_bnw_patch->id);
    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Side<3>::south()).id, domain_bse_patch->id);
    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Side<3>::top()).id, domain_tne_patch->id);

    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Side<3>::east()).id, domain_tse_patch->id);
    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Side<3>::north()).id, domain_tnw_patch->id);
    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Side<3>::bottom()).id, domain_bsw_patch->id);

    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Side<3>::west()).id, domain_tsw_patch->id);
    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Side<3>::north()).id, domain_tne_patch->id);
    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Side<3>::bottom()).id, domain_bse_patch->id);

    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Side<3>::east()).id, domain_tne_patch->id);
    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Side<3>::south()).id, domain_tsw_patch->id);
    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Side<3>::bottom()).id, domain_bnw_patch->id);

    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Side<3>::west()).id, domain_tnw_patch->id);
    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Side<3>::south()).id, domain_tse_patch->id);
    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Side<3>::bottom()).id, domain_bne_patch->id);

    // SIDE ranks
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Side<3>::east()).rank, domain_bse_patch->rank);
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Side<3>::north()).rank, domain_bnw_patch->rank);
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Side<3>::top()).rank, domain_tsw_patch->rank);

    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Side<3>::west()).rank, domain_bsw_patch->rank);
    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Side<3>::north()).rank, domain_bne_patch->rank);
    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Side<3>::top()).rank, domain_tse_patch->rank);

    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Side<3>::east()).rank, domain_bne_patch->rank);
    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Side<3>::south()).rank, domain_bsw_patch->rank);
    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Side<3>::top()).rank, domain_tnw_patch->rank);

    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Side<3>::west()).rank, domain_bnw_patch->rank);
    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Side<3>::south()).rank, domain_bse_patch->rank);
    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Side<3>::top()).rank, domain_tne_patch->rank);

    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Side<3>::east()).rank, domain_tse_patch->rank);
    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Side<3>::north()).rank, domain_tnw_patch->rank);
    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Side<3>::bottom()).rank, domain_bsw_patch->rank);

    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Side<3>::west()).rank, domain_tsw_patch->rank);
    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Side<3>::north()).rank, domain_tne_patch->rank);
    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Side<3>::bottom()).rank, domain_bse_patch->rank);

    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Side<3>::east()).rank, domain_tne_patch->rank);
    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Side<3>::south()).rank, domain_tsw_patch->rank);
    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Side<3>::bottom()).rank, domain_bnw_patch->rank);

    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Side<3>::west()).rank, domain_tnw_patch->rank);
    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Side<3>::south()).rank, domain_tse_patch->rank);
    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Side<3>::bottom()).rank, domain_bne_patch->rank);

    // Edge hasnbr
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Edge::bs()));
    CHECK_UNARY(domain_bsw_patch->hasNbr(Edge::tn()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Edge::bn()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Edge::ts()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Edge::bw()));
    CHECK_UNARY(domain_bsw_patch->hasNbr(Edge::te()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Edge::be()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Edge::tw()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Edge::sw()));
    CHECK_UNARY(domain_bsw_patch->hasNbr(Edge::ne()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Edge::se()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Edge::nw()));

    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Edge::bs()));
    CHECK_UNARY(domain_bse_patch->hasNbr(Edge::tn()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Edge::bn()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Edge::ts()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Edge::bw()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Edge::te()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Edge::be()));
    CHECK_UNARY(domain_bse_patch->hasNbr(Edge::tw()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Edge::sw()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Edge::ne()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Edge::se()));
    CHECK_UNARY(domain_bse_patch->hasNbr(Edge::nw()));

    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Edge::bs()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Edge::tn()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Edge::bn()));
    CHECK_UNARY(domain_bnw_patch->hasNbr(Edge::ts()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Edge::bw()));
    CHECK_UNARY(domain_bnw_patch->hasNbr(Edge::te()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Edge::be()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Edge::tw()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Edge::sw()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Edge::ne()));
    CHECK_UNARY(domain_bnw_patch->hasNbr(Edge::se()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Edge::nw()));

    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Edge::bs()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Edge::tn()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Edge::bn()));
    CHECK_UNARY(domain_bne_patch->hasNbr(Edge::ts()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Edge::bw()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Edge::te()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Edge::be()));
    CHECK_UNARY(domain_bne_patch->hasNbr(Edge::tw()));
    CHECK_UNARY(domain_bne_patch->hasNbr(Edge::sw()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Edge::ne()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Edge::se()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Edge::nw()));

    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Edge::bs()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Edge::tn()));
    CHECK_UNARY(domain_tsw_patch->hasNbr(Edge::bn()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Edge::ts()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Edge::bw()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Edge::te()));
    CHECK_UNARY(domain_tsw_patch->hasNbr(Edge::be()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Edge::tw()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Edge::sw()));
    CHECK_UNARY(domain_tsw_patch->hasNbr(Edge::ne()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Edge::se()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Edge::nw()));

    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Edge::bs()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Edge::tn()));
    CHECK_UNARY(domain_tse_patch->hasNbr(Edge::bn()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Edge::ts()));
    CHECK_UNARY(domain_tse_patch->hasNbr(Edge::bw()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Edge::te()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Edge::be()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Edge::tw()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Edge::sw()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Edge::ne()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Edge::se()));
    CHECK_UNARY(domain_tse_patch->hasNbr(Edge::nw()));

    CHECK_UNARY(domain_tnw_patch->hasNbr(Edge::bs()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Edge::tn()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Edge::bn()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Edge::ts()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Edge::bw()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Edge::te()));
    CHECK_UNARY(domain_tnw_patch->hasNbr(Edge::be()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Edge::tw()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Edge::sw()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Edge::ne()));
    CHECK_UNARY(domain_tnw_patch->hasNbr(Edge::se()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Edge::nw()));

    CHECK_UNARY(domain_tne_patch->hasNbr(Edge::bs()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Edge::tn()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Edge::bn()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Edge::ts()));
    CHECK_UNARY(domain_tne_patch->hasNbr(Edge::bw()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Edge::te()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Edge::be()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Edge::tw()));
    CHECK_UNARY(domain_tne_patch->hasNbr(Edge::sw()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Edge::ne()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Edge::se()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Edge::nw()));

    // Edge nbr id
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Edge::tn()).id, domain_tnw_patch->id);
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Edge::te()).id, domain_tse_patch->id);
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Edge::ne()).id, domain_bne_patch->id);

    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Edge::tn()).id, domain_tne_patch->id);
    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Edge::tw()).id, domain_tsw_patch->id);
    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Edge::nw()).id, domain_bnw_patch->id);

    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Edge::ts()).id, domain_tsw_patch->id);
    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Edge::te()).id, domain_tne_patch->id);
    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Edge::se()).id, domain_bse_patch->id);

    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Edge::ts()).id, domain_tse_patch->id);
    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Edge::tw()).id, domain_tnw_patch->id);
    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Edge::sw()).id, domain_bsw_patch->id);

    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Edge::bn()).id, domain_bnw_patch->id);
    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Edge::be()).id, domain_bse_patch->id);
    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Edge::ne()).id, domain_tne_patch->id);

    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Edge::bn()).id, domain_bne_patch->id);
    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Edge::bw()).id, domain_bsw_patch->id);
    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Edge::nw()).id, domain_tnw_patch->id);

    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Edge::bs()).id, domain_bsw_patch->id);
    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Edge::be()).id, domain_bne_patch->id);
    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Edge::se()).id, domain_tse_patch->id);

    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Edge::bs()).id, domain_bse_patch->id);
    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Edge::bw()).id, domain_bnw_patch->id);
    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Edge::sw()).id, domain_tsw_patch->id);

    // Edge nbr rank
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Edge::tn()).rank, domain_tnw_patch->rank);
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Edge::te()).rank, domain_tse_patch->rank);
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Edge::ne()).rank, domain_bne_patch->rank);

    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Edge::tn()).rank, domain_tne_patch->rank);
    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Edge::tw()).rank, domain_tsw_patch->rank);
    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Edge::nw()).rank, domain_bnw_patch->rank);

    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Edge::ts()).rank, domain_tsw_patch->rank);
    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Edge::te()).rank, domain_tne_patch->rank);
    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Edge::se()).rank, domain_bse_patch->rank);

    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Edge::ts()).rank, domain_tse_patch->rank);
    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Edge::tw()).rank, domain_tnw_patch->rank);
    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Edge::sw()).rank, domain_bsw_patch->rank);

    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Edge::bn()).rank, domain_bnw_patch->rank);
    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Edge::be()).rank, domain_bse_patch->rank);
    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Edge::ne()).rank, domain_tne_patch->rank);

    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Edge::bn()).rank, domain_bne_patch->rank);
    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Edge::bw()).rank, domain_bsw_patch->rank);
    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Edge::nw()).rank, domain_tnw_patch->rank);

    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Edge::bs()).rank, domain_bsw_patch->rank);
    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Edge::be()).rank, domain_bne_patch->rank);
    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Edge::se()).rank, domain_tse_patch->rank);

    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Edge::bs()).rank, domain_bse_patch->rank);
    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Edge::bw()).rank, domain_bnw_patch->rank);
    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Edge::sw()).rank, domain_tsw_patch->rank);

    // Corner hasNbr
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Corner<3>::bsw()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Corner<3>::bse()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Corner<3>::bnw()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Corner<3>::bne()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Corner<3>::tsw()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Corner<3>::tse()));
    CHECK_UNARY_FALSE(domain_bsw_patch->hasNbr(Corner<3>::tnw()));
    CHECK_UNARY(domain_bsw_patch->hasNbr(Corner<3>::tne()));

    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Corner<3>::bsw()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Corner<3>::bse()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Corner<3>::bnw()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Corner<3>::bne()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Corner<3>::tsw()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Corner<3>::tse()));
    CHECK_UNARY(domain_bse_patch->hasNbr(Corner<3>::tnw()));
    CHECK_UNARY_FALSE(domain_bse_patch->hasNbr(Corner<3>::tne()));

    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Corner<3>::bsw()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Corner<3>::bse()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Corner<3>::bnw()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Corner<3>::bne()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Corner<3>::tsw()));
    CHECK_UNARY(domain_bnw_patch->hasNbr(Corner<3>::tse()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Corner<3>::tnw()));
    CHECK_UNARY_FALSE(domain_bnw_patch->hasNbr(Corner<3>::tne()));

    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Corner<3>::bsw()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Corner<3>::bse()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Corner<3>::bnw()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Corner<3>::bne()));
    CHECK_UNARY(domain_bne_patch->hasNbr(Corner<3>::tsw()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Corner<3>::tse()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Corner<3>::tnw()));
    CHECK_UNARY_FALSE(domain_bne_patch->hasNbr(Corner<3>::tne()));

    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Corner<3>::bsw()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Corner<3>::bse()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Corner<3>::bnw()));
    CHECK_UNARY(domain_tsw_patch->hasNbr(Corner<3>::bne()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Corner<3>::tsw()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Corner<3>::tse()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Corner<3>::tnw()));
    CHECK_UNARY_FALSE(domain_tsw_patch->hasNbr(Corner<3>::tne()));

    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Corner<3>::bsw()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Corner<3>::bse()));
    CHECK_UNARY(domain_tse_patch->hasNbr(Corner<3>::bnw()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Corner<3>::bne()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Corner<3>::tsw()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Corner<3>::tse()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Corner<3>::tnw()));
    CHECK_UNARY_FALSE(domain_tse_patch->hasNbr(Corner<3>::tne()));

    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Corner<3>::bsw()));
    CHECK_UNARY(domain_tnw_patch->hasNbr(Corner<3>::bse()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Corner<3>::bnw()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Corner<3>::bne()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Corner<3>::tsw()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Corner<3>::tse()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Corner<3>::tnw()));
    CHECK_UNARY_FALSE(domain_tnw_patch->hasNbr(Corner<3>::tne()));

    CHECK_UNARY(domain_tne_patch->hasNbr(Corner<3>::bsw()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Corner<3>::bse()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Corner<3>::bnw()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Corner<3>::bne()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Corner<3>::tsw()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Corner<3>::tse()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Corner<3>::tnw()));
    CHECK_UNARY_FALSE(domain_tne_patch->hasNbr(Corner<3>::tne()));

    // Corner nbr id
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Corner<3>::tne()).id, domain_tne_patch->id);

    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Corner<3>::tnw()).id, domain_tnw_patch->id);

    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Corner<3>::tse()).id, domain_tse_patch->id);

    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Corner<3>::tsw()).id, domain_tsw_patch->id);

    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Corner<3>::bne()).id, domain_bne_patch->id);

    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Corner<3>::bnw()).id, domain_bnw_patch->id);

    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Corner<3>::bse()).id, domain_bse_patch->id);

    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Corner<3>::bsw()).id, domain_bsw_patch->id);

    // Corner nbr id
    CHECK_EQ(domain_bsw_patch->getNormalNbrInfo(Corner<3>::tne()).rank, domain_tne_patch->rank);

    CHECK_EQ(domain_bse_patch->getNormalNbrInfo(Corner<3>::tnw()).rank, domain_tnw_patch->rank);

    CHECK_EQ(domain_bnw_patch->getNormalNbrInfo(Corner<3>::tse()).rank, domain_tse_patch->rank);

    CHECK_EQ(domain_bne_patch->getNormalNbrInfo(Corner<3>::tsw()).rank, domain_tsw_patch->rank);

    CHECK_EQ(domain_tsw_patch->getNormalNbrInfo(Corner<3>::bne()).rank, domain_bne_patch->rank);

    CHECK_EQ(domain_tse_patch->getNormalNbrInfo(Corner<3>::bnw()).rank, domain_bnw_patch->rank);

    CHECK_EQ(domain_tnw_patch->getNormalNbrInfo(Corner<3>::bse()).rank, domain_bse_patch->rank);

    CHECK_EQ(domain_tne_patch->getNormalNbrInfo(Corner<3>::bsw()).rank, domain_bsw_patch->rank);
  }
}
