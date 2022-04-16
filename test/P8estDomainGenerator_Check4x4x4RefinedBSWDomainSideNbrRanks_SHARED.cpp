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

void
Check4x4x4RefinedBSWDomainSideNeighborRanks(const PatchVector& domain)
{
  CHECK(domain["bsw_bsw_bsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bsw_bsw_bse"]->rank);
  CHECK(domain["bsw_bsw_bsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bsw_bsw_bnw"]->rank);
  CHECK(domain["bsw_bsw_bsw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bsw_bsw_tsw"]->rank);

  CHECK(domain["bsw_bsw_bse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bsw_bsw_bsw"]->rank);
  CHECK(domain["bsw_bsw_bse"]->getCoarseNbrInfo(Side<3>::east()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bsw_bsw_bse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bsw_bsw_bne"]->rank);
  CHECK(domain["bsw_bsw_bse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bsw_bsw_tse"]->rank);

  CHECK(domain["bsw_bsw_bnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bsw_bsw_bne"]->rank);
  CHECK(domain["bsw_bsw_bnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bsw_bsw_bsw"]->rank);
  CHECK(domain["bsw_bsw_bnw"]->getCoarseNbrInfo(Side<3>::north()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bsw_bsw_bnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bsw_bsw_tnw"]->rank);

  CHECK(domain["bsw_bsw_bne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bsw_bsw_bnw"]->rank);
  CHECK(domain["bsw_bsw_bne"]->getCoarseNbrInfo(Side<3>::east()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bsw_bsw_bne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bsw_bsw_bse"]->rank);
  CHECK(domain["bsw_bsw_bne"]->getCoarseNbrInfo(Side<3>::north()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bsw_bsw_bne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bsw_bsw_tne"]->rank);

  CHECK(domain["bsw_bsw_tsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bsw_bsw_tse"]->rank);
  CHECK(domain["bsw_bsw_tsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bsw_bsw_tnw"]->rank);
  CHECK(domain["bsw_bsw_tsw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bsw_bsw_bsw"]->rank);
  CHECK(domain["bsw_bsw_tsw"]->getCoarseNbrInfo(Side<3>::top()).rank == domain["bsw_tsw"]->rank);

  CHECK(domain["bsw_bsw_tse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bsw_bsw_tsw"]->rank);
  CHECK(domain["bsw_bsw_tse"]->getCoarseNbrInfo(Side<3>::east()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bsw_bsw_tse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bsw_bsw_tne"]->rank);
  CHECK(domain["bsw_bsw_tse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bsw_bsw_bse"]->rank);
  CHECK(domain["bsw_bsw_tse"]->getCoarseNbrInfo(Side<3>::top()).rank == domain["bsw_tsw"]->rank);

  CHECK(domain["bsw_bsw_tnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bsw_bsw_tne"]->rank);
  CHECK(domain["bsw_bsw_tnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bsw_bsw_tsw"]->rank);
  CHECK(domain["bsw_bsw_tnw"]->getCoarseNbrInfo(Side<3>::north()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bsw_bsw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bsw_bsw_bnw"]->rank);
  CHECK(domain["bsw_bsw_tnw"]->getCoarseNbrInfo(Side<3>::top()).rank == domain["bsw_tsw"]->rank);

  CHECK(domain["bsw_bsw_tne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bsw_bsw_tnw"]->rank);
  CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Side<3>::east()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bsw_bsw_tne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bsw_bsw_tse"]->rank);
  CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Side<3>::north()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bsw_bsw_tne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bsw_bsw_bne"]->rank);
  CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Side<3>::top()).rank == domain["bsw_tsw"]->rank);

  CHECK(domain["bsw_bse"]->getFineNbrInfo(Side<3>::west()).ranks[0] == domain["bsw_bsw_bse"]->rank);
  CHECK(domain["bsw_bse"]->getFineNbrInfo(Side<3>::west()).ranks[1] == domain["bsw_bsw_bne"]->rank);
  CHECK(domain["bsw_bse"]->getFineNbrInfo(Side<3>::west()).ranks[2] == domain["bsw_bsw_tse"]->rank);
  CHECK(domain["bsw_bse"]->getFineNbrInfo(Side<3>::west()).ranks[3] == domain["bsw_bsw_tne"]->rank);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bse_bsw"]->rank);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bsw_tse"]->rank);

  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bsw_bnw"]->getFineNbrInfo(Side<3>::south()).ranks[0] == domain["bsw_bsw_bnw"]->rank);
  CHECK(domain["bsw_bnw"]->getFineNbrInfo(Side<3>::south()).ranks[1] == domain["bsw_bsw_bne"]->rank);
  CHECK(domain["bsw_bnw"]->getFineNbrInfo(Side<3>::south()).ranks[2] == domain["bsw_bsw_tnw"]->rank);
  CHECK(domain["bsw_bnw"]->getFineNbrInfo(Side<3>::south()).ranks[3] == domain["bsw_bsw_tne"]->rank);
  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bnw_bsw"]->rank);
  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bsw_tnw"]->rank);

  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bsw_tne"]->rank);

  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["bsw_tsw"]->getFineNbrInfo(Side<3>::bottom()).ranks[0] == domain["bsw_bsw_tsw"]->rank);
  CHECK(domain["bsw_tsw"]->getFineNbrInfo(Side<3>::bottom()).ranks[1] == domain["bsw_bsw_tse"]->rank);
  CHECK(domain["bsw_tsw"]->getFineNbrInfo(Side<3>::bottom()).ranks[2] == domain["bsw_bsw_tnw"]->rank);
  CHECK(domain["bsw_tsw"]->getFineNbrInfo(Side<3>::bottom()).ranks[3] == domain["bsw_bsw_tne"]->rank);
  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tsw_bsw"]->rank);

  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bsw_tsw"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tsw_bse"]->rank);

  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bsw_tsw"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tsw_bnw"]->rank);

  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tsw_bne"]->rank);

  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bse_bse"]->rank);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bse_tsw"]->rank);

  CHECK(domain["bse_bse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bse_bsw"]->rank);
  CHECK(domain["bse_bse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bse_bne"]->rank);
  CHECK(domain["bse_bse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bse_tse"]->rank);

  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bse_bne"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bse_bsw"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bse_tnw"]->rank);

  CHECK(domain["bse_bne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bse_bse"]->rank);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bne_bse"]->rank);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bse_tne"]->rank);

  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bse_tse"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bse_bsw"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tse_bsw"]->rank);

  CHECK(domain["bse_tse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bse_tne"]->rank);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bse_bse"]->rank);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tse_bse"]->rank);

  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bse_tne"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tse_bnw"]->rank);

  CHECK(domain["bse_tne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bse_tse"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bne_tse"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bse_bne"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tse_bne"]->rank);

  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bnw_bnw"]->rank);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bnw_tsw"]->rank);

  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bnw_bsw"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bnw_bne"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bnw_tse"]->rank);

  CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bnw_bne"]->rank);
  CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bnw_bsw"]->rank);
  CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bnw_tnw"]->rank);

  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bnw_bnw"]->rank);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bne_bnw"]->rank);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bnw_tne"]->rank);

  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bnw_tnw"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bnw_bsw"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tnw_bsw"]->rank);

  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tnw_bse"]->rank);

  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bnw_bnw"]->rank);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tnw_bnw"]->rank);

  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bnw_tnw"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bnw_bne"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tnw_bne"]->rank);

  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bne_bse"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bne_bnw"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bne_tsw"]->rank);

  CHECK(domain["bne_bse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bse_bne"]->rank);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bne_bne"]->rank);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bne_tse"]->rank);

  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bnw_bne"]->rank);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bne_bne"]->rank);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bne_tnw"]->rank);

  CHECK(domain["bne_bne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bne_bnw"]->rank);
  CHECK(domain["bne_bne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bne_bse"]->rank);
  CHECK(domain["bne_bne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["bne_tne"]->rank);

  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bne_tse"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tne_bsw"]->rank);

  CHECK(domain["bne_tse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bse_tne"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["bne_tne"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bne_bse"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tne_bse"]->rank);

  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["bne_tne"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bne_bnw"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tne_bnw"]->rank);

  CHECK(domain["bne_tne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["bne_tse"]->rank);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bne_bne"]->rank);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tne_bne"]->rank);

  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bsw_tsw"]->rank);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tsw_tsw"]->rank);

  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tsw_bsw"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tsw_tse"]->rank);

  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tsw_bsw"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tsw_tnw"]->rank);

  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tsw_tne"]->rank);

  CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tsw_tse"]->rank);
  CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tsw_tnw"]->rank);
  CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tsw_bsw"]->rank);

  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tsw_tsw"]->rank);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tse_tsw"]->rank);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tsw_tne"]->rank);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tsw_bse"]->rank);

  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tsw_tne"]->rank);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tsw_tsw"]->rank);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tnw_tsw"]->rank);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tsw_bnw"]->rank);

  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tsw_tnw"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tsw_tse"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tnw_tse"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tsw_bne"]->rank);

  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tse_bse"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tse_tsw"]->rank);

  CHECK(domain["tse_bse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tse_bne"]->rank);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bse_tse"]->rank);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tse_tse"]->rank);

  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tse_bne"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tse_tnw"]->rank);

  CHECK(domain["tse_bne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tse_bse"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tne_bse"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bse_tne"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tse_tne"]->rank);

  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tsw_tse"]->rank);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tse_tse"]->rank);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tse_bsw"]->rank);

  CHECK(domain["tse_tse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tse_tsw"]->rank);
  CHECK(domain["tse_tse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tse_tne"]->rank);
  CHECK(domain["tse_tse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tse_bse"]->rank);

  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tsw_tne"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tse_tne"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tse_tsw"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tne_tsw"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tse_bnw"]->rank);

  CHECK(domain["tse_tne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tse_tse"]->rank);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tne_tse"]->rank);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tse_bne"]->rank);

  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tnw_bnw"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tnw_tsw"]->rank);

  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tnw_bne"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tnw_tse"]->rank);

  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tnw_bne"]->rank);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bnw_tnw"]->rank);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tnw_tnw"]->rank);

  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tnw_bnw"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tne_bnw"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tnw_tne"]->rank);

  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tnw_tse"]->rank);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tsw_tnw"]->rank);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tnw_tnw"]->rank);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tnw_bsw"]->rank);

  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tnw_tsw"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tne_tsw"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tsw_tne"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tnw_tne"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tnw_bse"]->rank);

  CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tnw_tne"]->rank);
  CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tnw_tsw"]->rank);
  CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tnw_bnw"]->rank);

  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tnw_tnw"]->rank);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tne_tnw"]->rank);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tnw_tse"]->rank);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tnw_bne"]->rank);

  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tne_bse"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tne_bnw"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tne_tsw"]->rank);

  CHECK(domain["tne_bse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tse_bne"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tne_bne"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bne_tse"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tne_tse"]->rank);

  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tnw_bne"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tne_bne"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tne_tnw"]->rank);

  CHECK(domain["tne_bne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tne_bnw"]->rank);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tne_bse"]->rank);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["bne_tne"]->rank);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Side<3>::top()).rank == domain["tne_tne"]->rank);

  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tnw_tse"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tne_tse"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tne_tnw"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tne_bsw"]->rank);

  CHECK(domain["tne_tse"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tne_tsw"]->rank);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tse_tne"]->rank);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Side<3>::north()).rank == domain["tne_tne"]->rank);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tne_bse"]->rank);

  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tnw_tne"]->rank);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Side<3>::east()).rank == domain["tne_tne"]->rank);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tne_tsw"]->rank);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tne_bnw"]->rank);

  CHECK(domain["tne_tne"]->getNormalNbrInfo(Side<3>::west()).rank == domain["tne_tnw"]->rank);
  CHECK(domain["tne_tne"]->getNormalNbrInfo(Side<3>::south()).rank == domain["tne_tse"]->rank);
  CHECK(domain["tne_tne"]->getNormalNbrInfo(Side<3>::bottom()).rank == domain["tne_bne"]->rank);
}
