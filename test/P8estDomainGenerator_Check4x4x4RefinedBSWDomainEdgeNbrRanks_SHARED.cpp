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
Check4x4x4RefinedBSWDomainEdgeNeighborRanks(const PatchVector& domain)
{
  CHECK(domain["bsw_bsw_bsw"]->getNormalNbrInfo(Edge::tn()).rank == domain["bsw_bsw_tnw"]->rank);
  CHECK(domain["bsw_bsw_bsw"]->getNormalNbrInfo(Edge::te()).rank == domain["bsw_bsw_tse"]->rank);
  CHECK(domain["bsw_bsw_bsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bsw_bsw_bne"]->rank);

  CHECK(domain["bsw_bsw_bse"]->getNormalNbrInfo(Edge::tn()).rank == domain["bsw_bsw_tne"]->rank);
  CHECK(domain["bsw_bsw_bse"]->getNormalNbrInfo(Edge::tw()).rank == domain["bsw_bsw_tsw"]->rank);
  CHECK(domain["bsw_bsw_bse"]->getNormalNbrInfo(Edge::nw()).rank == domain["bsw_bsw_bnw"]->rank);

  CHECK(domain["bsw_bsw_bnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["bsw_bsw_tsw"]->rank);
  CHECK(domain["bsw_bsw_bnw"]->getNormalNbrInfo(Edge::te()).rank == domain["bsw_bsw_tne"]->rank);
  CHECK(domain["bsw_bsw_bnw"]->getNormalNbrInfo(Edge::se()).rank == domain["bsw_bsw_bse"]->rank);

  CHECK(domain["bsw_bsw_bne"]->getNormalNbrInfo(Edge::ts()).rank == domain["bsw_bsw_tse"]->rank);
  CHECK(domain["bsw_bsw_bne"]->getNormalNbrInfo(Edge::tw()).rank == domain["bsw_bsw_tnw"]->rank);
  CHECK(domain["bsw_bsw_bne"]->getNormalNbrInfo(Edge::sw()).rank == domain["bsw_bsw_bsw"]->rank);
  CHECK(domain["bsw_bsw_bne"]->getCoarseNbrInfo(Edge::ne()).rank == domain["bsw_bne"]->rank);

  CHECK(domain["bsw_bsw_tsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bsw_bsw_bnw"]->rank);
  CHECK(domain["bsw_bsw_tsw"]->getNormalNbrInfo(Edge::be()).rank == domain["bsw_bsw_bse"]->rank);
  CHECK(domain["bsw_bsw_tsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bsw_bsw_tne"]->rank);

  CHECK(domain["bsw_bsw_tse"]->getNormalNbrInfo(Edge::bn()).rank == domain["bsw_bsw_bne"]->rank);
  CHECK(domain["bsw_bsw_tse"]->getNormalNbrInfo(Edge::bw()).rank == domain["bsw_bsw_bsw"]->rank);
  CHECK(domain["bsw_bsw_tse"]->getCoarseNbrInfo(Edge::te()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["bsw_bsw_tse"]->getNormalNbrInfo(Edge::nw()).rank == domain["bsw_bsw_tnw"]->rank);

  CHECK(domain["bsw_bsw_tnw"]->getNormalNbrInfo(Edge::bs()).rank == domain["bsw_bsw_bsw"]->rank);
  CHECK(domain["bsw_bsw_tnw"]->getCoarseNbrInfo(Edge::tn()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["bsw_bsw_tnw"]->getNormalNbrInfo(Edge::be()).rank == domain["bsw_bsw_bne"]->rank);
  CHECK(domain["bsw_bsw_tnw"]->getNormalNbrInfo(Edge::se()).rank == domain["bsw_bsw_tse"]->rank);

  CHECK(domain["bsw_bsw_tne"]->getNormalNbrInfo(Edge::bs()).rank == domain["bsw_bsw_bse"]->rank);
  CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Edge::tn()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["bsw_bsw_tne"]->getNormalNbrInfo(Edge::bw()).rank == domain["bsw_bsw_bnw"]->rank);
  CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Edge::te()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["bsw_bsw_tne"]->getNormalNbrInfo(Edge::sw()).rank == domain["bsw_bsw_tsw"]->rank);
  CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Edge::ne()).rank == domain["bsw_bne"]->rank);

  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Edge::tn()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Edge::te()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Edge::tw()).rank == domain["bsw_tsw"]->rank);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Edge::ne()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Edge::nw()).rank == domain["bsw_bnw"]->rank);

  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Edge::tn()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["bsw_tsw"]->rank);
  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Edge::te()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Edge::se()).rank == domain["bsw_bse"]->rank);

  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::tn()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::ts()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::te()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::tw()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["bsw_bne"]->getFineNbrInfo(Edge::sw()).ranks[0] == domain["bsw_bsw_bne"]->rank);
  CHECK(domain["bsw_bne"]->getFineNbrInfo(Edge::sw()).ranks[1] == domain["bsw_bsw_tne"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::ne()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::se()).rank == domain["bse_bsw"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::nw()).rank == domain["bnw_bsw"]->rank);

  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Edge::tn()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Edge::te()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Edge::be()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bsw_tne"]->rank);

  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::tn()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::bn()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bsw_tse"]->getFineNbrInfo(Edge::bw()).ranks[0] == domain["bsw_bsw_tse"]->rank);
  CHECK(domain["bsw_tse"]->getFineNbrInfo(Edge::bw()).ranks[1] == domain["bsw_bsw_tne"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::te()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::be()).rank == domain["bse_bsw"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::tw()).rank == domain["tsw_bsw"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::ne()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::nw()).rank == domain["bsw_tnw"]->rank);

  CHECK(domain["bsw_tnw"]->getFineNbrInfo(Edge::bs()).ranks[0] == domain["bsw_bsw_tnw"]->rank);
  CHECK(domain["bsw_tnw"]->getFineNbrInfo(Edge::bs()).ranks[1] == domain["bsw_bsw_tne"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::tn()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bnw_bsw"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["tsw_bsw"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::te()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::be()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::se()).rank == domain["bsw_tse"]->rank);

  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::bs()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::tn()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::bn()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::ts()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::bw()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::te()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::be()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::tw()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::sw()).rank == domain["bsw_tsw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::ne()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::se()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::nw()).rank == domain["bnw_tsw"]->rank);

  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Edge::tn()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Edge::te()).rank == domain["bse_tse"]->rank);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Edge::tw()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bse_bne"]->rank);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Edge::nw()).rank == domain["bsw_bne"]->rank);

  CHECK(domain["bse_bse"]->getNormalNbrInfo(Edge::tn()).rank == domain["bse_tne"]->rank);
  CHECK(domain["bse_bse"]->getNormalNbrInfo(Edge::tw()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["bse_bse"]->getNormalNbrInfo(Edge::nw()).rank == domain["bse_bnw"]->rank);

  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::tn()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::te()).rank == domain["bse_tne"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::tw()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::sw()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bne_bse"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::se()).rank == domain["bse_bse"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::nw()).rank == domain["bnw_bse"]->rank);

  CHECK(domain["bse_bne"]->getNormalNbrInfo(Edge::tn()).rank == domain["bne_tse"]->rank);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Edge::ts()).rank == domain["bse_tse"]->rank);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Edge::tw()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Edge::sw()).rank == domain["bse_bsw"]->rank);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Edge::nw()).rank == domain["bne_bsw"]->rank);

  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::tn()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::bw()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::te()).rank == domain["tse_bse"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::be()).rank == domain["bse_bse"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::tw()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bse_tne"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::nw()).rank == domain["bsw_tne"]->rank);

  CHECK(domain["bse_tse"]->getNormalNbrInfo(Edge::tn()).rank == domain["tse_bne"]->rank);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Edge::bn()).rank == domain["bse_bne"]->rank);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Edge::bw()).rank == domain["bse_bsw"]->rank);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Edge::tw()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Edge::nw()).rank == domain["bse_tnw"]->rank);

  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::bs()).rank == domain["bse_bsw"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::tn()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::bw()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::te()).rank == domain["tse_bne"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::be()).rank == domain["bse_bne"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::tw()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::sw()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bne_tse"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::se()).rank == domain["bse_tse"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::nw()).rank == domain["bnw_tse"]->rank);

  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::bs()).rank == domain["bse_bse"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::tn()).rank == domain["tne_bse"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::bn()).rank == domain["bne_bse"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::ts()).rank == domain["tse_bse"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::bw()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::tw()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::sw()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::nw()).rank == domain["bne_tsw"]->rank);

  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Edge::tn()).rank == domain["bnw_tnw"]->rank);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Edge::ts()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Edge::te()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bnw_bne"]->rank);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Edge::se()).rank == domain["bsw_bne"]->rank);

  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::tn()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::ts()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::te()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::tw()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::sw()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::ne()).rank == domain["bne_bnw"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::se()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::nw()).rank == domain["bnw_bnw"]->rank);

  CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Edge::te()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Edge::se()).rank == domain["bnw_bse"]->rank);

  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Edge::ts()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Edge::te()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Edge::tw()).rank == domain["bnw_tnw"]->rank);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Edge::sw()).rank == domain["bnw_bsw"]->rank);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Edge::se()).rank == domain["bne_bsw"]->rank);

  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::bs()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::tn()).rank == domain["tnw_bnw"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bnw_bnw"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::ts()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::te()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::be()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::se()).rank == domain["bsw_tne"]->rank);

  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::bs()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::tn()).rank == domain["tnw_bne"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::bn()).rank == domain["bnw_bne"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::ts()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::bw()).rank == domain["bnw_bsw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::te()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::be()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::tw()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::sw()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::ne()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::se()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::nw()).rank == domain["bnw_tnw"]->rank);

  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Edge::bs()).rank == domain["bnw_bsw"]->rank);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Edge::te()).rank == domain["tnw_bne"]->rank);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Edge::be()).rank == domain["bnw_bne"]->rank);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Edge::se()).rank == domain["bnw_tse"]->rank);

  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::bs()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::ts()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::bw()).rank == domain["bnw_bnw"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::te()).rank == domain["tne_bnw"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::be()).rank == domain["bne_bnw"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::tw()).rank == domain["tnw_bnw"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::sw()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::se()).rank == domain["bne_tsw"]->rank);

  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::tn()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::ts()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::te()).rank == domain["bne_tse"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::tw()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::sw()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bne_bne"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::se()).rank == domain["bse_bne"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::nw()).rank == domain["bnw_bne"]->rank);

  CHECK(domain["bne_bse"]->getNormalNbrInfo(Edge::tn()).rank == domain["bne_tne"]->rank);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Edge::ts()).rank == domain["bse_tne"]->rank);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Edge::tw()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Edge::sw()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Edge::nw()).rank == domain["bne_bnw"]->rank);

  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Edge::te()).rank == domain["bne_tne"]->rank);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Edge::tw()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Edge::sw()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Edge::se()).rank == domain["bne_bse"]->rank);

  CHECK(domain["bne_bne"]->getNormalNbrInfo(Edge::ts()).rank == domain["bne_tse"]->rank);
  CHECK(domain["bne_bne"]->getNormalNbrInfo(Edge::tw()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["bne_bne"]->getNormalNbrInfo(Edge::sw()).rank == domain["bne_bsw"]->rank);

  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::bs()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::tn()).rank == domain["tne_bnw"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bne_bnw"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::ts()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::bw()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::te()).rank == domain["tne_bse"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::be()).rank == domain["bne_bse"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::tw()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::sw()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["bne_tne"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::se()).rank == domain["bse_tne"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::nw()).rank == domain["bnw_tne"]->rank);

  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::bs()).rank == domain["bse_bne"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::tn()).rank == domain["tne_bne"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::bn()).rank == domain["bne_bne"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::ts()).rank == domain["tse_bne"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::bw()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::tw()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::sw()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::nw()).rank == domain["bne_tnw"]->rank);

  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::bs()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::bw()).rank == domain["bnw_bne"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::te()).rank == domain["tne_bne"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::be()).rank == domain["bne_bne"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::tw()).rank == domain["tnw_bne"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::sw()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::se()).rank == domain["bne_tse"]->rank);

  CHECK(domain["bne_tne"]->getNormalNbrInfo(Edge::bs()).rank == domain["bne_bse"]->rank);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Edge::ts()).rank == domain["tne_bse"]->rank);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Edge::bw()).rank == domain["bne_bnw"]->rank);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Edge::tw()).rank == domain["tne_bnw"]->rank);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Edge::sw()).rank == domain["bne_tsw"]->rank);

  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Edge::tn()).rank == domain["tsw_tnw"]->rank);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Edge::te()).rank == domain["tsw_tse"]->rank);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Edge::be()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["tsw_bne"]->rank);

  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::tn()).rank == domain["tsw_tne"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::bn()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::bw()).rank == domain["bsw_tsw"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::te()).rank == domain["tse_tsw"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::be()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::tw()).rank == domain["tsw_tsw"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::ne()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::nw()).rank == domain["tsw_bnw"]->rank);

  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::bs()).rank == domain["bsw_tsw"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::tn()).rank == domain["tnw_tsw"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["tsw_tsw"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::te()).rank == domain["tsw_tne"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::be()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::ne()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::se()).rank == domain["tsw_bse"]->rank);

  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::bs()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::tn()).rank == domain["tnw_tse"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::bn()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::ts()).rank == domain["tsw_tse"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::bw()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::te()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::be()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::tw()).rank == domain["tsw_tnw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::sw()).rank == domain["tsw_bsw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::ne()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::se()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::nw()).rank == domain["tnw_bsw"]->rank);

  CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Edge::be()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["tsw_tne"]->rank);

  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Edge::bn()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Edge::bw()).rank == domain["tsw_bsw"]->rank);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Edge::be()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Edge::ne()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Edge::nw()).rank == domain["tsw_tnw"]->rank);

  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Edge::bs()).rank == domain["tsw_bsw"]->rank);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Edge::bn()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Edge::be()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Edge::ne()).rank == domain["tnw_tse"]->rank);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Edge::se()).rank == domain["tsw_tse"]->rank);

  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::bs()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::bn()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::bw()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::be()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::sw()).rank == domain["tsw_tsw"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::ne()).rank == domain["tne_tsw"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::se()).rank == domain["tse_tsw"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::nw()).rank == domain["tnw_tsw"]->rank);

  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::tn()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::bw()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::te()).rank == domain["tse_tse"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::be()).rank == domain["bse_tse"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::tw()).rank == domain["tsw_tse"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["tse_bne"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::nw()).rank == domain["tsw_bne"]->rank);

  CHECK(domain["tse_bse"]->getNormalNbrInfo(Edge::tn()).rank == domain["tse_tne"]->rank);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Edge::bn()).rank == domain["bse_tne"]->rank);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Edge::bw()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Edge::tw()).rank == domain["tse_tsw"]->rank);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Edge::nw()).rank == domain["tse_bnw"]->rank);

  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::bs()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::tn()).rank == domain["tne_tsw"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["tse_tsw"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::bw()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::te()).rank == domain["tse_tne"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::be()).rank == domain["bse_tne"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::tw()).rank == domain["tsw_tne"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::sw()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::ne()).rank == domain["tne_bse"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::se()).rank == domain["tse_bse"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::nw()).rank == domain["tnw_bse"]->rank);

  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::bs()).rank == domain["bse_tse"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::tn()).rank == domain["tne_tse"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::bn()).rank == domain["bne_tse"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::ts()).rank == domain["tse_tse"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::bw()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::tw()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::sw()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::nw()).rank == domain["tne_bsw"]->rank);

  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Edge::bw()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Edge::be()).rank == domain["tse_bse"]->rank);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["tse_tne"]->rank);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Edge::nw()).rank == domain["tsw_tne"]->rank);

  CHECK(domain["tse_tse"]->getNormalNbrInfo(Edge::bn()).rank == domain["tse_bne"]->rank);
  CHECK(domain["tse_tse"]->getNormalNbrInfo(Edge::bw()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["tse_tse"]->getNormalNbrInfo(Edge::nw()).rank == domain["tse_tnw"]->rank);

  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::bs()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::bn()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::bw()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::be()).rank == domain["tse_bne"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::sw()).rank == domain["tsw_tse"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::ne()).rank == domain["tne_tse"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::se()).rank == domain["tse_tse"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::nw()).rank == domain["tnw_tse"]->rank);

  CHECK(domain["tse_tne"]->getNormalNbrInfo(Edge::bs()).rank == domain["tse_bse"]->rank);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Edge::bn()).rank == domain["tne_bse"]->rank);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Edge::bw()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Edge::sw()).rank == domain["tse_tsw"]->rank);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Edge::nw()).rank == domain["tne_tsw"]->rank);

  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::bs()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::tn()).rank == domain["tnw_tnw"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bnw_tnw"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::ts()).rank == domain["tsw_tnw"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::te()).rank == domain["tnw_tse"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::be()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["tnw_bne"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::se()).rank == domain["tsw_bne"]->rank);

  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::bs()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::tn()).rank == domain["tnw_tne"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::bn()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::ts()).rank == domain["tsw_tne"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::bw()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::te()).rank == domain["tne_tsw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::be()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::tw()).rank == domain["tnw_tsw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::sw()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::ne()).rank == domain["tne_bnw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::se()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::nw()).rank == domain["tnw_bnw"]->rank);

  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Edge::bs()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["tnw_tsw"]->rank);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Edge::te()).rank == domain["tnw_tne"]->rank);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Edge::be()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Edge::se()).rank == domain["tnw_bse"]->rank);

  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::bs()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::ts()).rank == domain["tnw_tse"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::bw()).rank == domain["bnw_tnw"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::te()).rank == domain["tne_tnw"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::be()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::tw()).rank == domain["tnw_tnw"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::sw()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::se()).rank == domain["tne_bsw"]->rank);

  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Edge::bs()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["tnw_bnw"]->rank);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Edge::be()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["tnw_tne"]->rank);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Edge::se()).rank == domain["tsw_tne"]->rank);

  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::bs()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::bn()).rank == domain["tnw_bne"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::bw()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::be()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::sw()).rank == domain["tsw_tnw"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::ne()).rank == domain["tne_tnw"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::se()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::nw()).rank == domain["tnw_tnw"]->rank);

  CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Edge::bs()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Edge::be()).rank == domain["tnw_bne"]->rank);
  CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Edge::se()).rank == domain["tnw_tse"]->rank);

  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Edge::bs()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Edge::bw()).rank == domain["tnw_bnw"]->rank);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Edge::be()).rank == domain["tne_bnw"]->rank);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Edge::sw()).rank == domain["tnw_tsw"]->rank);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Edge::se()).rank == domain["tne_tsw"]->rank);

  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::bs()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::tn()).rank == domain["tne_tnw"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::ts()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::bw()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::te()).rank == domain["tne_tse"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::be()).rank == domain["bne_tse"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::tw()).rank == domain["tnw_tse"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::sw()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["tne_bne"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::se()).rank == domain["tse_bne"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::nw()).rank == domain["tnw_bne"]->rank);

  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::bs()).rank == domain["bse_tne"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::tn()).rank == domain["tne_tne"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::bn()).rank == domain["bne_tne"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::ts()).rank == domain["tse_tne"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::bw()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::tw()).rank == domain["tne_tsw"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::sw()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::nw()).rank == domain["tne_bnw"]->rank);

  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::bs()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::ts()).rank == domain["tne_tsw"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::bw()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::te()).rank == domain["tne_tne"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::be()).rank == domain["bne_tne"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::tw()).rank == domain["tnw_tne"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::sw()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::se()).rank == domain["tne_bse"]->rank);

  CHECK(domain["tne_bne"]->getNormalNbrInfo(Edge::bs()).rank == domain["bne_tse"]->rank);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Edge::ts()).rank == domain["tne_tse"]->rank);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Edge::bw()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Edge::tw()).rank == domain["tne_tnw"]->rank);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Edge::sw()).rank == domain["tne_bsw"]->rank);

  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::bs()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::bn()).rank == domain["tne_bnw"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::bw()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::be()).rank == domain["tne_bse"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::sw()).rank == domain["tsw_tne"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::ne()).rank == domain["tne_tne"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::se()).rank == domain["tse_tne"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::nw()).rank == domain["tnw_tne"]->rank);

  CHECK(domain["tne_tse"]->getNormalNbrInfo(Edge::bs()).rank == domain["tse_bne"]->rank);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Edge::bn()).rank == domain["tne_bne"]->rank);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Edge::bw()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Edge::sw()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Edge::nw()).rank == domain["tne_tnw"]->rank);

  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Edge::bs()).rank == domain["tne_bsw"]->rank);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Edge::bw()).rank == domain["tnw_bne"]->rank);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Edge::be()).rank == domain["tne_bne"]->rank);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Edge::sw()).rank == domain["tnw_tse"]->rank);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Edge::se()).rank == domain["tne_tse"]->rank);

  CHECK(domain["tne_tne"]->getNormalNbrInfo(Edge::bs()).rank == domain["tne_bse"]->rank);
  CHECK(domain["tne_tne"]->getNormalNbrInfo(Edge::bw()).rank == domain["tne_bnw"]->rank);
  CHECK(domain["tne_tne"]->getNormalNbrInfo(Edge::sw()).rank == domain["tne_tsw"]->rank);
}
