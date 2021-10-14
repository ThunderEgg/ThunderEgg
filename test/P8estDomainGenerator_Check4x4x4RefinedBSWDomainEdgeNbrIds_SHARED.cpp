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
Check4x4x4RefinedBSWDomainEdgeNeighborIds(const PatchVector& domain)
{
  CHECK(domain["bsw_bsw_bsw"]->getNormalNbrInfo(Edge::tn()).id == domain["bsw_bsw_tnw"]->id);
  CHECK(domain["bsw_bsw_bsw"]->getNormalNbrInfo(Edge::te()).id == domain["bsw_bsw_tse"]->id);
  CHECK(domain["bsw_bsw_bsw"]->getNormalNbrInfo(Edge::ne()).id == domain["bsw_bsw_bne"]->id);

  CHECK(domain["bsw_bsw_bse"]->getNormalNbrInfo(Edge::tn()).id == domain["bsw_bsw_tne"]->id);
  CHECK(domain["bsw_bsw_bse"]->getNormalNbrInfo(Edge::tw()).id == domain["bsw_bsw_tsw"]->id);
  CHECK(domain["bsw_bsw_bse"]->getNormalNbrInfo(Edge::nw()).id == domain["bsw_bsw_bnw"]->id);

  CHECK(domain["bsw_bsw_bnw"]->getNormalNbrInfo(Edge::ts()).id == domain["bsw_bsw_tsw"]->id);
  CHECK(domain["bsw_bsw_bnw"]->getNormalNbrInfo(Edge::te()).id == domain["bsw_bsw_tne"]->id);
  CHECK(domain["bsw_bsw_bnw"]->getNormalNbrInfo(Edge::se()).id == domain["bsw_bsw_bse"]->id);

  CHECK(domain["bsw_bsw_bne"]->getNormalNbrInfo(Edge::ts()).id == domain["bsw_bsw_tse"]->id);
  CHECK(domain["bsw_bsw_bne"]->getNormalNbrInfo(Edge::tw()).id == domain["bsw_bsw_tnw"]->id);
  CHECK(domain["bsw_bsw_bne"]->getNormalNbrInfo(Edge::sw()).id == domain["bsw_bsw_bsw"]->id);
  CHECK(domain["bsw_bsw_bne"]->getCoarseNbrInfo(Edge::ne()).id == domain["bsw_bne"]->id);

  CHECK(domain["bsw_bsw_tsw"]->getNormalNbrInfo(Edge::bn()).id == domain["bsw_bsw_bnw"]->id);
  CHECK(domain["bsw_bsw_tsw"]->getNormalNbrInfo(Edge::be()).id == domain["bsw_bsw_bse"]->id);
  CHECK(domain["bsw_bsw_tsw"]->getNormalNbrInfo(Edge::ne()).id == domain["bsw_bsw_tne"]->id);

  CHECK(domain["bsw_bsw_tse"]->getNormalNbrInfo(Edge::bn()).id == domain["bsw_bsw_bne"]->id);
  CHECK(domain["bsw_bsw_tse"]->getNormalNbrInfo(Edge::bw()).id == domain["bsw_bsw_bsw"]->id);
  CHECK(domain["bsw_bsw_tse"]->getCoarseNbrInfo(Edge::te()).id == domain["bsw_tse"]->id);
  CHECK(domain["bsw_bsw_tse"]->getNormalNbrInfo(Edge::nw()).id == domain["bsw_bsw_tnw"]->id);

  CHECK(domain["bsw_bsw_tnw"]->getNormalNbrInfo(Edge::bs()).id == domain["bsw_bsw_bsw"]->id);
  CHECK(domain["bsw_bsw_tnw"]->getCoarseNbrInfo(Edge::tn()).id == domain["bsw_tnw"]->id);
  CHECK(domain["bsw_bsw_tnw"]->getNormalNbrInfo(Edge::be()).id == domain["bsw_bsw_bne"]->id);
  CHECK(domain["bsw_bsw_tnw"]->getNormalNbrInfo(Edge::se()).id == domain["bsw_bsw_tse"]->id);

  CHECK(domain["bsw_bsw_tne"]->getNormalNbrInfo(Edge::bs()).id == domain["bsw_bsw_bse"]->id);
  CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Edge::tn()).id == domain["bsw_tnw"]->id);
  CHECK(domain["bsw_bsw_tne"]->getNormalNbrInfo(Edge::bw()).id == domain["bsw_bsw_bnw"]->id);
  CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Edge::te()).id == domain["bsw_tse"]->id);
  CHECK(domain["bsw_bsw_tne"]->getNormalNbrInfo(Edge::sw()).id == domain["bsw_bsw_tsw"]->id);
  CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Edge::ne()).id == domain["bsw_bne"]->id);

  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Edge::tn()).id == domain["bsw_tne"]->id);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Edge::te()).id == domain["bse_tsw"]->id);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Edge::tw()).id == domain["bsw_tsw"]->id);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Edge::ne()).id == domain["bse_bnw"]->id);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Edge::nw()).id == domain["bsw_bnw"]->id);

  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Edge::tn()).id == domain["bnw_tsw"]->id);
  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Edge::ts()).id == domain["bsw_tsw"]->id);
  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Edge::te()).id == domain["bsw_tne"]->id);
  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Edge::ne()).id == domain["bnw_bse"]->id);
  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Edge::se()).id == domain["bsw_bse"]->id);

  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::tn()).id == domain["bnw_tse"]->id);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::ts()).id == domain["bsw_tse"]->id);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::te()).id == domain["bse_tnw"]->id);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::tw()).id == domain["bsw_tnw"]->id);
  CHECK(domain["bsw_bne"]->getFineNbrInfo(Edge::sw()).ids[0] == domain["bsw_bsw_bne"]->id);
  CHECK(domain["bsw_bne"]->getFineNbrInfo(Edge::sw()).ids[1] == domain["bsw_bsw_tne"]->id);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::ne()).id == domain["bne_bsw"]->id);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::se()).id == domain["bse_bsw"]->id);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Edge::nw()).id == domain["bnw_bsw"]->id);

  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Edge::tn()).id == domain["tsw_bnw"]->id);
  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Edge::bn()).id == domain["bsw_bnw"]->id);
  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Edge::te()).id == domain["tsw_bse"]->id);
  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Edge::be()).id == domain["bsw_bse"]->id);
  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Edge::ne()).id == domain["bsw_tne"]->id);

  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::tn()).id == domain["tsw_bne"]->id);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::bn()).id == domain["bsw_bne"]->id);
  CHECK(domain["bsw_tse"]->getFineNbrInfo(Edge::bw()).ids[0] == domain["bsw_bsw_tse"]->id);
  CHECK(domain["bsw_tse"]->getFineNbrInfo(Edge::bw()).ids[1] == domain["bsw_bsw_tne"]->id);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::te()).id == domain["tse_bsw"]->id);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::be()).id == domain["bse_bsw"]->id);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::tw()).id == domain["tsw_bsw"]->id);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::ne()).id == domain["bse_tnw"]->id);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Edge::nw()).id == domain["bsw_tnw"]->id);

  CHECK(domain["bsw_tnw"]->getFineNbrInfo(Edge::bs()).ids[0] == domain["bsw_bsw_tnw"]->id);
  CHECK(domain["bsw_tnw"]->getFineNbrInfo(Edge::bs()).ids[1] == domain["bsw_bsw_tne"]->id);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::tn()).id == domain["tnw_bsw"]->id);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::bn()).id == domain["bnw_bsw"]->id);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::ts()).id == domain["tsw_bsw"]->id);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::te()).id == domain["tsw_bne"]->id);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::be()).id == domain["bsw_bne"]->id);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::ne()).id == domain["bnw_tse"]->id);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Edge::se()).id == domain["bsw_tse"]->id);

  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::bs()).id == domain["bsw_bse"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::tn()).id == domain["tnw_bse"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::bn()).id == domain["bnw_bse"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::ts()).id == domain["tsw_bse"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::bw()).id == domain["bsw_bnw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::te()).id == domain["tse_bnw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::be()).id == domain["bse_bnw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::tw()).id == domain["tsw_bnw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::sw()).id == domain["bsw_tsw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::ne()).id == domain["bne_tsw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::se()).id == domain["bse_tsw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Edge::nw()).id == domain["bnw_tsw"]->id);

  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Edge::tn()).id == domain["bse_tnw"]->id);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Edge::te()).id == domain["bse_tse"]->id);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Edge::tw()).id == domain["bsw_tse"]->id);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Edge::ne()).id == domain["bse_bne"]->id);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Edge::nw()).id == domain["bsw_bne"]->id);

  CHECK(domain["bse_bse"]->getNormalNbrInfo(Edge::tn()).id == domain["bse_tne"]->id);
  CHECK(domain["bse_bse"]->getNormalNbrInfo(Edge::tw()).id == domain["bse_tsw"]->id);
  CHECK(domain["bse_bse"]->getNormalNbrInfo(Edge::nw()).id == domain["bse_bnw"]->id);

  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::tn()).id == domain["bne_tsw"]->id);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::ts()).id == domain["bse_tsw"]->id);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::te()).id == domain["bse_tne"]->id);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::tw()).id == domain["bsw_tne"]->id);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::sw()).id == domain["bsw_bse"]->id);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::ne()).id == domain["bne_bse"]->id);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::se()).id == domain["bse_bse"]->id);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Edge::nw()).id == domain["bnw_bse"]->id);

  CHECK(domain["bse_bne"]->getNormalNbrInfo(Edge::tn()).id == domain["bne_tse"]->id);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Edge::ts()).id == domain["bse_tse"]->id);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Edge::tw()).id == domain["bse_tnw"]->id);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Edge::sw()).id == domain["bse_bsw"]->id);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Edge::nw()).id == domain["bne_bsw"]->id);

  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::tn()).id == domain["tse_bnw"]->id);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::bn()).id == domain["bse_bnw"]->id);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::bw()).id == domain["bsw_bse"]->id);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::te()).id == domain["tse_bse"]->id);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::be()).id == domain["bse_bse"]->id);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::tw()).id == domain["tsw_bse"]->id);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::ne()).id == domain["bse_tne"]->id);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Edge::nw()).id == domain["bsw_tne"]->id);

  CHECK(domain["bse_tse"]->getNormalNbrInfo(Edge::tn()).id == domain["tse_bne"]->id);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Edge::bn()).id == domain["bse_bne"]->id);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Edge::bw()).id == domain["bse_bsw"]->id);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Edge::tw()).id == domain["tse_bsw"]->id);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Edge::nw()).id == domain["bse_tnw"]->id);

  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::bs()).id == domain["bse_bsw"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::tn()).id == domain["tne_bsw"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::bn()).id == domain["bne_bsw"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::ts()).id == domain["tse_bsw"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::bw()).id == domain["bsw_bne"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::te()).id == domain["tse_bne"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::be()).id == domain["bse_bne"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::tw()).id == domain["tsw_bne"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::sw()).id == domain["bsw_tse"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::ne()).id == domain["bne_tse"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::se()).id == domain["bse_tse"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Edge::nw()).id == domain["bnw_tse"]->id);

  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::bs()).id == domain["bse_bse"]->id);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::tn()).id == domain["tne_bse"]->id);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::bn()).id == domain["bne_bse"]->id);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::ts()).id == domain["tse_bse"]->id);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::bw()).id == domain["bse_bnw"]->id);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::tw()).id == domain["tse_bnw"]->id);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::sw()).id == domain["bse_tsw"]->id);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Edge::nw()).id == domain["bne_tsw"]->id);

  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Edge::tn()).id == domain["bnw_tnw"]->id);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Edge::ts()).id == domain["bsw_tnw"]->id);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Edge::te()).id == domain["bnw_tse"]->id);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Edge::ne()).id == domain["bnw_bne"]->id);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Edge::se()).id == domain["bsw_bne"]->id);

  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::tn()).id == domain["bnw_tne"]->id);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::ts()).id == domain["bsw_tne"]->id);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::te()).id == domain["bne_tsw"]->id);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::tw()).id == domain["bnw_tsw"]->id);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::sw()).id == domain["bsw_bnw"]->id);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::ne()).id == domain["bne_bnw"]->id);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::se()).id == domain["bse_bnw"]->id);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Edge::nw()).id == domain["bnw_bnw"]->id);

  CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Edge::ts()).id == domain["bnw_tsw"]->id);
  CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Edge::te()).id == domain["bnw_tne"]->id);
  CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Edge::se()).id == domain["bnw_bse"]->id);

  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Edge::ts()).id == domain["bnw_tse"]->id);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Edge::te()).id == domain["bne_tnw"]->id);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Edge::tw()).id == domain["bnw_tnw"]->id);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Edge::sw()).id == domain["bnw_bsw"]->id);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Edge::se()).id == domain["bne_bsw"]->id);

  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::bs()).id == domain["bsw_bnw"]->id);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::tn()).id == domain["tnw_bnw"]->id);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::bn()).id == domain["bnw_bnw"]->id);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::ts()).id == domain["tsw_bnw"]->id);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::te()).id == domain["tnw_bse"]->id);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::be()).id == domain["bnw_bse"]->id);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::ne()).id == domain["bnw_tne"]->id);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Edge::se()).id == domain["bsw_tne"]->id);

  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::bs()).id == domain["bsw_bne"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::tn()).id == domain["tnw_bne"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::bn()).id == domain["bnw_bne"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::ts()).id == domain["tsw_bne"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::bw()).id == domain["bnw_bsw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::te()).id == domain["tne_bsw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::be()).id == domain["bne_bsw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::tw()).id == domain["tnw_bsw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::sw()).id == domain["bsw_tnw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::ne()).id == domain["bne_tnw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::se()).id == domain["bse_tnw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Edge::nw()).id == domain["bnw_tnw"]->id);

  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Edge::bs()).id == domain["bnw_bsw"]->id);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Edge::ts()).id == domain["tnw_bsw"]->id);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Edge::te()).id == domain["tnw_bne"]->id);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Edge::be()).id == domain["bnw_bne"]->id);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Edge::se()).id == domain["bnw_tse"]->id);

  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::bs()).id == domain["bnw_bse"]->id);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::ts()).id == domain["tnw_bse"]->id);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::bw()).id == domain["bnw_bnw"]->id);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::te()).id == domain["tne_bnw"]->id);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::be()).id == domain["bne_bnw"]->id);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::tw()).id == domain["tnw_bnw"]->id);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::sw()).id == domain["bnw_tsw"]->id);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Edge::se()).id == domain["bne_tsw"]->id);

  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::tn()).id == domain["bne_tnw"]->id);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::ts()).id == domain["bse_tnw"]->id);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::te()).id == domain["bne_tse"]->id);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::tw()).id == domain["bnw_tse"]->id);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::sw()).id == domain["bsw_bne"]->id);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::ne()).id == domain["bne_bne"]->id);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::se()).id == domain["bse_bne"]->id);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Edge::nw()).id == domain["bnw_bne"]->id);

  CHECK(domain["bne_bse"]->getNormalNbrInfo(Edge::tn()).id == domain["bne_tne"]->id);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Edge::ts()).id == domain["bse_tne"]->id);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Edge::tw()).id == domain["bne_tsw"]->id);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Edge::sw()).id == domain["bse_bnw"]->id);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Edge::nw()).id == domain["bne_bnw"]->id);

  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Edge::ts()).id == domain["bne_tsw"]->id);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Edge::te()).id == domain["bne_tne"]->id);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Edge::tw()).id == domain["bnw_tne"]->id);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Edge::sw()).id == domain["bnw_bse"]->id);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Edge::se()).id == domain["bne_bse"]->id);

  CHECK(domain["bne_bne"]->getNormalNbrInfo(Edge::ts()).id == domain["bne_tse"]->id);
  CHECK(domain["bne_bne"]->getNormalNbrInfo(Edge::tw()).id == domain["bne_tnw"]->id);
  CHECK(domain["bne_bne"]->getNormalNbrInfo(Edge::sw()).id == domain["bne_bsw"]->id);

  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::bs()).id == domain["bse_bnw"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::tn()).id == domain["tne_bnw"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::bn()).id == domain["bne_bnw"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::ts()).id == domain["tse_bnw"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::bw()).id == domain["bnw_bse"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::te()).id == domain["tne_bse"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::be()).id == domain["bne_bse"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::tw()).id == domain["tnw_bse"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::sw()).id == domain["bsw_tne"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::ne()).id == domain["bne_tne"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::se()).id == domain["bse_tne"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Edge::nw()).id == domain["bnw_tne"]->id);

  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::bs()).id == domain["bse_bne"]->id);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::tn()).id == domain["tne_bne"]->id);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::bn()).id == domain["bne_bne"]->id);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::ts()).id == domain["tse_bne"]->id);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::bw()).id == domain["bne_bsw"]->id);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::tw()).id == domain["tne_bsw"]->id);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::sw()).id == domain["bse_tnw"]->id);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Edge::nw()).id == domain["bne_tnw"]->id);

  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::bs()).id == domain["bne_bsw"]->id);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::ts()).id == domain["tne_bsw"]->id);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::bw()).id == domain["bnw_bne"]->id);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::te()).id == domain["tne_bne"]->id);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::be()).id == domain["bne_bne"]->id);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::tw()).id == domain["tnw_bne"]->id);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::sw()).id == domain["bnw_tse"]->id);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Edge::se()).id == domain["bne_tse"]->id);

  CHECK(domain["bne_tne"]->getNormalNbrInfo(Edge::bs()).id == domain["bne_bse"]->id);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Edge::ts()).id == domain["tne_bse"]->id);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Edge::bw()).id == domain["bne_bnw"]->id);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Edge::tw()).id == domain["tne_bnw"]->id);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Edge::sw()).id == domain["bne_tsw"]->id);

  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Edge::tn()).id == domain["tsw_tnw"]->id);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Edge::bn()).id == domain["bsw_tnw"]->id);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Edge::te()).id == domain["tsw_tse"]->id);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Edge::be()).id == domain["bsw_tse"]->id);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Edge::ne()).id == domain["tsw_bne"]->id);

  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::tn()).id == domain["tsw_tne"]->id);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::bn()).id == domain["bsw_tne"]->id);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::bw()).id == domain["bsw_tsw"]->id);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::te()).id == domain["tse_tsw"]->id);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::be()).id == domain["bse_tsw"]->id);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::tw()).id == domain["tsw_tsw"]->id);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::ne()).id == domain["tse_bnw"]->id);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Edge::nw()).id == domain["tsw_bnw"]->id);

  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::bs()).id == domain["bsw_tsw"]->id);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::tn()).id == domain["tnw_tsw"]->id);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::bn()).id == domain["bnw_tsw"]->id);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::ts()).id == domain["tsw_tsw"]->id);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::te()).id == domain["tsw_tne"]->id);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::be()).id == domain["bsw_tne"]->id);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::ne()).id == domain["tnw_bse"]->id);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Edge::se()).id == domain["tsw_bse"]->id);

  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::bs()).id == domain["bsw_tse"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::tn()).id == domain["tnw_tse"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::bn()).id == domain["bnw_tse"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::ts()).id == domain["tsw_tse"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::bw()).id == domain["bsw_tnw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::te()).id == domain["tse_tnw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::be()).id == domain["bse_tnw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::tw()).id == domain["tsw_tnw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::sw()).id == domain["tsw_bsw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::ne()).id == domain["tne_bsw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::se()).id == domain["tse_bsw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Edge::nw()).id == domain["tnw_bsw"]->id);

  CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Edge::bn()).id == domain["tsw_bnw"]->id);
  CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Edge::be()).id == domain["tsw_bse"]->id);
  CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Edge::ne()).id == domain["tsw_tne"]->id);

  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Edge::bn()).id == domain["tsw_bne"]->id);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Edge::bw()).id == domain["tsw_bsw"]->id);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Edge::be()).id == domain["tse_bsw"]->id);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Edge::ne()).id == domain["tse_tnw"]->id);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Edge::nw()).id == domain["tsw_tnw"]->id);

  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Edge::bs()).id == domain["tsw_bsw"]->id);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Edge::bn()).id == domain["tnw_bsw"]->id);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Edge::be()).id == domain["tsw_bne"]->id);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Edge::ne()).id == domain["tnw_tse"]->id);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Edge::se()).id == domain["tsw_tse"]->id);

  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::bs()).id == domain["tsw_bse"]->id);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::bn()).id == domain["tnw_bse"]->id);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::bw()).id == domain["tsw_bnw"]->id);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::be()).id == domain["tse_bnw"]->id);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::sw()).id == domain["tsw_tsw"]->id);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::ne()).id == domain["tne_tsw"]->id);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::se()).id == domain["tse_tsw"]->id);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Edge::nw()).id == domain["tnw_tsw"]->id);

  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::tn()).id == domain["tse_tnw"]->id);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::bn()).id == domain["bse_tnw"]->id);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::bw()).id == domain["bsw_tse"]->id);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::te()).id == domain["tse_tse"]->id);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::be()).id == domain["bse_tse"]->id);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::tw()).id == domain["tsw_tse"]->id);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::ne()).id == domain["tse_bne"]->id);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Edge::nw()).id == domain["tsw_bne"]->id);

  CHECK(domain["tse_bse"]->getNormalNbrInfo(Edge::tn()).id == domain["tse_tne"]->id);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Edge::bn()).id == domain["bse_tne"]->id);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Edge::bw()).id == domain["bse_tsw"]->id);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Edge::tw()).id == domain["tse_tsw"]->id);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Edge::nw()).id == domain["tse_bnw"]->id);

  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::bs()).id == domain["bse_tsw"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::tn()).id == domain["tne_tsw"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::bn()).id == domain["bne_tsw"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::ts()).id == domain["tse_tsw"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::bw()).id == domain["bsw_tne"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::te()).id == domain["tse_tne"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::be()).id == domain["bse_tne"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::tw()).id == domain["tsw_tne"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::sw()).id == domain["tsw_bse"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::ne()).id == domain["tne_bse"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::se()).id == domain["tse_bse"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Edge::nw()).id == domain["tnw_bse"]->id);

  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::bs()).id == domain["bse_tse"]->id);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::tn()).id == domain["tne_tse"]->id);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::bn()).id == domain["bne_tse"]->id);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::ts()).id == domain["tse_tse"]->id);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::bw()).id == domain["bse_tnw"]->id);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::tw()).id == domain["tse_tnw"]->id);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::sw()).id == domain["tse_bsw"]->id);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Edge::nw()).id == domain["tne_bsw"]->id);

  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Edge::bn()).id == domain["tse_bnw"]->id);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Edge::bw()).id == domain["tsw_bse"]->id);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Edge::be()).id == domain["tse_bse"]->id);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Edge::ne()).id == domain["tse_tne"]->id);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Edge::nw()).id == domain["tsw_tne"]->id);

  CHECK(domain["tse_tse"]->getNormalNbrInfo(Edge::bn()).id == domain["tse_bne"]->id);
  CHECK(domain["tse_tse"]->getNormalNbrInfo(Edge::bw()).id == domain["tse_bsw"]->id);
  CHECK(domain["tse_tse"]->getNormalNbrInfo(Edge::nw()).id == domain["tse_tnw"]->id);

  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::bs()).id == domain["tse_bsw"]->id);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::bn()).id == domain["tne_bsw"]->id);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::bw()).id == domain["tsw_bne"]->id);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::be()).id == domain["tse_bne"]->id);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::sw()).id == domain["tsw_tse"]->id);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::ne()).id == domain["tne_tse"]->id);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::se()).id == domain["tse_tse"]->id);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Edge::nw()).id == domain["tnw_tse"]->id);

  CHECK(domain["tse_tne"]->getNormalNbrInfo(Edge::bs()).id == domain["tse_bse"]->id);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Edge::bn()).id == domain["tne_bse"]->id);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Edge::bw()).id == domain["tse_bnw"]->id);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Edge::sw()).id == domain["tse_tsw"]->id);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Edge::nw()).id == domain["tne_tsw"]->id);

  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::bs()).id == domain["bsw_tnw"]->id);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::tn()).id == domain["tnw_tnw"]->id);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::bn()).id == domain["bnw_tnw"]->id);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::ts()).id == domain["tsw_tnw"]->id);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::te()).id == domain["tnw_tse"]->id);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::be()).id == domain["bnw_tse"]->id);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::ne()).id == domain["tnw_bne"]->id);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Edge::se()).id == domain["tsw_bne"]->id);

  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::bs()).id == domain["bsw_tne"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::tn()).id == domain["tnw_tne"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::bn()).id == domain["bnw_tne"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::ts()).id == domain["tsw_tne"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::bw()).id == domain["bnw_tsw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::te()).id == domain["tne_tsw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::be()).id == domain["bne_tsw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::tw()).id == domain["tnw_tsw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::sw()).id == domain["tsw_bnw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::ne()).id == domain["tne_bnw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::se()).id == domain["tse_bnw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Edge::nw()).id == domain["tnw_bnw"]->id);

  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Edge::bs()).id == domain["bnw_tsw"]->id);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Edge::ts()).id == domain["tnw_tsw"]->id);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Edge::te()).id == domain["tnw_tne"]->id);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Edge::be()).id == domain["bnw_tne"]->id);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Edge::se()).id == domain["tnw_bse"]->id);

  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::bs()).id == domain["bnw_tse"]->id);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::ts()).id == domain["tnw_tse"]->id);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::bw()).id == domain["bnw_tnw"]->id);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::te()).id == domain["tne_tnw"]->id);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::be()).id == domain["bne_tnw"]->id);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::tw()).id == domain["tnw_tnw"]->id);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::sw()).id == domain["tnw_bsw"]->id);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Edge::se()).id == domain["tne_bsw"]->id);

  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Edge::bs()).id == domain["tsw_bnw"]->id);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Edge::bn()).id == domain["tnw_bnw"]->id);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Edge::be()).id == domain["tnw_bse"]->id);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Edge::ne()).id == domain["tnw_tne"]->id);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Edge::se()).id == domain["tsw_tne"]->id);

  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::bs()).id == domain["tsw_bne"]->id);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::bn()).id == domain["tnw_bne"]->id);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::bw()).id == domain["tnw_bsw"]->id);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::be()).id == domain["tne_bsw"]->id);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::sw()).id == domain["tsw_tnw"]->id);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::ne()).id == domain["tne_tnw"]->id);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::se()).id == domain["tse_tnw"]->id);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Edge::nw()).id == domain["tnw_tnw"]->id);

  CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Edge::bs()).id == domain["tnw_bsw"]->id);
  CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Edge::be()).id == domain["tnw_bne"]->id);
  CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Edge::se()).id == domain["tnw_tse"]->id);

  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Edge::bs()).id == domain["tnw_bse"]->id);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Edge::bw()).id == domain["tnw_bnw"]->id);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Edge::be()).id == domain["tne_bnw"]->id);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Edge::sw()).id == domain["tnw_tsw"]->id);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Edge::se()).id == domain["tne_tsw"]->id);

  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::bs()).id == domain["bse_tnw"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::tn()).id == domain["tne_tnw"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::bn()).id == domain["bne_tnw"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::ts()).id == domain["tse_tnw"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::bw()).id == domain["bnw_tse"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::te()).id == domain["tne_tse"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::be()).id == domain["bne_tse"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::tw()).id == domain["tnw_tse"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::sw()).id == domain["tsw_bne"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::ne()).id == domain["tne_bne"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::se()).id == domain["tse_bne"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Edge::nw()).id == domain["tnw_bne"]->id);

  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::bs()).id == domain["bse_tne"]->id);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::tn()).id == domain["tne_tne"]->id);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::bn()).id == domain["bne_tne"]->id);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::ts()).id == domain["tse_tne"]->id);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::bw()).id == domain["bne_tsw"]->id);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::tw()).id == domain["tne_tsw"]->id);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::sw()).id == domain["tse_bnw"]->id);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Edge::nw()).id == domain["tne_bnw"]->id);

  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::bs()).id == domain["bne_tsw"]->id);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::ts()).id == domain["tne_tsw"]->id);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::bw()).id == domain["bnw_tne"]->id);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::te()).id == domain["tne_tne"]->id);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::be()).id == domain["bne_tne"]->id);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::tw()).id == domain["tnw_tne"]->id);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::sw()).id == domain["tnw_bse"]->id);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Edge::se()).id == domain["tne_bse"]->id);

  CHECK(domain["tne_bne"]->getNormalNbrInfo(Edge::bs()).id == domain["bne_tse"]->id);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Edge::ts()).id == domain["tne_tse"]->id);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Edge::bw()).id == domain["bne_tnw"]->id);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Edge::tw()).id == domain["tne_tnw"]->id);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Edge::sw()).id == domain["tne_bsw"]->id);

  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::bs()).id == domain["tse_bnw"]->id);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::bn()).id == domain["tne_bnw"]->id);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::bw()).id == domain["tnw_bse"]->id);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::be()).id == domain["tne_bse"]->id);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::sw()).id == domain["tsw_tne"]->id);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::ne()).id == domain["tne_tne"]->id);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::se()).id == domain["tse_tne"]->id);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Edge::nw()).id == domain["tnw_tne"]->id);

  CHECK(domain["tne_tse"]->getNormalNbrInfo(Edge::bs()).id == domain["tse_bne"]->id);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Edge::bn()).id == domain["tne_bne"]->id);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Edge::bw()).id == domain["tne_bsw"]->id);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Edge::sw()).id == domain["tse_tnw"]->id);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Edge::nw()).id == domain["tne_tnw"]->id);

  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Edge::bs()).id == domain["tne_bsw"]->id);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Edge::bw()).id == domain["tnw_bne"]->id);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Edge::be()).id == domain["tne_bne"]->id);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Edge::sw()).id == domain["tnw_tse"]->id);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Edge::se()).id == domain["tne_tse"]->id);

  CHECK(domain["tne_tne"]->getNormalNbrInfo(Edge::bs()).id == domain["tne_bse"]->id);
  CHECK(domain["tne_tne"]->getNormalNbrInfo(Edge::bw()).id == domain["tne_bnw"]->id);
  CHECK(domain["tne_tne"]->getNormalNbrInfo(Edge::sw()).id == domain["tne_tsw"]->id);
}
