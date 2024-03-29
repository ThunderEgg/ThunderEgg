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
Check4x4x4DomainEdgeNeighborIds(const PatchVector& domain)
{
  // edge nbr id
  CHECK_EQ(domain["bsw_bsw"]->getNormalNbrInfo(Edge::tn()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["bsw_bsw"]->getNormalNbrInfo(Edge::te()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["bsw_bsw"]->getNormalNbrInfo(Edge::ne()).id, domain["bsw_bne"]->id);

  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Edge::tn()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Edge::te()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Edge::tw()).id, domain["bsw_tsw"]->id);
  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Edge::ne()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Edge::nw()).id, domain["bsw_bnw"]->id);

  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Edge::tn()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Edge::ts()).id, domain["bsw_tsw"]->id);
  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Edge::te()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Edge::ne()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Edge::se()).id, domain["bsw_bse"]->id);

  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Edge::tn()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Edge::ts()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Edge::te()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Edge::tw()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Edge::sw()).id, domain["bsw_bsw"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Edge::ne()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Edge::se()).id, domain["bse_bsw"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Edge::nw()).id, domain["bnw_bsw"]->id);

  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Edge::tn()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Edge::bn()).id, domain["bsw_bnw"]->id);
  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Edge::te()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Edge::be()).id, domain["bsw_bse"]->id);
  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Edge::ne()).id, domain["bsw_tne"]->id);

  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Edge::tn()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Edge::bn()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Edge::bw()).id, domain["bsw_bsw"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Edge::te()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Edge::be()).id, domain["bse_bsw"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Edge::tw()).id, domain["tsw_bsw"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Edge::ne()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Edge::nw()).id, domain["bsw_tnw"]->id);

  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Edge::bs()).id, domain["bsw_bsw"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Edge::tn()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Edge::bn()).id, domain["bnw_bsw"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Edge::ts()).id, domain["tsw_bsw"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Edge::te()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Edge::be()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Edge::ne()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Edge::se()).id, domain["bsw_tse"]->id);

  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Edge::bs()).id, domain["bsw_bse"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Edge::tn()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Edge::bn()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Edge::ts()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Edge::bw()).id, domain["bsw_bnw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Edge::te()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Edge::be()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Edge::tw()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Edge::sw()).id, domain["bsw_tsw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Edge::ne()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Edge::se()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Edge::nw()).id, domain["bnw_tsw"]->id);

  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Edge::tn()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Edge::te()).id, domain["bse_tse"]->id);
  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Edge::tw()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Edge::ne()).id, domain["bse_bne"]->id);
  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Edge::nw()).id, domain["bsw_bne"]->id);

  CHECK_EQ(domain["bse_bse"]->getNormalNbrInfo(Edge::tn()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["bse_bse"]->getNormalNbrInfo(Edge::tw()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["bse_bse"]->getNormalNbrInfo(Edge::nw()).id, domain["bse_bnw"]->id);

  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Edge::tn()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Edge::ts()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Edge::te()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Edge::tw()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Edge::sw()).id, domain["bsw_bse"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Edge::ne()).id, domain["bne_bse"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Edge::se()).id, domain["bse_bse"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Edge::nw()).id, domain["bnw_bse"]->id);

  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Edge::tn()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Edge::ts()).id, domain["bse_tse"]->id);
  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Edge::tw()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Edge::sw()).id, domain["bse_bsw"]->id);
  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Edge::nw()).id, domain["bne_bsw"]->id);

  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Edge::tn()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Edge::bn()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Edge::bw()).id, domain["bsw_bse"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Edge::te()).id, domain["tse_bse"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Edge::be()).id, domain["bse_bse"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Edge::tw()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Edge::ne()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Edge::nw()).id, domain["bsw_tne"]->id);

  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Edge::tn()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Edge::bn()).id, domain["bse_bne"]->id);
  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Edge::bw()).id, domain["bse_bsw"]->id);
  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Edge::tw()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Edge::nw()).id, domain["bse_tnw"]->id);

  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Edge::bs()).id, domain["bse_bsw"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Edge::tn()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Edge::bn()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Edge::ts()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Edge::bw()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Edge::te()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Edge::be()).id, domain["bse_bne"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Edge::tw()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Edge::sw()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Edge::ne()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Edge::se()).id, domain["bse_tse"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Edge::nw()).id, domain["bnw_tse"]->id);

  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Edge::bs()).id, domain["bse_bse"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Edge::tn()).id, domain["tne_bse"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Edge::bn()).id, domain["bne_bse"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Edge::ts()).id, domain["tse_bse"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Edge::bw()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Edge::tw()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Edge::sw()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Edge::nw()).id, domain["bne_tsw"]->id);

  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Edge::tn()).id, domain["bnw_tnw"]->id);
  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Edge::ts()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Edge::te()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Edge::ne()).id, domain["bnw_bne"]->id);
  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Edge::se()).id, domain["bsw_bne"]->id);

  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Edge::tn()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Edge::ts()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Edge::te()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Edge::tw()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Edge::sw()).id, domain["bsw_bnw"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Edge::ne()).id, domain["bne_bnw"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Edge::se()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Edge::nw()).id, domain["bnw_bnw"]->id);

  CHECK_EQ(domain["bnw_bnw"]->getNormalNbrInfo(Edge::ts()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["bnw_bnw"]->getNormalNbrInfo(Edge::te()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["bnw_bnw"]->getNormalNbrInfo(Edge::se()).id, domain["bnw_bse"]->id);

  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Edge::ts()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Edge::te()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Edge::tw()).id, domain["bnw_tnw"]->id);
  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Edge::sw()).id, domain["bnw_bsw"]->id);
  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Edge::se()).id, domain["bne_bsw"]->id);

  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Edge::bs()).id, domain["bsw_bnw"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Edge::tn()).id, domain["tnw_bnw"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Edge::bn()).id, domain["bnw_bnw"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Edge::ts()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Edge::te()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Edge::be()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Edge::ne()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Edge::se()).id, domain["bsw_tne"]->id);

  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Edge::bs()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Edge::tn()).id, domain["tnw_bne"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Edge::bn()).id, domain["bnw_bne"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Edge::ts()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Edge::bw()).id, domain["bnw_bsw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Edge::te()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Edge::be()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Edge::tw()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Edge::sw()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Edge::ne()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Edge::se()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Edge::nw()).id, domain["bnw_tnw"]->id);

  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Edge::bs()).id, domain["bnw_bsw"]->id);
  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Edge::ts()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Edge::te()).id, domain["tnw_bne"]->id);
  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Edge::be()).id, domain["bnw_bne"]->id);
  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Edge::se()).id, domain["bnw_tse"]->id);

  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Edge::bs()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Edge::ts()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Edge::bw()).id, domain["bnw_bnw"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Edge::te()).id, domain["tne_bnw"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Edge::be()).id, domain["bne_bnw"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Edge::tw()).id, domain["tnw_bnw"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Edge::sw()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Edge::se()).id, domain["bne_tsw"]->id);

  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Edge::tn()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Edge::ts()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Edge::te()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Edge::tw()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Edge::sw()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Edge::ne()).id, domain["bne_bne"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Edge::se()).id, domain["bse_bne"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Edge::nw()).id, domain["bnw_bne"]->id);

  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Edge::tn()).id, domain["bne_tne"]->id);
  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Edge::ts()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Edge::tw()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Edge::sw()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Edge::nw()).id, domain["bne_bnw"]->id);

  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Edge::ts()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Edge::te()).id, domain["bne_tne"]->id);
  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Edge::tw()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Edge::sw()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Edge::se()).id, domain["bne_bse"]->id);

  CHECK_EQ(domain["bne_bne"]->getNormalNbrInfo(Edge::ts()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["bne_bne"]->getNormalNbrInfo(Edge::tw()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["bne_bne"]->getNormalNbrInfo(Edge::sw()).id, domain["bne_bsw"]->id);

  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Edge::bs()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Edge::tn()).id, domain["tne_bnw"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Edge::bn()).id, domain["bne_bnw"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Edge::ts()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Edge::bw()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Edge::te()).id, domain["tne_bse"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Edge::be()).id, domain["bne_bse"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Edge::tw()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Edge::sw()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Edge::ne()).id, domain["bne_tne"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Edge::se()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Edge::nw()).id, domain["bnw_tne"]->id);

  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Edge::bs()).id, domain["bse_bne"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Edge::tn()).id, domain["tne_bne"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Edge::bn()).id, domain["bne_bne"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Edge::ts()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Edge::bw()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Edge::tw()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Edge::sw()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Edge::nw()).id, domain["bne_tnw"]->id);

  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Edge::bs()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Edge::ts()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Edge::bw()).id, domain["bnw_bne"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Edge::te()).id, domain["tne_bne"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Edge::be()).id, domain["bne_bne"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Edge::tw()).id, domain["tnw_bne"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Edge::sw()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Edge::se()).id, domain["bne_tse"]->id);

  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Edge::bs()).id, domain["bne_bse"]->id);
  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Edge::ts()).id, domain["tne_bse"]->id);
  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Edge::bw()).id, domain["bne_bnw"]->id);
  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Edge::tw()).id, domain["tne_bnw"]->id);
  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Edge::sw()).id, domain["bne_tsw"]->id);

  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Edge::tn()).id, domain["tsw_tnw"]->id);
  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Edge::bn()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Edge::te()).id, domain["tsw_tse"]->id);
  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Edge::be()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Edge::ne()).id, domain["tsw_bne"]->id);

  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Edge::tn()).id, domain["tsw_tne"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Edge::bn()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Edge::bw()).id, domain["bsw_tsw"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Edge::te()).id, domain["tse_tsw"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Edge::be()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Edge::tw()).id, domain["tsw_tsw"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Edge::ne()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Edge::nw()).id, domain["tsw_bnw"]->id);

  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Edge::bs()).id, domain["bsw_tsw"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Edge::tn()).id, domain["tnw_tsw"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Edge::bn()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Edge::ts()).id, domain["tsw_tsw"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Edge::te()).id, domain["tsw_tne"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Edge::be()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Edge::ne()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Edge::se()).id, domain["tsw_bse"]->id);

  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Edge::bs()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Edge::tn()).id, domain["tnw_tse"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Edge::bn()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Edge::ts()).id, domain["tsw_tse"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Edge::bw()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Edge::te()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Edge::be()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Edge::tw()).id, domain["tsw_tnw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Edge::sw()).id, domain["tsw_bsw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Edge::ne()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Edge::se()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Edge::nw()).id, domain["tnw_bsw"]->id);

  CHECK_EQ(domain["tsw_tsw"]->getNormalNbrInfo(Edge::bn()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["tsw_tsw"]->getNormalNbrInfo(Edge::be()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["tsw_tsw"]->getNormalNbrInfo(Edge::ne()).id, domain["tsw_tne"]->id);

  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Edge::bn()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Edge::bw()).id, domain["tsw_bsw"]->id);
  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Edge::be()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Edge::ne()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Edge::nw()).id, domain["tsw_tnw"]->id);

  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Edge::bs()).id, domain["tsw_bsw"]->id);
  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Edge::bn()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Edge::be()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Edge::ne()).id, domain["tnw_tse"]->id);
  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Edge::se()).id, domain["tsw_tse"]->id);

  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Edge::bs()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Edge::bn()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Edge::bw()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Edge::be()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Edge::sw()).id, domain["tsw_tsw"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Edge::ne()).id, domain["tne_tsw"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Edge::se()).id, domain["tse_tsw"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Edge::nw()).id, domain["tnw_tsw"]->id);

  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Edge::tn()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Edge::bn()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Edge::bw()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Edge::te()).id, domain["tse_tse"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Edge::be()).id, domain["bse_tse"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Edge::tw()).id, domain["tsw_tse"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Edge::ne()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Edge::nw()).id, domain["tsw_bne"]->id);

  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Edge::tn()).id, domain["tse_tne"]->id);
  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Edge::bn()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Edge::bw()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Edge::tw()).id, domain["tse_tsw"]->id);
  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Edge::nw()).id, domain["tse_bnw"]->id);

  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Edge::bs()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Edge::tn()).id, domain["tne_tsw"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Edge::bn()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Edge::ts()).id, domain["tse_tsw"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Edge::bw()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Edge::te()).id, domain["tse_tne"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Edge::be()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Edge::tw()).id, domain["tsw_tne"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Edge::sw()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Edge::ne()).id, domain["tne_bse"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Edge::se()).id, domain["tse_bse"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Edge::nw()).id, domain["tnw_bse"]->id);

  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Edge::bs()).id, domain["bse_tse"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Edge::tn()).id, domain["tne_tse"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Edge::bn()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Edge::ts()).id, domain["tse_tse"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Edge::bw()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Edge::tw()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Edge::sw()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Edge::nw()).id, domain["tne_bsw"]->id);

  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Edge::bn()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Edge::bw()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Edge::be()).id, domain["tse_bse"]->id);
  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Edge::ne()).id, domain["tse_tne"]->id);
  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Edge::nw()).id, domain["tsw_tne"]->id);

  CHECK_EQ(domain["tse_tse"]->getNormalNbrInfo(Edge::bn()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["tse_tse"]->getNormalNbrInfo(Edge::bw()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["tse_tse"]->getNormalNbrInfo(Edge::nw()).id, domain["tse_tnw"]->id);

  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Edge::bs()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Edge::bn()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Edge::bw()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Edge::be()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Edge::sw()).id, domain["tsw_tse"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Edge::ne()).id, domain["tne_tse"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Edge::se()).id, domain["tse_tse"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Edge::nw()).id, domain["tnw_tse"]->id);

  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Edge::bs()).id, domain["tse_bse"]->id);
  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Edge::bn()).id, domain["tne_bse"]->id);
  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Edge::bw()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Edge::sw()).id, domain["tse_tsw"]->id);
  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Edge::nw()).id, domain["tne_tsw"]->id);

  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Edge::bs()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Edge::tn()).id, domain["tnw_tnw"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Edge::bn()).id, domain["bnw_tnw"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Edge::ts()).id, domain["tsw_tnw"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Edge::te()).id, domain["tnw_tse"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Edge::be()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Edge::ne()).id, domain["tnw_bne"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Edge::se()).id, domain["tsw_bne"]->id);

  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Edge::bs()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Edge::tn()).id, domain["tnw_tne"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Edge::bn()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Edge::ts()).id, domain["tsw_tne"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Edge::bw()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Edge::te()).id, domain["tne_tsw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Edge::be()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Edge::tw()).id, domain["tnw_tsw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Edge::sw()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Edge::ne()).id, domain["tne_bnw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Edge::se()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Edge::nw()).id, domain["tnw_bnw"]->id);

  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Edge::bs()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Edge::ts()).id, domain["tnw_tsw"]->id);
  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Edge::te()).id, domain["tnw_tne"]->id);
  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Edge::be()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Edge::se()).id, domain["tnw_bse"]->id);

  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Edge::bs()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Edge::ts()).id, domain["tnw_tse"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Edge::bw()).id, domain["bnw_tnw"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Edge::te()).id, domain["tne_tnw"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Edge::be()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Edge::tw()).id, domain["tnw_tnw"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Edge::sw()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Edge::se()).id, domain["tne_bsw"]->id);

  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Edge::bs()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Edge::bn()).id, domain["tnw_bnw"]->id);
  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Edge::be()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Edge::ne()).id, domain["tnw_tne"]->id);
  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Edge::se()).id, domain["tsw_tne"]->id);

  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Edge::bs()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Edge::bn()).id, domain["tnw_bne"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Edge::bw()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Edge::be()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Edge::sw()).id, domain["tsw_tnw"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Edge::ne()).id, domain["tne_tnw"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Edge::se()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Edge::nw()).id, domain["tnw_tnw"]->id);

  CHECK_EQ(domain["tnw_tnw"]->getNormalNbrInfo(Edge::bs()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["tnw_tnw"]->getNormalNbrInfo(Edge::be()).id, domain["tnw_bne"]->id);
  CHECK_EQ(domain["tnw_tnw"]->getNormalNbrInfo(Edge::se()).id, domain["tnw_tse"]->id);

  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Edge::bs()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Edge::bw()).id, domain["tnw_bnw"]->id);
  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Edge::be()).id, domain["tne_bnw"]->id);
  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Edge::sw()).id, domain["tnw_tsw"]->id);
  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Edge::se()).id, domain["tne_tsw"]->id);

  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Edge::bs()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Edge::tn()).id, domain["tne_tnw"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Edge::bn()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Edge::ts()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Edge::bw()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Edge::te()).id, domain["tne_tse"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Edge::be()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Edge::tw()).id, domain["tnw_tse"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Edge::sw()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Edge::ne()).id, domain["tne_bne"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Edge::se()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Edge::nw()).id, domain["tnw_bne"]->id);

  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Edge::bs()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Edge::tn()).id, domain["tne_tne"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Edge::bn()).id, domain["bne_tne"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Edge::ts()).id, domain["tse_tne"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Edge::bw()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Edge::tw()).id, domain["tne_tsw"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Edge::sw()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Edge::nw()).id, domain["tne_bnw"]->id);

  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Edge::bs()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Edge::ts()).id, domain["tne_tsw"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Edge::bw()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Edge::te()).id, domain["tne_tne"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Edge::be()).id, domain["bne_tne"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Edge::tw()).id, domain["tnw_tne"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Edge::sw()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Edge::se()).id, domain["tne_bse"]->id);

  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Edge::bs()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Edge::ts()).id, domain["tne_tse"]->id);
  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Edge::bw()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Edge::tw()).id, domain["tne_tnw"]->id);
  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Edge::sw()).id, domain["tne_bsw"]->id);

  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Edge::bs()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Edge::bn()).id, domain["tne_bnw"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Edge::bw()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Edge::be()).id, domain["tne_bse"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Edge::sw()).id, domain["tsw_tne"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Edge::ne()).id, domain["tne_tne"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Edge::se()).id, domain["tse_tne"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Edge::nw()).id, domain["tnw_tne"]->id);

  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Edge::bs()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Edge::bn()).id, domain["tne_bne"]->id);
  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Edge::bw()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Edge::sw()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Edge::nw()).id, domain["tne_tnw"]->id);

  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Edge::bs()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Edge::bw()).id, domain["tnw_bne"]->id);
  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Edge::be()).id, domain["tne_bne"]->id);
  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Edge::sw()).id, domain["tnw_tse"]->id);
  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Edge::se()).id, domain["tne_tse"]->id);

  CHECK_EQ(domain["tne_tne"]->getNormalNbrInfo(Edge::bs()).id, domain["tne_bse"]->id);
  CHECK_EQ(domain["tne_tne"]->getNormalNbrInfo(Edge::bw()).id, domain["tne_bnw"]->id);
  CHECK_EQ(domain["tne_tne"]->getNormalNbrInfo(Edge::sw()).id, domain["tne_tsw"]->id);
}
