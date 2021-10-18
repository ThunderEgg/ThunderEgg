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
Check4x4x4DomainEdgeHasNeighbors(const PatchVector& domain)
{
  // edge hasnbr
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Edge::bs()));
  CHECK(domain["bsw_bsw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Edge::bw()));
  CHECK(domain["bsw_bsw"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Edge::sw()));
  CHECK(domain["bsw_bsw"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Edge::bs()));
  CHECK(domain["bsw_bse"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Edge::bw()));
  CHECK(domain["bsw_bse"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Edge::be()));
  CHECK(domain["bsw_bse"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Edge::sw()));
  CHECK(domain["bsw_bse"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Edge::se()));
  CHECK(domain["bsw_bse"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Edge::bs()));
  CHECK(domain["bsw_bnw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Edge::bn()));
  CHECK(domain["bsw_bnw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Edge::bw()));
  CHECK(domain["bsw_bnw"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Edge::sw()));
  CHECK(domain["bsw_bnw"]->hasNbr(Edge::ne()));
  CHECK(domain["bsw_bnw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bsw_bne"]->hasNbr(Edge::bs()));
  CHECK(domain["bsw_bne"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bsw_bne"]->hasNbr(Edge::bn()));
  CHECK(domain["bsw_bne"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bsw_bne"]->hasNbr(Edge::bw()));
  CHECK(domain["bsw_bne"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bsw_bne"]->hasNbr(Edge::be()));
  CHECK(domain["bsw_bne"]->hasNbr(Edge::tw()));
  CHECK(domain["bsw_bne"]->hasNbr(Edge::sw()));
  CHECK(domain["bsw_bne"]->hasNbr(Edge::ne()));
  CHECK(domain["bsw_bne"]->hasNbr(Edge::se()));
  CHECK(domain["bsw_bne"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Edge::bs()));
  CHECK(domain["bsw_tsw"]->hasNbr(Edge::tn()));
  CHECK(domain["bsw_tsw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Edge::bw()));
  CHECK(domain["bsw_tsw"]->hasNbr(Edge::te()));
  CHECK(domain["bsw_tsw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Edge::sw()));
  CHECK(domain["bsw_tsw"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bsw_tse"]->hasNbr(Edge::bs()));
  CHECK(domain["bsw_tse"]->hasNbr(Edge::tn()));
  CHECK(domain["bsw_tse"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["bsw_tse"]->hasNbr(Edge::ts()));
  CHECK(domain["bsw_tse"]->hasNbr(Edge::bw()));
  CHECK(domain["bsw_tse"]->hasNbr(Edge::te()));
  CHECK(domain["bsw_tse"]->hasNbr(Edge::be()));
  CHECK(domain["bsw_tse"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bsw_tse"]->hasNbr(Edge::sw()));
  CHECK(domain["bsw_tse"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bsw_tse"]->hasNbr(Edge::se()));
  CHECK(domain["bsw_tse"]->hasNbr(Edge::nw()));

  CHECK(domain["bsw_tnw"]->hasNbr(Edge::bs()));
  CHECK(domain["bsw_tnw"]->hasNbr(Edge::tn()));
  CHECK(domain["bsw_tnw"]->hasNbr(Edge::bn()));
  CHECK(domain["bsw_tnw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bsw_tnw"]->hasNbr(Edge::bw()));
  CHECK(domain["bsw_tnw"]->hasNbr(Edge::te()));
  CHECK(domain["bsw_tnw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["bsw_tnw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bsw_tnw"]->hasNbr(Edge::sw()));
  CHECK(domain["bsw_tnw"]->hasNbr(Edge::ne()));
  CHECK(domain["bsw_tnw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bsw_tnw"]->hasNbr(Edge::nw()));

  CHECK(domain["bsw_tne"]->hasNbr(Edge::bs()));
  CHECK(domain["bsw_tne"]->hasNbr(Edge::tn()));
  CHECK(domain["bsw_tne"]->hasNbr(Edge::bn()));
  CHECK(domain["bsw_tne"]->hasNbr(Edge::ts()));
  CHECK(domain["bsw_tne"]->hasNbr(Edge::bw()));
  CHECK(domain["bsw_tne"]->hasNbr(Edge::te()));
  CHECK(domain["bsw_tne"]->hasNbr(Edge::be()));
  CHECK(domain["bsw_tne"]->hasNbr(Edge::tw()));
  CHECK(domain["bsw_tne"]->hasNbr(Edge::sw()));
  CHECK(domain["bsw_tne"]->hasNbr(Edge::ne()));
  CHECK(domain["bsw_tne"]->hasNbr(Edge::se()));
  CHECK(domain["bsw_tne"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Edge::bs()));
  CHECK(domain["bse_bsw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Edge::bw()));
  CHECK(domain["bse_bsw"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Edge::be()));
  CHECK(domain["bse_bsw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Edge::sw()));
  CHECK(domain["bse_bsw"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Edge::se()));
  CHECK(domain["bse_bsw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bse_bse"]->hasNbr(Edge::bs()));
  CHECK(domain["bse_bse"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Edge::be()));
  CHECK(domain["bse_bse"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Edge::se()));
  CHECK(domain["bse_bse"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bse_bnw"]->hasNbr(Edge::bs()));
  CHECK(domain["bse_bnw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bse_bnw"]->hasNbr(Edge::bn()));
  CHECK(domain["bse_bnw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bse_bnw"]->hasNbr(Edge::bw()));
  CHECK(domain["bse_bnw"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bse_bnw"]->hasNbr(Edge::be()));
  CHECK(domain["bse_bnw"]->hasNbr(Edge::tw()));
  CHECK(domain["bse_bnw"]->hasNbr(Edge::sw()));
  CHECK(domain["bse_bnw"]->hasNbr(Edge::ne()));
  CHECK(domain["bse_bnw"]->hasNbr(Edge::se()));
  CHECK(domain["bse_bnw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bse_bne"]->hasNbr(Edge::bs()));
  CHECK(domain["bse_bne"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bse_bne"]->hasNbr(Edge::bn()));
  CHECK(domain["bse_bne"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bse_bne"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["bse_bne"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bse_bne"]->hasNbr(Edge::be()));
  CHECK(domain["bse_bne"]->hasNbr(Edge::tw()));
  CHECK(domain["bse_bne"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bse_bne"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bse_bne"]->hasNbr(Edge::se()));
  CHECK(domain["bse_bne"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bse_tsw"]->hasNbr(Edge::bs()));
  CHECK(domain["bse_tsw"]->hasNbr(Edge::tn()));
  CHECK(domain["bse_tsw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["bse_tsw"]->hasNbr(Edge::ts()));
  CHECK(domain["bse_tsw"]->hasNbr(Edge::bw()));
  CHECK(domain["bse_tsw"]->hasNbr(Edge::te()));
  CHECK(domain["bse_tsw"]->hasNbr(Edge::be()));
  CHECK(domain["bse_tsw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bse_tsw"]->hasNbr(Edge::sw()));
  CHECK(domain["bse_tsw"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bse_tsw"]->hasNbr(Edge::se()));
  CHECK(domain["bse_tsw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bse_tse"]->hasNbr(Edge::bs()));
  CHECK(domain["bse_tse"]->hasNbr(Edge::tn()));
  CHECK(domain["bse_tse"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["bse_tse"]->hasNbr(Edge::ts()));
  CHECK(domain["bse_tse"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["bse_tse"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bse_tse"]->hasNbr(Edge::be()));
  CHECK(domain["bse_tse"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bse_tse"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bse_tse"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bse_tse"]->hasNbr(Edge::se()));
  CHECK(domain["bse_tse"]->hasNbr(Edge::nw()));

  CHECK(domain["bse_tnw"]->hasNbr(Edge::bs()));
  CHECK(domain["bse_tnw"]->hasNbr(Edge::tn()));
  CHECK(domain["bse_tnw"]->hasNbr(Edge::bn()));
  CHECK(domain["bse_tnw"]->hasNbr(Edge::ts()));
  CHECK(domain["bse_tnw"]->hasNbr(Edge::bw()));
  CHECK(domain["bse_tnw"]->hasNbr(Edge::te()));
  CHECK(domain["bse_tnw"]->hasNbr(Edge::be()));
  CHECK(domain["bse_tnw"]->hasNbr(Edge::tw()));
  CHECK(domain["bse_tnw"]->hasNbr(Edge::sw()));
  CHECK(domain["bse_tnw"]->hasNbr(Edge::ne()));
  CHECK(domain["bse_tnw"]->hasNbr(Edge::se()));
  CHECK(domain["bse_tnw"]->hasNbr(Edge::nw()));

  CHECK(domain["bse_tne"]->hasNbr(Edge::bs()));
  CHECK(domain["bse_tne"]->hasNbr(Edge::tn()));
  CHECK(domain["bse_tne"]->hasNbr(Edge::bn()));
  CHECK(domain["bse_tne"]->hasNbr(Edge::ts()));
  CHECK(domain["bse_tne"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["bse_tne"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bse_tne"]->hasNbr(Edge::be()));
  CHECK(domain["bse_tne"]->hasNbr(Edge::tw()));
  CHECK(domain["bse_tne"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bse_tne"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bse_tne"]->hasNbr(Edge::se()));
  CHECK(domain["bse_tne"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Edge::bs()));
  CHECK(domain["bnw_bsw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Edge::bn()));
  CHECK(domain["bnw_bsw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Edge::bw()));
  CHECK(domain["bnw_bsw"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Edge::sw()));
  CHECK(domain["bnw_bsw"]->hasNbr(Edge::ne()));
  CHECK(domain["bnw_bsw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bnw_bse"]->hasNbr(Edge::bs()));
  CHECK(domain["bnw_bse"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bnw_bse"]->hasNbr(Edge::bn()));
  CHECK(domain["bnw_bse"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bnw_bse"]->hasNbr(Edge::bw()));
  CHECK(domain["bnw_bse"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bnw_bse"]->hasNbr(Edge::be()));
  CHECK(domain["bnw_bse"]->hasNbr(Edge::tw()));
  CHECK(domain["bnw_bse"]->hasNbr(Edge::sw()));
  CHECK(domain["bnw_bse"]->hasNbr(Edge::ne()));
  CHECK(domain["bnw_bse"]->hasNbr(Edge::se()));
  CHECK(domain["bnw_bse"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Edge::bn()));
  CHECK(domain["bnw_bnw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Edge::bw()));
  CHECK(domain["bnw_bnw"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Edge::ne()));
  CHECK(domain["bnw_bnw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Edge::bn()));
  CHECK(domain["bnw_bne"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Edge::bw()));
  CHECK(domain["bnw_bne"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Edge::be()));
  CHECK(domain["bnw_bne"]->hasNbr(Edge::tw()));
  CHECK(domain["bnw_bne"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Edge::ne()));
  CHECK(domain["bnw_bne"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Edge::nw()));

  CHECK(domain["bnw_tsw"]->hasNbr(Edge::bs()));
  CHECK(domain["bnw_tsw"]->hasNbr(Edge::tn()));
  CHECK(domain["bnw_tsw"]->hasNbr(Edge::bn()));
  CHECK(domain["bnw_tsw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bnw_tsw"]->hasNbr(Edge::bw()));
  CHECK(domain["bnw_tsw"]->hasNbr(Edge::te()));
  CHECK(domain["bnw_tsw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["bnw_tsw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bnw_tsw"]->hasNbr(Edge::sw()));
  CHECK(domain["bnw_tsw"]->hasNbr(Edge::ne()));
  CHECK(domain["bnw_tsw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bnw_tsw"]->hasNbr(Edge::nw()));

  CHECK(domain["bnw_tse"]->hasNbr(Edge::bs()));
  CHECK(domain["bnw_tse"]->hasNbr(Edge::tn()));
  CHECK(domain["bnw_tse"]->hasNbr(Edge::bn()));
  CHECK(domain["bnw_tse"]->hasNbr(Edge::ts()));
  CHECK(domain["bnw_tse"]->hasNbr(Edge::bw()));
  CHECK(domain["bnw_tse"]->hasNbr(Edge::te()));
  CHECK(domain["bnw_tse"]->hasNbr(Edge::be()));
  CHECK(domain["bnw_tse"]->hasNbr(Edge::tw()));
  CHECK(domain["bnw_tse"]->hasNbr(Edge::sw()));
  CHECK(domain["bnw_tse"]->hasNbr(Edge::ne()));
  CHECK(domain["bnw_tse"]->hasNbr(Edge::se()));
  CHECK(domain["bnw_tse"]->hasNbr(Edge::nw()));

  CHECK(domain["bnw_tnw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Edge::bn()));
  CHECK(domain["bnw_tnw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Edge::bw()));
  CHECK(domain["bnw_tnw"]->hasNbr(Edge::te()));
  CHECK(domain["bnw_tnw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Edge::ne()));
  CHECK(domain["bnw_tnw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Edge::nw()));

  CHECK(domain["bnw_tne"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["bnw_tne"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bnw_tne"]->hasNbr(Edge::bn()));
  CHECK(domain["bnw_tne"]->hasNbr(Edge::ts()));
  CHECK(domain["bnw_tne"]->hasNbr(Edge::bw()));
  CHECK(domain["bnw_tne"]->hasNbr(Edge::te()));
  CHECK(domain["bnw_tne"]->hasNbr(Edge::be()));
  CHECK(domain["bnw_tne"]->hasNbr(Edge::tw()));
  CHECK(domain["bnw_tne"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bnw_tne"]->hasNbr(Edge::ne()));
  CHECK(domain["bnw_tne"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bnw_tne"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bne_bsw"]->hasNbr(Edge::bs()));
  CHECK(domain["bne_bsw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bne_bsw"]->hasNbr(Edge::bn()));
  CHECK(domain["bne_bsw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bne_bsw"]->hasNbr(Edge::bw()));
  CHECK(domain["bne_bsw"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bne_bsw"]->hasNbr(Edge::be()));
  CHECK(domain["bne_bsw"]->hasNbr(Edge::tw()));
  CHECK(domain["bne_bsw"]->hasNbr(Edge::sw()));
  CHECK(domain["bne_bsw"]->hasNbr(Edge::ne()));
  CHECK(domain["bne_bsw"]->hasNbr(Edge::se()));
  CHECK(domain["bne_bsw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bne_bse"]->hasNbr(Edge::bs()));
  CHECK(domain["bne_bse"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bne_bse"]->hasNbr(Edge::bn()));
  CHECK(domain["bne_bse"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bne_bse"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["bne_bse"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bne_bse"]->hasNbr(Edge::be()));
  CHECK(domain["bne_bse"]->hasNbr(Edge::tw()));
  CHECK(domain["bne_bse"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bne_bse"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bne_bse"]->hasNbr(Edge::se()));
  CHECK(domain["bne_bse"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Edge::bn()));
  CHECK(domain["bne_bnw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Edge::bw()));
  CHECK(domain["bne_bnw"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Edge::be()));
  CHECK(domain["bne_bnw"]->hasNbr(Edge::tw()));
  CHECK(domain["bne_bnw"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Edge::ne()));
  CHECK(domain["bne_bnw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["bne_bne"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Edge::bn()));
  CHECK(domain["bne_bne"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Edge::be()));
  CHECK(domain["bne_bne"]->hasNbr(Edge::tw()));
  CHECK(domain["bne_bne"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Edge::nw()));

  CHECK(domain["bne_tsw"]->hasNbr(Edge::bs()));
  CHECK(domain["bne_tsw"]->hasNbr(Edge::tn()));
  CHECK(domain["bne_tsw"]->hasNbr(Edge::bn()));
  CHECK(domain["bne_tsw"]->hasNbr(Edge::ts()));
  CHECK(domain["bne_tsw"]->hasNbr(Edge::bw()));
  CHECK(domain["bne_tsw"]->hasNbr(Edge::te()));
  CHECK(domain["bne_tsw"]->hasNbr(Edge::be()));
  CHECK(domain["bne_tsw"]->hasNbr(Edge::tw()));
  CHECK(domain["bne_tsw"]->hasNbr(Edge::sw()));
  CHECK(domain["bne_tsw"]->hasNbr(Edge::ne()));
  CHECK(domain["bne_tsw"]->hasNbr(Edge::se()));
  CHECK(domain["bne_tsw"]->hasNbr(Edge::nw()));

  CHECK(domain["bne_tse"]->hasNbr(Edge::bs()));
  CHECK(domain["bne_tse"]->hasNbr(Edge::tn()));
  CHECK(domain["bne_tse"]->hasNbr(Edge::bn()));
  CHECK(domain["bne_tse"]->hasNbr(Edge::ts()));
  CHECK(domain["bne_tse"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["bne_tse"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bne_tse"]->hasNbr(Edge::be()));
  CHECK(domain["bne_tse"]->hasNbr(Edge::tw()));
  CHECK(domain["bne_tse"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bne_tse"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bne_tse"]->hasNbr(Edge::se()));
  CHECK(domain["bne_tse"]->hasNbr(Edge::nw()));

  CHECK(domain["bne_tnw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["bne_tnw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bne_tnw"]->hasNbr(Edge::bn()));
  CHECK(domain["bne_tnw"]->hasNbr(Edge::ts()));
  CHECK(domain["bne_tnw"]->hasNbr(Edge::bw()));
  CHECK(domain["bne_tnw"]->hasNbr(Edge::te()));
  CHECK(domain["bne_tnw"]->hasNbr(Edge::be()));
  CHECK(domain["bne_tnw"]->hasNbr(Edge::tw()));
  CHECK(domain["bne_tnw"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bne_tnw"]->hasNbr(Edge::ne()));
  CHECK(domain["bne_tnw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bne_tnw"]->hasNbr(Edge::nw()));

  CHECK(domain["bne_tne"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Edge::bn()));
  CHECK(domain["bne_tne"]->hasNbr(Edge::ts()));
  CHECK(domain["bne_tne"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Edge::be()));
  CHECK(domain["bne_tne"]->hasNbr(Edge::tw()));
  CHECK(domain["bne_tne"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Edge::bs()));
  CHECK(domain["tsw_bsw"]->hasNbr(Edge::tn()));
  CHECK(domain["tsw_bsw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Edge::bw()));
  CHECK(domain["tsw_bsw"]->hasNbr(Edge::te()));
  CHECK(domain["tsw_bsw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Edge::sw()));
  CHECK(domain["tsw_bsw"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["tsw_bse"]->hasNbr(Edge::bs()));
  CHECK(domain["tsw_bse"]->hasNbr(Edge::tn()));
  CHECK(domain["tsw_bse"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tsw_bse"]->hasNbr(Edge::ts()));
  CHECK(domain["tsw_bse"]->hasNbr(Edge::bw()));
  CHECK(domain["tsw_bse"]->hasNbr(Edge::te()));
  CHECK(domain["tsw_bse"]->hasNbr(Edge::be()));
  CHECK(domain["tsw_bse"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tsw_bse"]->hasNbr(Edge::sw()));
  CHECK(domain["tsw_bse"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tsw_bse"]->hasNbr(Edge::se()));
  CHECK(domain["tsw_bse"]->hasNbr(Edge::nw()));

  CHECK(domain["tsw_bnw"]->hasNbr(Edge::bs()));
  CHECK(domain["tsw_bnw"]->hasNbr(Edge::tn()));
  CHECK(domain["tsw_bnw"]->hasNbr(Edge::bn()));
  CHECK(domain["tsw_bnw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["tsw_bnw"]->hasNbr(Edge::bw()));
  CHECK(domain["tsw_bnw"]->hasNbr(Edge::te()));
  CHECK(domain["tsw_bnw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tsw_bnw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tsw_bnw"]->hasNbr(Edge::sw()));
  CHECK(domain["tsw_bnw"]->hasNbr(Edge::ne()));
  CHECK(domain["tsw_bnw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tsw_bnw"]->hasNbr(Edge::nw()));

  CHECK(domain["tsw_bne"]->hasNbr(Edge::bs()));
  CHECK(domain["tsw_bne"]->hasNbr(Edge::tn()));
  CHECK(domain["tsw_bne"]->hasNbr(Edge::bn()));
  CHECK(domain["tsw_bne"]->hasNbr(Edge::ts()));
  CHECK(domain["tsw_bne"]->hasNbr(Edge::bw()));
  CHECK(domain["tsw_bne"]->hasNbr(Edge::te()));
  CHECK(domain["tsw_bne"]->hasNbr(Edge::be()));
  CHECK(domain["tsw_bne"]->hasNbr(Edge::tw()));
  CHECK(domain["tsw_bne"]->hasNbr(Edge::sw()));
  CHECK(domain["tsw_bne"]->hasNbr(Edge::ne()));
  CHECK(domain["tsw_bne"]->hasNbr(Edge::se()));
  CHECK(domain["tsw_bne"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Edge::tn()));
  CHECK(domain["tsw_tsw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Edge::te()));
  CHECK(domain["tsw_tsw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Edge::sw()));
  CHECK(domain["tsw_tsw"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Edge::tn()));
  CHECK(domain["tsw_tse"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Edge::ts()));
  CHECK(domain["tsw_tse"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Edge::te()));
  CHECK(domain["tsw_tse"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Edge::sw()));
  CHECK(domain["tsw_tse"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Edge::se()));
  CHECK(domain["tsw_tse"]->hasNbr(Edge::nw()));

  CHECK(domain["tsw_tnw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Edge::tn()));
  CHECK(domain["tsw_tnw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Edge::te()));
  CHECK(domain["tsw_tnw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Edge::sw()));
  CHECK(domain["tsw_tnw"]->hasNbr(Edge::ne()));
  CHECK(domain["tsw_tnw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Edge::nw()));

  CHECK(domain["tsw_tne"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tsw_tne"]->hasNbr(Edge::tn()));
  CHECK(domain["tsw_tne"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tsw_tne"]->hasNbr(Edge::ts()));
  CHECK(domain["tsw_tne"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tsw_tne"]->hasNbr(Edge::te()));
  CHECK(domain["tsw_tne"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tsw_tne"]->hasNbr(Edge::tw()));
  CHECK(domain["tsw_tne"]->hasNbr(Edge::sw()));
  CHECK(domain["tsw_tne"]->hasNbr(Edge::ne()));
  CHECK(domain["tsw_tne"]->hasNbr(Edge::se()));
  CHECK(domain["tsw_tne"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["tse_bsw"]->hasNbr(Edge::bs()));
  CHECK(domain["tse_bsw"]->hasNbr(Edge::tn()));
  CHECK(domain["tse_bsw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tse_bsw"]->hasNbr(Edge::ts()));
  CHECK(domain["tse_bsw"]->hasNbr(Edge::bw()));
  CHECK(domain["tse_bsw"]->hasNbr(Edge::te()));
  CHECK(domain["tse_bsw"]->hasNbr(Edge::be()));
  CHECK(domain["tse_bsw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tse_bsw"]->hasNbr(Edge::sw()));
  CHECK(domain["tse_bsw"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tse_bsw"]->hasNbr(Edge::se()));
  CHECK(domain["tse_bsw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["tse_bse"]->hasNbr(Edge::bs()));
  CHECK(domain["tse_bse"]->hasNbr(Edge::tn()));
  CHECK(domain["tse_bse"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tse_bse"]->hasNbr(Edge::ts()));
  CHECK(domain["tse_bse"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tse_bse"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["tse_bse"]->hasNbr(Edge::be()));
  CHECK(domain["tse_bse"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tse_bse"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tse_bse"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tse_bse"]->hasNbr(Edge::se()));
  CHECK(domain["tse_bse"]->hasNbr(Edge::nw()));

  CHECK(domain["tse_bnw"]->hasNbr(Edge::bs()));
  CHECK(domain["tse_bnw"]->hasNbr(Edge::tn()));
  CHECK(domain["tse_bnw"]->hasNbr(Edge::bn()));
  CHECK(domain["tse_bnw"]->hasNbr(Edge::ts()));
  CHECK(domain["tse_bnw"]->hasNbr(Edge::bw()));
  CHECK(domain["tse_bnw"]->hasNbr(Edge::te()));
  CHECK(domain["tse_bnw"]->hasNbr(Edge::be()));
  CHECK(domain["tse_bnw"]->hasNbr(Edge::tw()));
  CHECK(domain["tse_bnw"]->hasNbr(Edge::sw()));
  CHECK(domain["tse_bnw"]->hasNbr(Edge::ne()));
  CHECK(domain["tse_bnw"]->hasNbr(Edge::se()));
  CHECK(domain["tse_bnw"]->hasNbr(Edge::nw()));

  CHECK(domain["tse_bne"]->hasNbr(Edge::bs()));
  CHECK(domain["tse_bne"]->hasNbr(Edge::tn()));
  CHECK(domain["tse_bne"]->hasNbr(Edge::bn()));
  CHECK(domain["tse_bne"]->hasNbr(Edge::ts()));
  CHECK(domain["tse_bne"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tse_bne"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["tse_bne"]->hasNbr(Edge::be()));
  CHECK(domain["tse_bne"]->hasNbr(Edge::tw()));
  CHECK(domain["tse_bne"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tse_bne"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tse_bne"]->hasNbr(Edge::se()));
  CHECK(domain["tse_bne"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Edge::tn()));
  CHECK(domain["tse_tsw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Edge::ts()));
  CHECK(domain["tse_tsw"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Edge::te()));
  CHECK(domain["tse_tsw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Edge::sw()));
  CHECK(domain["tse_tsw"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Edge::se()));
  CHECK(domain["tse_tsw"]->hasNbr(Edge::nw()));

  CHECK_FALSE(domain["tse_tse"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Edge::tn()));
  CHECK(domain["tse_tse"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Edge::ts()));
  CHECK(domain["tse_tse"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Edge::se()));
  CHECK(domain["tse_tse"]->hasNbr(Edge::nw()));

  CHECK(domain["tse_tnw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tse_tnw"]->hasNbr(Edge::tn()));
  CHECK(domain["tse_tnw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tse_tnw"]->hasNbr(Edge::ts()));
  CHECK(domain["tse_tnw"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tse_tnw"]->hasNbr(Edge::te()));
  CHECK(domain["tse_tnw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tse_tnw"]->hasNbr(Edge::tw()));
  CHECK(domain["tse_tnw"]->hasNbr(Edge::sw()));
  CHECK(domain["tse_tnw"]->hasNbr(Edge::ne()));
  CHECK(domain["tse_tnw"]->hasNbr(Edge::se()));
  CHECK(domain["tse_tnw"]->hasNbr(Edge::nw()));

  CHECK(domain["tse_tne"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Edge::tn()));
  CHECK(domain["tse_tne"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Edge::ts()));
  CHECK(domain["tse_tne"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Edge::tw()));
  CHECK(domain["tse_tne"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Edge::se()));
  CHECK(domain["tse_tne"]->hasNbr(Edge::nw()));

  CHECK(domain["tnw_bsw"]->hasNbr(Edge::bs()));
  CHECK(domain["tnw_bsw"]->hasNbr(Edge::tn()));
  CHECK(domain["tnw_bsw"]->hasNbr(Edge::bn()));
  CHECK(domain["tnw_bsw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["tnw_bsw"]->hasNbr(Edge::bw()));
  CHECK(domain["tnw_bsw"]->hasNbr(Edge::te()));
  CHECK(domain["tnw_bsw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tnw_bsw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tnw_bsw"]->hasNbr(Edge::sw()));
  CHECK(domain["tnw_bsw"]->hasNbr(Edge::ne()));
  CHECK(domain["tnw_bsw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tnw_bsw"]->hasNbr(Edge::nw()));

  CHECK(domain["tnw_bse"]->hasNbr(Edge::bs()));
  CHECK(domain["tnw_bse"]->hasNbr(Edge::tn()));
  CHECK(domain["tnw_bse"]->hasNbr(Edge::bn()));
  CHECK(domain["tnw_bse"]->hasNbr(Edge::ts()));
  CHECK(domain["tnw_bse"]->hasNbr(Edge::bw()));
  CHECK(domain["tnw_bse"]->hasNbr(Edge::te()));
  CHECK(domain["tnw_bse"]->hasNbr(Edge::be()));
  CHECK(domain["tnw_bse"]->hasNbr(Edge::tw()));
  CHECK(domain["tnw_bse"]->hasNbr(Edge::sw()));
  CHECK(domain["tnw_bse"]->hasNbr(Edge::ne()));
  CHECK(domain["tnw_bse"]->hasNbr(Edge::se()));
  CHECK(domain["tnw_bse"]->hasNbr(Edge::nw()));

  CHECK(domain["tnw_bnw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Edge::bn()));
  CHECK(domain["tnw_bnw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Edge::bw()));
  CHECK(domain["tnw_bnw"]->hasNbr(Edge::te()));
  CHECK(domain["tnw_bnw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Edge::ne()));
  CHECK(domain["tnw_bnw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Edge::nw()));

  CHECK(domain["tnw_bne"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tnw_bne"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["tnw_bne"]->hasNbr(Edge::bn()));
  CHECK(domain["tnw_bne"]->hasNbr(Edge::ts()));
  CHECK(domain["tnw_bne"]->hasNbr(Edge::bw()));
  CHECK(domain["tnw_bne"]->hasNbr(Edge::te()));
  CHECK(domain["tnw_bne"]->hasNbr(Edge::be()));
  CHECK(domain["tnw_bne"]->hasNbr(Edge::tw()));
  CHECK(domain["tnw_bne"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tnw_bne"]->hasNbr(Edge::ne()));
  CHECK(domain["tnw_bne"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tnw_bne"]->hasNbr(Edge::nw()));

  CHECK(domain["tnw_tsw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Edge::tn()));
  CHECK(domain["tnw_tsw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Edge::te()));
  CHECK(domain["tnw_tsw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Edge::sw()));
  CHECK(domain["tnw_tsw"]->hasNbr(Edge::ne()));
  CHECK(domain["tnw_tsw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Edge::nw()));

  CHECK(domain["tnw_tse"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tnw_tse"]->hasNbr(Edge::tn()));
  CHECK(domain["tnw_tse"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tnw_tse"]->hasNbr(Edge::ts()));
  CHECK(domain["tnw_tse"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tnw_tse"]->hasNbr(Edge::te()));
  CHECK(domain["tnw_tse"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tnw_tse"]->hasNbr(Edge::tw()));
  CHECK(domain["tnw_tse"]->hasNbr(Edge::sw()));
  CHECK(domain["tnw_tse"]->hasNbr(Edge::ne()));
  CHECK(domain["tnw_tse"]->hasNbr(Edge::se()));
  CHECK(domain["tnw_tse"]->hasNbr(Edge::nw()));

  CHECK(domain["tnw_tnw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Edge::ts()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Edge::te()));
  CHECK(domain["tnw_tnw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Edge::tw()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Edge::ne()));
  CHECK(domain["tnw_tnw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Edge::nw()));

  CHECK(domain["tnw_tne"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Edge::ts()));
  CHECK(domain["tnw_tne"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Edge::te()));
  CHECK(domain["tnw_tne"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Edge::tw()));
  CHECK(domain["tnw_tne"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Edge::ne()));
  CHECK(domain["tnw_tne"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Edge::nw()));

  CHECK(domain["tne_bsw"]->hasNbr(Edge::bs()));
  CHECK(domain["tne_bsw"]->hasNbr(Edge::tn()));
  CHECK(domain["tne_bsw"]->hasNbr(Edge::bn()));
  CHECK(domain["tne_bsw"]->hasNbr(Edge::ts()));
  CHECK(domain["tne_bsw"]->hasNbr(Edge::bw()));
  CHECK(domain["tne_bsw"]->hasNbr(Edge::te()));
  CHECK(domain["tne_bsw"]->hasNbr(Edge::be()));
  CHECK(domain["tne_bsw"]->hasNbr(Edge::tw()));
  CHECK(domain["tne_bsw"]->hasNbr(Edge::sw()));
  CHECK(domain["tne_bsw"]->hasNbr(Edge::ne()));
  CHECK(domain["tne_bsw"]->hasNbr(Edge::se()));
  CHECK(domain["tne_bsw"]->hasNbr(Edge::nw()));

  CHECK(domain["tne_bse"]->hasNbr(Edge::bs()));
  CHECK(domain["tne_bse"]->hasNbr(Edge::tn()));
  CHECK(domain["tne_bse"]->hasNbr(Edge::bn()));
  CHECK(domain["tne_bse"]->hasNbr(Edge::ts()));
  CHECK(domain["tne_bse"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tne_bse"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["tne_bse"]->hasNbr(Edge::be()));
  CHECK(domain["tne_bse"]->hasNbr(Edge::tw()));
  CHECK(domain["tne_bse"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tne_bse"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tne_bse"]->hasNbr(Edge::se()));
  CHECK(domain["tne_bse"]->hasNbr(Edge::nw()));

  CHECK(domain["tne_bnw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tne_bnw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["tne_bnw"]->hasNbr(Edge::bn()));
  CHECK(domain["tne_bnw"]->hasNbr(Edge::ts()));
  CHECK(domain["tne_bnw"]->hasNbr(Edge::bw()));
  CHECK(domain["tne_bnw"]->hasNbr(Edge::te()));
  CHECK(domain["tne_bnw"]->hasNbr(Edge::be()));
  CHECK(domain["tne_bnw"]->hasNbr(Edge::tw()));
  CHECK(domain["tne_bnw"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tne_bnw"]->hasNbr(Edge::ne()));
  CHECK(domain["tne_bnw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tne_bnw"]->hasNbr(Edge::nw()));

  CHECK(domain["tne_bne"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Edge::bn()));
  CHECK(domain["tne_bne"]->hasNbr(Edge::ts()));
  CHECK(domain["tne_bne"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Edge::be()));
  CHECK(domain["tne_bne"]->hasNbr(Edge::tw()));
  CHECK(domain["tne_bne"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Edge::nw()));

  CHECK(domain["tne_tsw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tne_tsw"]->hasNbr(Edge::tn()));
  CHECK(domain["tne_tsw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tne_tsw"]->hasNbr(Edge::ts()));
  CHECK(domain["tne_tsw"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tne_tsw"]->hasNbr(Edge::te()));
  CHECK(domain["tne_tsw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tne_tsw"]->hasNbr(Edge::tw()));
  CHECK(domain["tne_tsw"]->hasNbr(Edge::sw()));
  CHECK(domain["tne_tsw"]->hasNbr(Edge::ne()));
  CHECK(domain["tne_tsw"]->hasNbr(Edge::se()));
  CHECK(domain["tne_tsw"]->hasNbr(Edge::nw()));

  CHECK(domain["tne_tse"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Edge::tn()));
  CHECK(domain["tne_tse"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Edge::ts()));
  CHECK(domain["tne_tse"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Edge::tw()));
  CHECK(domain["tne_tse"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Edge::se()));
  CHECK(domain["tne_tse"]->hasNbr(Edge::nw()));

  CHECK(domain["tne_tnw"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Edge::ts()));
  CHECK(domain["tne_tnw"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Edge::te()));
  CHECK(domain["tne_tnw"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Edge::tw()));
  CHECK(domain["tne_tnw"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Edge::ne()));
  CHECK(domain["tne_tnw"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Edge::nw()));

  CHECK(domain["tne_tne"]->hasNbr(Edge::bs()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Edge::tn()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Edge::bn()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Edge::ts()));
  CHECK(domain["tne_tne"]->hasNbr(Edge::bw()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Edge::te()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Edge::be()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Edge::tw()));
  CHECK(domain["tne_tne"]->hasNbr(Edge::sw()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Edge::ne()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Edge::se()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Edge::nw()));
}