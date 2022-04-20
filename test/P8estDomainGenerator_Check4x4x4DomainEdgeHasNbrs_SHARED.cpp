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
Check4x4x4DomainEdgeHasNeighbors(const PatchVector& domain)
{
  // edge hasnbr
  CHECK_UNARY_FALSE(domain["bsw_bsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bsw_bsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bsw_bsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["bsw_bsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bsw_bsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bsw_bsw"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bsw_bsw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["bsw_bsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bsw_bsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bsw_bsw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bsw_bsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bsw_bse"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bsw_bse"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bsw_bse"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bsw_bse"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bsw_bse"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bsw_bnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bsw_bnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bsw_bnw"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bsw_bnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bsw_bnw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bsw_bne"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bsw_bne"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bsw_bne"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bsw_bne"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bsw_bne"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bsw_bne"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bsw_bne"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bsw_bne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bsw_bne"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bsw_bne"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bsw_bne"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bsw_bne"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bsw_tsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["bsw_tsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bsw_tsw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["bsw_tsw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bsw_tsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bsw_tse"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bsw_tse"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["bsw_tse"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["bsw_tse"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["bsw_tse"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bsw_tse"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["bsw_tse"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bsw_tse"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bsw_tse"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bsw_tse"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bsw_tse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bsw_tse"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["bsw_tnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bsw_tnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["bsw_tnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bsw_tnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bsw_tnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bsw_tnw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["bsw_tnw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["bsw_tnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bsw_tnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bsw_tnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bsw_tnw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bsw_tnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bse_bsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bse_bsw"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bse_bsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bse_bsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bse_bsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bse_bse"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bse_bse"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bse_bse"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bse_bnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bse_bnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bse_bnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bse_bnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bse_bnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bse_bnw"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bse_bnw"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bse_bnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bse_bnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bse_bnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bse_bnw"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bse_bnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bse_bne"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bse_bne"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bse_bne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bse_bne"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bse_bne"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bse_tsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bse_tsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["bse_tsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["bse_tsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["bse_tsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bse_tsw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["bse_tsw"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bse_tsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bse_tsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bse_tsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bse_tsw"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bse_tsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bse_tse"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["bse_tse"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["bse_tse"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bse_tse"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bse_tse"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["bse_tne"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bse_tne"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["bse_tne"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bse_tne"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["bse_tne"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["bse_tne"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bse_tne"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bse_tne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bse_tne"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bse_tne"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bse_tne"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bse_tne"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bnw_bsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bnw_bsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bnw_bsw"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bnw_bsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bnw_bsw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bnw_bse"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bnw_bse"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bnw_bse"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bnw_bse"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bnw_bse"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bnw_bse"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bnw_bse"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bnw_bse"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bnw_bse"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bnw_bse"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bnw_bse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bnw_bse"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bnw_bnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bnw_bnw"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bnw_bnw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bnw_bne"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bnw_bne"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bnw_bne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bnw_bne"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bnw_bne"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["bnw_tsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bnw_tsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["bnw_tsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bnw_tsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bnw_tsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bnw_tsw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["bnw_tsw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["bnw_tsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bnw_tsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bnw_tsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bnw_tsw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bnw_tsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["bnw_tnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bnw_tnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bnw_tnw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["bnw_tnw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bnw_tnw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["bnw_tne"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["bnw_tne"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bnw_tne"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bnw_tne"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["bnw_tne"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bnw_tne"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["bnw_tne"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bnw_tne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bnw_tne"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bnw_tne"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bnw_tne"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bnw_tne"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bne_bsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bne_bsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bne_bsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bne_bsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bne_bsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bne_bsw"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bne_bsw"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bne_bsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bne_bsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bne_bsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bne_bsw"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bne_bsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bne_bse"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bne_bse"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bne_bse"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bne_bse"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bne_bse"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bne_bnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bne_bnw"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bne_bnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bne_bnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bne_bnw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bne_bne"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bne_bne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bne_bne"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["bne_tse"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["bne_tse"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["bne_tse"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bne_tse"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["bne_tse"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["bne_tse"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bne_tse"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bne_tse"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bne_tse"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bne_tse"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bne_tse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["bne_tse"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["bne_tnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["bne_tnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bne_tnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bne_tnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["bne_tnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["bne_tnw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["bne_tnw"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bne_tnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bne_tnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bne_tnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["bne_tnw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bne_tnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["bne_tne"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["bne_tne"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["bne_tne"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["bne_tne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["bne_tne"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["tsw_bsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tsw_bsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["tsw_bsw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tsw_bsw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tsw_bsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["tsw_bse"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["tsw_bse"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tsw_bse"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tsw_bse"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tsw_bse"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["tsw_bse"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tsw_bse"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["tsw_bse"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tsw_bse"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tsw_bse"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tsw_bse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tsw_bse"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tsw_bnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["tsw_bnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tsw_bnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["tsw_bnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["tsw_bnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["tsw_bnw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tsw_bnw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tsw_bnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tsw_bnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tsw_bnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tsw_bnw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tsw_bnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tsw_tsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tsw_tsw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tsw_tsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tsw_tse"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tsw_tse"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tsw_tse"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tsw_tse"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tsw_tse"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tsw_tnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tsw_tnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tsw_tnw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tsw_tnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tsw_tnw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tsw_tne"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tsw_tne"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tsw_tne"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tsw_tne"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tsw_tne"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tsw_tne"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tsw_tne"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tsw_tne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tsw_tne"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tsw_tne"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tsw_tne"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tsw_tne"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["tse_bsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["tse_bsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tse_bsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tse_bsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tse_bsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["tse_bsw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tse_bsw"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["tse_bsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tse_bsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tse_bsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tse_bsw"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tse_bsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["tse_bse"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tse_bse"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tse_bse"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["tse_bse"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tse_bse"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tse_bne"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["tse_bne"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tse_bne"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["tse_bne"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tse_bne"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tse_bne"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["tse_bne"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["tse_bne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tse_bne"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tse_bne"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tse_bne"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tse_bne"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tse_tsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tse_tsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tse_tsw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tse_tsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tse_tsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tse_tse"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tse_tse"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tse_tse"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tse_tnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tse_tnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tse_tnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tse_tnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tse_tnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tse_tnw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tse_tnw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tse_tnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tse_tnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tse_tnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tse_tnw"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tse_tnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tse_tne"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tse_tne"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tse_tne"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tse_tne"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tse_tne"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tnw_bsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["tnw_bsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tnw_bsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["tnw_bsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["tnw_bsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["tnw_bsw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tnw_bsw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tnw_bsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tnw_bsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tnw_bsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tnw_bsw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tnw_bsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tnw_bnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["tnw_bnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["tnw_bnw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tnw_bnw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tnw_bnw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tnw_bne"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tnw_bne"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["tnw_bne"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["tnw_bne"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tnw_bne"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["tnw_bne"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tnw_bne"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["tnw_bne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tnw_bne"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tnw_bne"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tnw_bne"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tnw_bne"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tnw_tsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tnw_tsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tnw_tsw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tnw_tsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tnw_tsw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tnw_tse"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tnw_tse"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tnw_tse"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tnw_tse"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tnw_tse"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tnw_tse"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tnw_tse"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tnw_tse"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tnw_tse"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tnw_tse"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tnw_tse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tnw_tse"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tnw_tnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tnw_tnw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tnw_tnw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tnw_tne"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tnw_tne"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tnw_tne"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tnw_tne"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tnw_tne"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tne_bse"]->hasNbr(Edge::bs()));
  CHECK_UNARY(domain["tne_bse"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tne_bse"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["tne_bse"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tne_bse"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tne_bse"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["tne_bse"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["tne_bse"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tne_bse"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tne_bse"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tne_bse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tne_bse"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tne_bnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tne_bnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["tne_bnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["tne_bnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tne_bnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY(domain["tne_bnw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tne_bnw"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["tne_bnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tne_bnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tne_bnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tne_bnw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tne_bnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tne_bne"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Edge::bn()));
  CHECK_UNARY(domain["tne_bne"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tne_bne"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Edge::be()));
  CHECK_UNARY(domain["tne_bne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tne_bne"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tne_tsw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tne_tsw"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tne_tsw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tne_tsw"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tne_tsw"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tne_tsw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tne_tsw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tne_tsw"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tne_tsw"]->hasNbr(Edge::sw()));
  CHECK_UNARY(domain["tne_tsw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tne_tsw"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tne_tsw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tne_tse"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Edge::tn()));
  CHECK_UNARY(domain["tne_tse"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tne_tse"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tne_tse"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Edge::se()));
  CHECK_UNARY(domain["tne_tse"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tne_tnw"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tne_tnw"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Edge::te()));
  CHECK_UNARY(domain["tne_tnw"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tne_tnw"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Edge::ne()));
  CHECK_UNARY(domain["tne_tnw"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Edge::nw()));

  CHECK_UNARY(domain["tne_tne"]->hasNbr(Edge::bs()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Edge::tn()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Edge::bn()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Edge::ts()));
  CHECK_UNARY(domain["tne_tne"]->hasNbr(Edge::bw()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Edge::te()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Edge::be()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Edge::tw()));
  CHECK_UNARY(domain["tne_tne"]->hasNbr(Edge::sw()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Edge::ne()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Edge::se()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Edge::nw()));
}
