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
Check4x4x4RefinedBSWDomainCornerHasNeighbors(const PatchVector& domain)
{
  CHECK_UNARY_FALSE(domain["bsw_bsw_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bsw_bsw_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bsw_bsw_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bsw_bsw_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bsw_bsw_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bsw_bsw_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bsw_bsw_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bsw_bsw_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bsw_bsw_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["bsw_bsw_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bsw_bsw_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["bsw_bsw_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bsw_bsw_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["bsw_bsw_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["bsw_bsw_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bsw_bsw_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bsw_bsw_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bsw_bse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bsw_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bsw_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bsw_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bsw_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bsw_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bsw_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bsw_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bsw_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bsw_bne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bsw_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bsw_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bsw_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bsw_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["bsw_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bsw_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bsw_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bsw_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bsw_tse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["bsw_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["bsw_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bsw_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bsw_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bsw_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bsw_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bsw_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["bsw_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bsw_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["bsw_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bsw_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bsw_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bsw_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bsw_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bsw_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bse_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bse_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bse_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bse_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bse_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bse_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bse_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bse_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bse_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bse_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bse_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bse_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bse_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bse_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bse_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bse_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bse_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["bse_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["bse_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bse_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bse_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bse_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bse_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["bse_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bse_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bse_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bse_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["bse_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bse_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["bse_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bse_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bse_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bse_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bse_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bse_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bnw_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bnw_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bnw_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bnw_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bnw_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bnw_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bnw_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bnw_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bnw_bse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bnw_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bnw_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bnw_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bnw_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bnw_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bnw_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bnw_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["bnw_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bnw_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["bnw_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bnw_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bnw_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bnw_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bnw_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bnw_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["bnw_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bnw_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bnw_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["bnw_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["bnw_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bnw_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bnw_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bnw_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bnw_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bnw_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bnw_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bne_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bne_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bne_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bne_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bne_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bne_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bne_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bne_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bne_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bne_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bne_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bne_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bne_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bne_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bne_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["bne_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["bne_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bne_tse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["bne_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bne_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bne_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bne_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["bne_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bne_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["bne_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["bne_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bne_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bne_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bne_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["bne_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bne_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bne_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["bne_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["bne_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["bne_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tsw_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tsw_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["tsw_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tsw_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tsw_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tsw_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tsw_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tsw_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tsw_bse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["tsw_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["tsw_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tsw_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tsw_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tsw_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tsw_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tsw_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["tsw_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tsw_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["tsw_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["tsw_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tsw_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tsw_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tsw_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tsw_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tsw_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tsw_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tsw_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tsw_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tsw_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tsw_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tsw_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tsw_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tsw_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tsw_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tsw_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tse_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tse_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tse_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tse_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tse_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tse_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["tse_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["tse_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tse_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["tse_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tse_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["tse_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tse_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tse_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tse_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tse_bne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["tse_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tse_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["tse_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tse_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tse_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tse_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tse_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tse_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tse_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tse_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tse_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tse_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tse_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tse_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tse_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tse_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tse_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tse_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tse_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tnw_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tnw_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tnw_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tnw_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tnw_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["tnw_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tnw_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["tnw_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["tnw_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tnw_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["tnw_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tnw_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tnw_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tnw_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tnw_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tnw_bne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["tnw_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["tnw_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tnw_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tnw_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tnw_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tnw_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tnw_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tnw_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tnw_tse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tnw_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tnw_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tnw_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tnw_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tnw_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tnw_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tnw_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tnw_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tnw_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tnw_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY(domain["tne_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tne_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tne_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tne_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tne_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["tne_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tne_bse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY(domain["tne_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tne_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tne_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tne_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tne_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tne_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["tne_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY(domain["tne_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tne_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tne_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tne_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY(domain["tne_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tne_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tne_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tne_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tne_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY(domain["tne_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tne_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tne_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tne_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tne_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tne_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY(domain["tne_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tne_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tne_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY(domain["tne_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tne_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK_UNARY(domain["tne_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_UNARY_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::tne()));
}
