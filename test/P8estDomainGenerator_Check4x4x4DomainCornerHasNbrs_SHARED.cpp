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
Check4x4x4DomainCornerHasNeighbors(const PatchVector& domain)
{
  // corner hasNbr
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bsw_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bsw_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bsw_bse"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bsw_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bsw_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bsw_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bsw_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bsw_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bsw_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bsw_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bsw_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bsw_bne"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bsw_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bsw_bne"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bsw_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bsw_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["bsw_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bsw_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bsw_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bsw_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bsw_tse"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["bsw_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["bsw_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bsw_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bsw_tse"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bsw_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bsw_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bsw_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["bsw_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bsw_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["bsw_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bsw_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bsw_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bsw_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bsw_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["bsw_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["bsw_tne"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["bsw_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["bsw_tne"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bsw_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bsw_tne"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bsw_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bsw_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bse_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bse_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bse_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bse_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bse_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bse_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bse_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bse_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bse_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bse_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bse_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bse_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bse_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bse_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bse_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bse_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bse_bne"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bse_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bse_bne"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bse_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bse_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bse_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bse_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["bse_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["bse_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bse_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bse_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bse_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bse_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bse_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bse_tse"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["bse_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bse_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bse_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bse_tse"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bse_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bse_tse"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["bse_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["bse_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["bse_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["bse_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bse_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bse_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bse_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bse_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["bse_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bse_tne"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["bse_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bse_tne"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bse_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bse_tne"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bse_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bse_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bnw_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bnw_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bnw_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bnw_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bnw_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bnw_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bnw_bse"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bnw_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bnw_bse"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bnw_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bnw_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bnw_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bnw_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bnw_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bnw_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bnw_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bnw_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["bnw_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bnw_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["bnw_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bnw_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bnw_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bnw_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bnw_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["bnw_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["bnw_tse"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["bnw_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["bnw_tse"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bnw_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bnw_tse"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bnw_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bnw_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["bnw_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bnw_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bnw_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["bnw_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["bnw_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bnw_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bnw_tne"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bnw_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bnw_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bnw_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bnw_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bne_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bne_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bne_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bne_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bne_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bne_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bne_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bne_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bne_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bne_bse"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bne_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bne_bse"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bne_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bne_bse"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bne_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bne_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bne_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bne_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bne_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bne_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bne_bne"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["bne_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["bne_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["bne_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["bne_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bne_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bne_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bne_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["bne_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["bne_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bne_tse"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["bne_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bne_tse"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bne_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bne_tse"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["bne_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bne_tse"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["bne_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["bne_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bne_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bne_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bne_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["bne_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bne_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bne_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["bne_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["bne_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["bne_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tsw_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tsw_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["tsw_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tsw_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tsw_bse"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tsw_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tsw_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tsw_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tsw_bse"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["tsw_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["tsw_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tsw_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tsw_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tsw_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tsw_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tsw_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["tsw_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tsw_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["tsw_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tsw_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tsw_bne"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tsw_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tsw_bne"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["tsw_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["tsw_bne"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["tsw_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["tsw_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tsw_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tsw_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tsw_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tsw_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tsw_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tsw_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tsw_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tsw_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tsw_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tsw_tne"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tsw_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tsw_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tsw_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tsw_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tsw_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tsw_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tse_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tse_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tse_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tse_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tse_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tse_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["tse_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["tse_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tse_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tse_bse"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tse_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tse_bse"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tse_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tse_bse"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["tse_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tse_bse"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tse_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tse_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tse_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tse_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["tse_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["tse_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["tse_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["tse_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tse_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tse_bne"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tse_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tse_bne"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["tse_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tse_bne"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["tse_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tse_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tse_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tse_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tse_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tse_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tse_tse"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tse_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tse_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tse_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tse_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tse_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tse_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tse_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tse_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tse_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tse_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tse_tne"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tnw_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tnw_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tnw_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tnw_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tnw_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["tnw_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tnw_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["tnw_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tnw_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tnw_bse"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tnw_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tnw_bse"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["tnw_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["tnw_bse"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["tnw_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["tnw_bse"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tnw_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["tnw_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tnw_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tnw_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tnw_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tnw_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tnw_bne"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["tnw_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["tnw_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tnw_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tnw_bne"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tnw_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tnw_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tnw_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tnw_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tnw_tse"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tnw_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tnw_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tnw_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tnw_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tnw_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tnw_tse"]->hasNbr(Corner<3>::tne()));

  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tnw_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tnw_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tnw_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tnw_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tnw_tne"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tne_bsw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tne_bsw"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tne_bsw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tne_bsw"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["tne_bsw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["tne_bsw"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["tne_bsw"]->hasNbr(Corner<3>::tnw()));
  CHECK(domain["tne_bsw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tne_bse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tne_bse"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tne_bse"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tne_bse"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["tne_bse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tne_bse"]->hasNbr(Corner<3>::tse()));
  CHECK(domain["tne_bse"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tne_bse"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tne_bnw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tne_bnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tne_bnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tne_bnw"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["tne_bnw"]->hasNbr(Corner<3>::tsw()));
  CHECK(domain["tne_bnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tne_bnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tne_bnw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tne_bne"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Corner<3>::bne()));
  CHECK(domain["tne_bne"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tne_bne"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tne_tsw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tne_tsw"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tne_tsw"]->hasNbr(Corner<3>::bnw()));
  CHECK(domain["tne_tsw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tne_tsw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tne_tsw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tne_tsw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tne_tsw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tne_tse"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Corner<3>::bse()));
  CHECK(domain["tne_tse"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tne_tse"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tne_tnw"]->hasNbr(Corner<3>::bsw()));
  CHECK(domain["tne_tnw"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tne_tnw"]->hasNbr(Corner<3>::tne()));

  CHECK(domain["tne_tne"]->hasNbr(Corner<3>::bsw()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::bse()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::bnw()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::bne()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::tsw()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::tse()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::tnw()));
  CHECK_FALSE(domain["tne_tne"]->hasNbr(Corner<3>::tne()));
}