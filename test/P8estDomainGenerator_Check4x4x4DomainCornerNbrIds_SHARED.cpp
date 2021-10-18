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
Check4x4x4DomainCornerNeighborIds(const PatchVector& domain)
{
  // corner nbr ids
  CHECK(domain["bsw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["bsw_tne"]->id);

  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["bsw_tnw"]->id);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["bse_tnw"]->id);

  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["bsw_tse"]->id);
  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["bnw_tse"]->id);

  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["bsw_tsw"]->id);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["bse_tsw"]->id);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["bnw_tsw"]->id);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["bne_tsw"]->id);

  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bsw_bne"]->id);
  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tsw_bne"]->id);

  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bsw_bnw"]->id);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bse_bnw"]->id);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tsw_bnw"]->id);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tse_bnw"]->id);

  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bsw_bse"]->id);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bnw_bse"]->id);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tsw_bse"]->id);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tnw_bse"]->id);

  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bsw_bsw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bse_bsw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bnw_bsw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bne_bsw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tsw_bsw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tse_bsw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tnw_bsw"]->id);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tne_bsw"]->id);

  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["bsw_tne"]->id);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["bse_tne"]->id);

  CHECK(domain["bse_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["bse_tnw"]->id);

  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["bsw_tse"]->id);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["bse_tse"]->id);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["bnw_tse"]->id);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["bne_tse"]->id);

  CHECK(domain["bse_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["bse_tsw"]->id);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["bne_tsw"]->id);

  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bsw_bne"]->id);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bse_bne"]->id);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tsw_bne"]->id);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tse_bne"]->id);

  CHECK(domain["bse_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bse_bnw"]->id);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tse_bnw"]->id);

  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bsw_bse"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bse_bse"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bnw_bse"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bne_bse"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tsw_bse"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tse_bse"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tnw_bse"]->id);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tne_bse"]->id);

  CHECK(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bse_bsw"]->id);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bne_bsw"]->id);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tse_bsw"]->id);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tne_bsw"]->id);

  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["bsw_tne"]->id);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["bnw_tne"]->id);

  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["bsw_tnw"]->id);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["bse_tnw"]->id);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["bnw_tnw"]->id);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["bne_tnw"]->id);

  CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["bnw_tse"]->id);

  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["bnw_tsw"]->id);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["bne_tsw"]->id);

  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bsw_bne"]->id);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bnw_bne"]->id);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tsw_bne"]->id);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tnw_bne"]->id);

  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bsw_bnw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bse_bnw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bnw_bnw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bne_bnw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tsw_bnw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tse_bnw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tnw_bnw"]->id);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tne_bnw"]->id);

  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bnw_bse"]->id);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tnw_bse"]->id);

  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bnw_bsw"]->id);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bne_bsw"]->id);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tnw_bsw"]->id);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tne_bsw"]->id);

  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["bsw_tne"]->id);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["bse_tne"]->id);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["bnw_tne"]->id);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["bne_tne"]->id);

  CHECK(domain["bne_bse"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["bse_tnw"]->id);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["bne_tnw"]->id);

  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["bnw_tse"]->id);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["bne_tse"]->id);

  CHECK(domain["bne_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["bne_tsw"]->id);

  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bsw_bne"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bse_bne"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bnw_bne"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bne_bne"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tsw_bne"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tse_bne"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tnw_bne"]->id);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tne_bne"]->id);

  CHECK(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bse_bnw"]->id);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bne_bnw"]->id);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tse_bnw"]->id);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tne_bnw"]->id);

  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bnw_bse"]->id);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bne_bse"]->id);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tnw_bse"]->id);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tne_bse"]->id);

  CHECK(domain["bne_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bne_bsw"]->id);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tne_bsw"]->id);

  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bsw_tne"]->id);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tsw_tne"]->id);

  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bsw_tnw"]->id);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bse_tnw"]->id);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tsw_tnw"]->id);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tse_tnw"]->id);

  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bsw_tse"]->id);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bnw_tse"]->id);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tsw_tse"]->id);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tnw_tse"]->id);

  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bsw_tsw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bse_tsw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bnw_tsw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bne_tsw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tsw_tsw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tse_tsw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tnw_tsw"]->id);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tne_tsw"]->id);

  CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["tsw_bne"]->id);

  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["tsw_bnw"]->id);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["tse_bnw"]->id);

  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["tsw_bse"]->id);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["tnw_bse"]->id);

  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["tsw_bsw"]->id);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["tse_bsw"]->id);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["tnw_bsw"]->id);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["tne_bsw"]->id);

  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bsw_tne"]->id);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bse_tne"]->id);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tsw_tne"]->id);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tse_tne"]->id);

  CHECK(domain["tse_bse"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bse_tnw"]->id);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tse_tnw"]->id);

  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bsw_tse"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bse_tse"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bnw_tse"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bne_tse"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tsw_tse"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tse_tse"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tnw_tse"]->id);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tne_tse"]->id);

  CHECK(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bse_tsw"]->id);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bne_tsw"]->id);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tse_tsw"]->id);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tne_tsw"]->id);

  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["tsw_bne"]->id);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["tse_bne"]->id);

  CHECK(domain["tse_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["tse_bnw"]->id);

  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["tsw_bse"]->id);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["tse_bse"]->id);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["tnw_bse"]->id);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["tne_bse"]->id);

  CHECK(domain["tse_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["tse_bsw"]->id);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["tne_bsw"]->id);

  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bsw_tne"]->id);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bnw_tne"]->id);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tsw_tne"]->id);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tnw_tne"]->id);

  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bsw_tnw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bse_tnw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bnw_tnw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bne_tnw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tsw_tnw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tse_tnw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tnw_tnw"]->id);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tne_tnw"]->id);

  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bnw_tse"]->id);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tnw_tse"]->id);

  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bnw_tsw"]->id);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bne_tsw"]->id);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tnw_tsw"]->id);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tne_tsw"]->id);

  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["tsw_bne"]->id);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["tnw_bne"]->id);

  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["tsw_bnw"]->id);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["tse_bnw"]->id);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["tnw_bnw"]->id);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["tne_bnw"]->id);

  CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["tnw_bse"]->id);

  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["tnw_bsw"]->id);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["tne_bsw"]->id);

  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bsw_tne"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bse_tne"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bnw_tne"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["bne_tne"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tsw_tne"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tse_tne"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tnw_tne"]->id);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id == domain["tne_tne"]->id);

  CHECK(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bse_tnw"]->id);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["bne_tnw"]->id);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tse_tnw"]->id);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id == domain["tne_tnw"]->id);

  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bnw_tse"]->id);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["bne_tse"]->id);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tnw_tse"]->id);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id == domain["tne_tse"]->id);

  CHECK(domain["tne_bne"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["bne_tsw"]->id);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id == domain["tne_tsw"]->id);

  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["tsw_bne"]->id);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["tse_bne"]->id);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["tnw_bne"]->id);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id == domain["tne_bne"]->id);

  CHECK(domain["tne_tse"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["tse_bnw"]->id);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id == domain["tne_bnw"]->id);

  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["tnw_bse"]->id);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id == domain["tne_bse"]->id);

  CHECK(domain["tne_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id == domain["tne_bsw"]->id);
}
void
Check4x4x4DomainCornerNeighborRanks(const PatchVector& domain)
{
  // corner nbr ranks
  CHECK(domain["bsw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["bsw_tne"]->rank);

  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["bsw_bse"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["bse_tnw"]->rank);

  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["bnw_tse"]->rank);

  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["bsw_tsw"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["bne_tsw"]->rank);

  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tsw_bne"]->rank);

  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tse_bnw"]->rank);

  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tnw_bse"]->rank);

  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bsw_bsw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bse_bsw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bnw_bsw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tsw_bsw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tne_bsw"]->rank);

  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["bse_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["bse_tne"]->rank);

  CHECK(domain["bse_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["bse_tnw"]->rank);

  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["bse_tse"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["bne_tse"]->rank);

  CHECK(domain["bse_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["bse_bne"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["bne_tsw"]->rank);

  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bse_bne"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tse_bne"]->rank);

  CHECK(domain["bse_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bse_tse"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tse_bnw"]->rank);

  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bsw_bse"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bse_bse"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bne_bse"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tse_bse"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tne_bse"]->rank);

  CHECK(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bse_bsw"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tne_bsw"]->rank);

  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["bnw_tne"]->rank);

  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["bnw_tnw"]->rank);
  CHECK(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["bne_tnw"]->rank);

  CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["bnw_tse"]->rank);

  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["bnw_bne"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["bne_tsw"]->rank);

  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bnw_bne"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tnw_bne"]->rank);

  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bsw_bnw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bnw_bnw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bne_bnw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tnw_bnw"]->rank);
  CHECK(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tne_bnw"]->rank);

  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tnw_bse"]->rank);

  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bnw_bsw"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tne_bsw"]->rank);

  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["bse_tne"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["bne_tne"]->rank);

  CHECK(domain["bne_bse"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["bne_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["bne_tnw"]->rank);

  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["bne_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["bne_tse"]->rank);

  CHECK(domain["bne_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["bne_tsw"]->rank);

  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bsw_bne"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bse_bne"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bnw_bne"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bne_bne"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tse_bne"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tnw_bne"]->rank);
  CHECK(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tne_bne"]->rank);

  CHECK(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bse_bnw"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bne_bnw"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tne_bnw"]->rank);

  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bnw_bse"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bne_bse"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tne_bse"]->rank);

  CHECK(domain["bne_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bne_bsw"]->rank);
  CHECK(domain["bne_tne"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tne_bsw"]->rank);

  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tsw_tne"]->rank);

  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tsw_tnw"]->rank);
  CHECK(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tse_tnw"]->rank);

  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tsw_tse"]->rank);
  CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tnw_tse"]->rank);

  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bsw_tsw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tsw_tsw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tse_tsw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tnw_tsw"]->rank);
  CHECK(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tne_tsw"]->rank);

  CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["tsw_bne"]->rank);

  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["tsw_tse"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["tse_bnw"]->rank);

  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["tnw_bse"]->rank);

  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["tsw_bsw"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["tne_bsw"]->rank);

  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bse_tne"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tsw_tne"]->rank);
  CHECK(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tse_tne"]->rank);

  CHECK(domain["tse_bse"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["tse_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tse_tnw"]->rank);

  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bsw_tse"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bse_tse"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bne_tse"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tsw_tse"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tse_tse"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tnw_tse"]->rank);
  CHECK(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tne_tse"]->rank);

  CHECK(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bse_tsw"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tse_tsw"]->rank);
  CHECK(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tne_tsw"]->rank);

  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["tse_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["tse_bne"]->rank);

  CHECK(domain["tse_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["tse_bnw"]->rank);

  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["tsw_bse"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["tse_bse"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["tne_bse"]->rank);

  CHECK(domain["tse_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["tse_bsw"]->rank);
  CHECK(domain["tse_tne"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["tne_bsw"]->rank);

  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tsw_tne"]->rank);
  CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tnw_tne"]->rank);

  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bsw_tnw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bnw_tnw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tsw_tnw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tnw_tnw"]->rank);
  CHECK(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tne_tnw"]->rank);

  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tnw_tse"]->rank);

  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bnw_tsw"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tnw_tsw"]->rank);
  CHECK(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tne_tsw"]->rank);

  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["tnw_bne"]->rank);

  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["tsw_bnw"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["tnw_bnw"]->rank);
  CHECK(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["tne_bnw"]->rank);

  CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["tnw_bse"]->rank);

  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["tnw_bsw"]->rank);
  CHECK(domain["tnw_tne"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["tne_bsw"]->rank);

  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bsw_tne"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bse_tne"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bnw_tne"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["bne_tne"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tsw_tne"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tse_tne"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tnw_tne"]->rank);
  CHECK(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank == domain["tne_tne"]->rank);

  CHECK(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bse_tnw"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["bne_tnw"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tse_tnw"]->rank);
  CHECK(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank == domain["tne_tnw"]->rank);

  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bnw_tse"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["bne_tse"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tnw_tse"]->rank);
  CHECK(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank == domain["tne_tse"]->rank);

  CHECK(domain["tne_bne"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["bne_tsw"]->rank);
  CHECK(domain["tne_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank == domain["tne_tsw"]->rank);

  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["tsw_bne"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["tse_bne"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["tnw_bne"]->rank);
  CHECK(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank == domain["tne_bne"]->rank);

  CHECK(domain["tne_tse"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["tse_bnw"]->rank);
  CHECK(domain["tne_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank == domain["tne_bnw"]->rank);

  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["tnw_bse"]->rank);
  CHECK(domain["tne_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank == domain["tne_bse"]->rank);

  CHECK(domain["tne_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank == domain["tne_bsw"]->rank);
}
