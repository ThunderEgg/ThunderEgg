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
Check4x4x4RefinedBSWDomainCornerNeighborIds(const PatchVector& domain)
{
  CHECK_EQ(domain["bsw_bsw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["bsw_bsw_tne"]->id);

  CHECK_EQ(domain["bsw_bsw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["bsw_bsw_tnw"]->id);

  CHECK_EQ(domain["bsw_bsw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["bsw_bsw_tse"]->id);

  CHECK_EQ(domain["bsw_bsw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["bsw_bsw_tsw"]->id);

  CHECK_EQ(domain["bsw_bsw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bsw_bsw_bne"]->id);

  CHECK_EQ(domain["bsw_bsw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bsw_bsw_bnw"]->id);

  CHECK_EQ(domain["bsw_bsw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bsw_bsw_bse"]->id);

  CHECK_EQ(domain["bsw_bsw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bsw_bsw_bsw"]->id);
  CHECK_EQ(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Corner<3>::tne()).id, domain["bsw_tne"]->id);

  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["bse_tnw"]->id);

  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["bnw_tse"]->id);

  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["bsw_tsw"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["bne_tsw"]->id);

  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tsw_bne"]->id);

  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bsw_bnw"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tse_bnw"]->id);

  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bsw_bse"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tnw_bse"]->id);

  CHECK_EQ(domain["bsw_tne"]->getFineNbrInfo(Corner<3>::bsw()).ids[0], domain["bsw_bsw_tne"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bse_bsw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bnw_bsw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tsw_bsw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tne_bsw"]->id);

  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["bse_tne"]->id);

  CHECK_EQ(domain["bse_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["bse_tnw"]->id);

  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["bse_tse"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["bne_tse"]->id);

  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["bne_tsw"]->id);

  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bse_bne"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tse_bne"]->id);

  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tse_bnw"]->id);

  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bsw_bse"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bse_bse"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bne_bse"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tse_bse"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tne_bse"]->id);

  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bse_bsw"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tne_bsw"]->id);

  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["bnw_tne"]->id);

  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["bnw_tnw"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["bne_tnw"]->id);

  CHECK_EQ(domain["bnw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["bnw_tse"]->id);

  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["bne_tsw"]->id);

  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bnw_bne"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tnw_bne"]->id);

  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bsw_bnw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bnw_bnw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bne_bnw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tnw_bnw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tne_bnw"]->id);

  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tnw_bse"]->id);

  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bnw_bsw"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tne_bsw"]->id);

  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["bne_tne"]->id);

  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["bne_tnw"]->id);

  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["bne_tse"]->id);

  CHECK_EQ(domain["bne_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["bne_tsw"]->id);

  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bse_bne"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bnw_bne"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bne_bne"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tnw_bne"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tne_bne"]->id);

  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bne_bnw"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tne_bnw"]->id);

  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bne_bse"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tne_bse"]->id);

  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tne_bsw"]->id);

  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tsw_tne"]->id);

  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tsw_tnw"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tse_tnw"]->id);

  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tsw_tse"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tnw_tse"]->id);

  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bsw_tsw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tsw_tsw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tse_tsw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tnw_tsw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tne_tsw"]->id);

  CHECK_EQ(domain["tsw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["tsw_bne"]->id);

  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["tse_bnw"]->id);

  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["tnw_bse"]->id);

  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["tsw_bsw"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["tne_bsw"]->id);

  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tsw_tne"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tse_tne"]->id);

  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tse_tnw"]->id);

  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bse_tse"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tsw_tse"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tse_tse"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tnw_tse"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tne_tse"]->id);

  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tse_tsw"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tne_tsw"]->id);

  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["tse_bne"]->id);

  CHECK_EQ(domain["tse_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["tse_bnw"]->id);

  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["tse_bse"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["tne_bse"]->id);

  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["tne_bsw"]->id);

  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tsw_tne"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tnw_tne"]->id);

  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bnw_tnw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tsw_tnw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tnw_tnw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tne_tnw"]->id);

  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tnw_tse"]->id);

  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tnw_tsw"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tne_tsw"]->id);

  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["tnw_bne"]->id);

  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["tnw_bnw"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["tne_bnw"]->id);

  CHECK_EQ(domain["tnw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["tnw_bse"]->id);

  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["tne_bsw"]->id);

  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["bne_tne"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tsw_tne"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tse_tne"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tnw_tne"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tne()).id, domain["tne_tne"]->id);

  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::tnw()).id, domain["tne_tnw"]->id);

  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tnw_tse"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::tse()).id, domain["tne_tse"]->id);

  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Corner<3>::tsw()).id, domain["tne_tsw"]->id);

  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["tnw_bne"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bne()).id, domain["tne_bne"]->id);

  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Corner<3>::bnw()).id, domain["tne_bnw"]->id);

  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Corner<3>::bse()).id, domain["tne_bse"]->id);

  CHECK_EQ(domain["tne_tne"]->getNormalNbrInfo(Corner<3>::bsw()).id, domain["tne_bsw"]->id);
}
