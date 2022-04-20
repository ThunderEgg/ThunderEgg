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
Check4x4x4RefinedBSWDomainCornerNeighborRanks(const PatchVector& domain)
{
  CHECK_EQ(domain["bsw_bsw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["bsw_bsw_tne"]->rank);

  CHECK_EQ(domain["bsw_bsw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["bsw_bsw_tnw"]->rank);

  CHECK_EQ(domain["bsw_bsw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["bsw_bsw_tse"]->rank);

  CHECK_EQ(domain["bsw_bsw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["bsw_bsw_tsw"]->rank);

  CHECK_EQ(domain["bsw_bsw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bsw_bsw_bne"]->rank);

  CHECK_EQ(domain["bsw_bsw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bsw_bsw_bnw"]->rank);

  CHECK_EQ(domain["bsw_bsw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bsw_bsw_bse"]->rank);

  CHECK_EQ(domain["bsw_bsw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bsw_bsw_bsw"]->rank);
  CHECK_EQ(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Corner<3>::tne()).rank, domain["bsw_tne"]->rank);

  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["bsw_tnw"]->rank);
  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["bse_tnw"]->rank);

  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["bsw_tse"]->rank);
  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["bnw_tse"]->rank);

  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["bsw_tsw"]->rank);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["bse_tsw"]->rank);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["bnw_tsw"]->rank);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["bne_tsw"]->rank);

  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bsw_bne"]->rank);
  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tsw_bne"]->rank);

  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bsw_bnw"]->rank);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bse_bnw"]->rank);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tsw_bnw"]->rank);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tse_bnw"]->rank);

  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bsw_bse"]->rank);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bnw_bse"]->rank);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tsw_bse"]->rank);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tnw_bse"]->rank);

  CHECK_EQ(domain["bsw_tne"]->getFineNbrInfo(Corner<3>::bsw()).ranks[0], domain["bsw_bsw_tne"]->rank);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bse_bsw"]->rank);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bnw_bsw"]->rank);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bne_bsw"]->rank);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tsw_bsw"]->rank);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tse_bsw"]->rank);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tnw_bsw"]->rank);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tne_bsw"]->rank);

  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["bsw_tne"]->rank);
  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["bse_tne"]->rank);

  CHECK_EQ(domain["bse_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["bse_tnw"]->rank);

  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["bsw_tse"]->rank);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["bse_tse"]->rank);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["bnw_tse"]->rank);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["bne_tse"]->rank);

  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["bse_tsw"]->rank);
  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["bne_tsw"]->rank);

  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bsw_bne"]->rank);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bse_bne"]->rank);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tsw_bne"]->rank);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tse_bne"]->rank);

  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bse_bnw"]->rank);
  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tse_bnw"]->rank);

  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bsw_bse"]->rank);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bse_bse"]->rank);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bnw_bse"]->rank);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bne_bse"]->rank);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tsw_bse"]->rank);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tse_bse"]->rank);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tnw_bse"]->rank);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tne_bse"]->rank);

  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bse_bsw"]->rank);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bne_bsw"]->rank);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tse_bsw"]->rank);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tne_bsw"]->rank);

  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["bsw_tne"]->rank);
  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["bnw_tne"]->rank);

  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["bsw_tnw"]->rank);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["bse_tnw"]->rank);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["bnw_tnw"]->rank);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["bne_tnw"]->rank);

  CHECK_EQ(domain["bnw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["bnw_tse"]->rank);

  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["bnw_tsw"]->rank);
  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["bne_tsw"]->rank);

  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bsw_bne"]->rank);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bnw_bne"]->rank);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tsw_bne"]->rank);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tnw_bne"]->rank);

  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bsw_bnw"]->rank);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bse_bnw"]->rank);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bnw_bnw"]->rank);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bne_bnw"]->rank);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tsw_bnw"]->rank);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tse_bnw"]->rank);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tnw_bnw"]->rank);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tne_bnw"]->rank);

  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bnw_bse"]->rank);
  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tnw_bse"]->rank);

  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bnw_bsw"]->rank);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bne_bsw"]->rank);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tnw_bsw"]->rank);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tne_bsw"]->rank);

  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["bsw_tne"]->rank);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["bse_tne"]->rank);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["bnw_tne"]->rank);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["bne_tne"]->rank);

  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["bse_tnw"]->rank);
  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["bne_tnw"]->rank);

  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["bnw_tse"]->rank);
  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["bne_tse"]->rank);

  CHECK_EQ(domain["bne_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["bne_tsw"]->rank);

  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bsw_bne"]->rank);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bse_bne"]->rank);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bnw_bne"]->rank);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bne_bne"]->rank);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tsw_bne"]->rank);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tse_bne"]->rank);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tnw_bne"]->rank);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tne_bne"]->rank);

  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bse_bnw"]->rank);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bne_bnw"]->rank);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tse_bnw"]->rank);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tne_bnw"]->rank);

  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bnw_bse"]->rank);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bne_bse"]->rank);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tnw_bse"]->rank);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tne_bse"]->rank);

  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bne_bsw"]->rank);
  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tne_bsw"]->rank);

  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bsw_tne"]->rank);
  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tsw_tne"]->rank);

  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bsw_tnw"]->rank);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bse_tnw"]->rank);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tsw_tnw"]->rank);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tse_tnw"]->rank);

  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bsw_tse"]->rank);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bnw_tse"]->rank);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tsw_tse"]->rank);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tnw_tse"]->rank);

  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bsw_tsw"]->rank);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bse_tsw"]->rank);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bnw_tsw"]->rank);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bne_tsw"]->rank);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tsw_tsw"]->rank);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tse_tsw"]->rank);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tnw_tsw"]->rank);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tne_tsw"]->rank);

  CHECK_EQ(domain["tsw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["tsw_bne"]->rank);

  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["tsw_bnw"]->rank);
  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["tse_bnw"]->rank);

  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["tsw_bse"]->rank);
  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["tnw_bse"]->rank);

  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["tsw_bsw"]->rank);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["tse_bsw"]->rank);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["tnw_bsw"]->rank);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["tne_bsw"]->rank);

  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bsw_tne"]->rank);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bse_tne"]->rank);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tsw_tne"]->rank);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tse_tne"]->rank);

  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bse_tnw"]->rank);
  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tse_tnw"]->rank);

  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bsw_tse"]->rank);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bse_tse"]->rank);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bnw_tse"]->rank);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bne_tse"]->rank);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tsw_tse"]->rank);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tse_tse"]->rank);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tnw_tse"]->rank);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tne_tse"]->rank);

  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bse_tsw"]->rank);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bne_tsw"]->rank);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tse_tsw"]->rank);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tne_tsw"]->rank);

  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["tsw_bne"]->rank);
  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["tse_bne"]->rank);

  CHECK_EQ(domain["tse_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["tse_bnw"]->rank);

  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["tsw_bse"]->rank);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["tse_bse"]->rank);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["tnw_bse"]->rank);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["tne_bse"]->rank);

  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["tse_bsw"]->rank);
  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["tne_bsw"]->rank);

  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bsw_tne"]->rank);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bnw_tne"]->rank);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tsw_tne"]->rank);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tnw_tne"]->rank);

  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bsw_tnw"]->rank);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bse_tnw"]->rank);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bnw_tnw"]->rank);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bne_tnw"]->rank);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tsw_tnw"]->rank);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tse_tnw"]->rank);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tnw_tnw"]->rank);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tne_tnw"]->rank);

  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bnw_tse"]->rank);
  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tnw_tse"]->rank);

  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bnw_tsw"]->rank);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bne_tsw"]->rank);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tnw_tsw"]->rank);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tne_tsw"]->rank);

  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["tsw_bne"]->rank);
  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["tnw_bne"]->rank);

  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["tsw_bnw"]->rank);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["tse_bnw"]->rank);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["tnw_bnw"]->rank);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["tne_bnw"]->rank);

  CHECK_EQ(domain["tnw_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["tnw_bse"]->rank);

  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["tnw_bsw"]->rank);
  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["tne_bsw"]->rank);

  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bsw_tne"]->rank);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bse_tne"]->rank);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bnw_tne"]->rank);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["bne_tne"]->rank);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tsw_tne"]->rank);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tse_tne"]->rank);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tnw_tne"]->rank);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Corner<3>::tne()).rank, domain["tne_tne"]->rank);

  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bse_tnw"]->rank);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["bne_tnw"]->rank);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tse_tnw"]->rank);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Corner<3>::tnw()).rank, domain["tne_tnw"]->rank);

  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bnw_tse"]->rank);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["bne_tse"]->rank);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tnw_tse"]->rank);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Corner<3>::tse()).rank, domain["tne_tse"]->rank);

  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["bne_tsw"]->rank);
  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Corner<3>::tsw()).rank, domain["tne_tsw"]->rank);

  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["tsw_bne"]->rank);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["tse_bne"]->rank);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["tnw_bne"]->rank);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Corner<3>::bne()).rank, domain["tne_bne"]->rank);

  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["tse_bnw"]->rank);
  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Corner<3>::bnw()).rank, domain["tne_bnw"]->rank);

  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["tnw_bse"]->rank);
  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Corner<3>::bse()).rank, domain["tne_bse"]->rank);

  CHECK_EQ(domain["tne_tne"]->getNormalNbrInfo(Corner<3>::bsw()).rank, domain["tne_bsw"]->rank);
}
