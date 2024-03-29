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
Check4x4x4DomainSideNeighborIds(const PatchVector& domain)
{
  // side nbr id
  CHECK_EQ(domain["bsw_bsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bsw_bse"]->id);
  CHECK_EQ(domain["bsw_bsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["bsw_bnw"]->id);
  CHECK_EQ(domain["bsw_bsw"]->getNormalNbrInfo(Side<3>::top()).id, domain["bsw_tsw"]->id);

  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Side<3>::west()).id, domain["bsw_bsw"]->id);
  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Side<3>::east()).id, domain["bse_bsw"]->id);
  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Side<3>::north()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bsw_bse"]->getNormalNbrInfo(Side<3>::top()).id, domain["bsw_tse"]->id);

  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["bsw_bsw"]->id);
  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Side<3>::north()).id, domain["bnw_bsw"]->id);
  CHECK_EQ(domain["bsw_bnw"]->getNormalNbrInfo(Side<3>::top()).id, domain["bsw_tnw"]->id);

  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::west()).id, domain["bsw_bnw"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::east()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::south()).id, domain["bsw_bse"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::north()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::top()).id, domain["bsw_tne"]->id);

  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bsw_bsw"]->id);
  CHECK_EQ(domain["bsw_tsw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tsw_bsw"]->id);

  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::west()).id, domain["bsw_tsw"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::east()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::north()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bsw_bse"]->id);
  CHECK_EQ(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::top()).id, domain["tsw_bse"]->id);

  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["bsw_tsw"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::north()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bsw_bnw"]->id);
  CHECK_EQ(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tsw_bnw"]->id);

  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::west()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::east()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::south()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::north()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::top()).id, domain["tsw_bne"]->id);

  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Side<3>::west()).id, domain["bsw_bse"]->id);
  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bse_bse"]->id);
  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bse_bsw"]->getNormalNbrInfo(Side<3>::top()).id, domain["bse_tsw"]->id);

  CHECK_EQ(domain["bse_bse"]->getNormalNbrInfo(Side<3>::west()).id, domain["bse_bsw"]->id);
  CHECK_EQ(domain["bse_bse"]->getNormalNbrInfo(Side<3>::north()).id, domain["bse_bne"]->id);
  CHECK_EQ(domain["bse_bse"]->getNormalNbrInfo(Side<3>::top()).id, domain["bse_tse"]->id);

  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::west()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bse_bne"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["bse_bsw"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::north()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::top()).id, domain["bse_tnw"]->id);

  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Side<3>::west()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Side<3>::south()).id, domain["bse_bse"]->id);
  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Side<3>::north()).id, domain["bne_bse"]->id);
  CHECK_EQ(domain["bse_bne"]->getNormalNbrInfo(Side<3>::top()).id, domain["bse_tne"]->id);

  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::west()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bse_tse"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bse_bsw"]->id);
  CHECK_EQ(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tse_bsw"]->id);

  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Side<3>::west()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Side<3>::north()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bse_bse"]->id);
  CHECK_EQ(domain["bse_tse"]->getNormalNbrInfo(Side<3>::top()).id, domain["tse_bse"]->id);

  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::west()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::north()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tse_bnw"]->id);

  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Side<3>::west()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Side<3>::south()).id, domain["bse_tse"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Side<3>::north()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bse_bne"]->id);
  CHECK_EQ(domain["bse_tne"]->getNormalNbrInfo(Side<3>::top()).id, domain["tse_bne"]->id);

  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Side<3>::south()).id, domain["bsw_bnw"]->id);
  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["bnw_bnw"]->id);
  CHECK_EQ(domain["bnw_bsw"]->getNormalNbrInfo(Side<3>::top()).id, domain["bnw_tsw"]->id);

  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::west()).id, domain["bnw_bsw"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::east()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::south()).id, domain["bsw_bne"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::north()).id, domain["bnw_bne"]->id);
  CHECK_EQ(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::top()).id, domain["bnw_tse"]->id);

  CHECK_EQ(domain["bnw_bnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bnw_bne"]->id);
  CHECK_EQ(domain["bnw_bnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["bnw_bsw"]->id);
  CHECK_EQ(domain["bnw_bnw"]->getNormalNbrInfo(Side<3>::top()).id, domain["bnw_tnw"]->id);

  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Side<3>::west()).id, domain["bnw_bnw"]->id);
  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Side<3>::east()).id, domain["bne_bnw"]->id);
  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Side<3>::south()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bnw_bne"]->getNormalNbrInfo(Side<3>::top()).id, domain["bnw_tne"]->id);

  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::south()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["bnw_tnw"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bnw_bsw"]->id);
  CHECK_EQ(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tnw_bsw"]->id);

  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::west()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::east()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::south()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::north()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::top()).id, domain["tnw_bse"]->id);

  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bnw_bnw"]->id);
  CHECK_EQ(domain["bnw_tnw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tnw_bnw"]->id);

  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::west()).id, domain["bnw_tnw"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::east()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::south()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bnw_bne"]->id);
  CHECK_EQ(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::top()).id, domain["tnw_bne"]->id);

  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::west()).id, domain["bnw_bse"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bne_bse"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::south()).id, domain["bse_bnw"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["bne_bnw"]->id);
  CHECK_EQ(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::top()).id, domain["bne_tsw"]->id);

  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Side<3>::west()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Side<3>::south()).id, domain["bse_bne"]->id);
  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Side<3>::north()).id, domain["bne_bne"]->id);
  CHECK_EQ(domain["bne_bse"]->getNormalNbrInfo(Side<3>::top()).id, domain["bne_tse"]->id);

  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Side<3>::west()).id, domain["bnw_bne"]->id);
  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bne_bne"]->id);
  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bne_bnw"]->getNormalNbrInfo(Side<3>::top()).id, domain["bne_tnw"]->id);

  CHECK_EQ(domain["bne_bne"]->getNormalNbrInfo(Side<3>::west()).id, domain["bne_bnw"]->id);
  CHECK_EQ(domain["bne_bne"]->getNormalNbrInfo(Side<3>::south()).id, domain["bne_bse"]->id);
  CHECK_EQ(domain["bne_bne"]->getNormalNbrInfo(Side<3>::top()).id, domain["bne_tne"]->id);

  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::west()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::south()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bne_bsw"]->id);
  CHECK_EQ(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tne_bsw"]->id);

  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Side<3>::west()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Side<3>::south()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Side<3>::north()).id, domain["bne_tne"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bne_bse"]->id);
  CHECK_EQ(domain["bne_tse"]->getNormalNbrInfo(Side<3>::top()).id, domain["tne_bse"]->id);

  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::west()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["bne_tne"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bne_bnw"]->id);
  CHECK_EQ(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tne_bnw"]->id);

  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Side<3>::west()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Side<3>::south()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bne_bne"]->id);
  CHECK_EQ(domain["bne_tne"]->getNormalNbrInfo(Side<3>::top()).id, domain["tne_bne"]->id);

  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bsw_tsw"]->id);
  CHECK_EQ(domain["tsw_bsw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tsw_tsw"]->id);

  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::west()).id, domain["tsw_bsw"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::east()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::north()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bsw_tse"]->id);
  CHECK_EQ(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::top()).id, domain["tsw_tse"]->id);

  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["tsw_bsw"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::north()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bsw_tnw"]->id);
  CHECK_EQ(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tsw_tnw"]->id);

  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::west()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::east()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::south()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::north()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bsw_tne"]->id);
  CHECK_EQ(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::top()).id, domain["tsw_tne"]->id);

  CHECK_EQ(domain["tsw_tsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tsw_tse"]->id);
  CHECK_EQ(domain["tsw_tsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["tsw_tnw"]->id);
  CHECK_EQ(domain["tsw_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tsw_bsw"]->id);

  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Side<3>::west()).id, domain["tsw_tsw"]->id);
  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Side<3>::east()).id, domain["tse_tsw"]->id);
  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Side<3>::north()).id, domain["tsw_tne"]->id);
  CHECK_EQ(domain["tsw_tse"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tsw_bse"]->id);

  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tsw_tne"]->id);
  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["tsw_tsw"]->id);
  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Side<3>::north()).id, domain["tnw_tsw"]->id);
  CHECK_EQ(domain["tsw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tsw_bnw"]->id);

  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::west()).id, domain["tsw_tnw"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::east()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::south()).id, domain["tsw_tse"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::north()).id, domain["tnw_tse"]->id);
  CHECK_EQ(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tsw_bne"]->id);

  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::west()).id, domain["tsw_bse"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tse_bse"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bse_tsw"]->id);
  CHECK_EQ(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tse_tsw"]->id);

  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Side<3>::west()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Side<3>::north()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bse_tse"]->id);
  CHECK_EQ(domain["tse_bse"]->getNormalNbrInfo(Side<3>::top()).id, domain["tse_tse"]->id);

  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::west()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["tse_bsw"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::north()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bse_tnw"]->id);
  CHECK_EQ(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tse_tnw"]->id);

  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Side<3>::west()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Side<3>::south()).id, domain["tse_bse"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Side<3>::north()).id, domain["tne_bse"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bse_tne"]->id);
  CHECK_EQ(domain["tse_bne"]->getNormalNbrInfo(Side<3>::top()).id, domain["tse_tne"]->id);

  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Side<3>::west()).id, domain["tsw_tse"]->id);
  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tse_tse"]->id);
  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tse_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tse_bsw"]->id);

  CHECK_EQ(domain["tse_tse"]->getNormalNbrInfo(Side<3>::west()).id, domain["tse_tsw"]->id);
  CHECK_EQ(domain["tse_tse"]->getNormalNbrInfo(Side<3>::north()).id, domain["tse_tne"]->id);
  CHECK_EQ(domain["tse_tse"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tse_bse"]->id);

  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::west()).id, domain["tsw_tne"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tse_tne"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["tse_tsw"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::north()).id, domain["tne_tsw"]->id);
  CHECK_EQ(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tse_bnw"]->id);

  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Side<3>::west()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Side<3>::south()).id, domain["tse_tse"]->id);
  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Side<3>::north()).id, domain["tne_tse"]->id);
  CHECK_EQ(domain["tse_tne"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tse_bne"]->id);

  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::south()).id, domain["tsw_bnw"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["tnw_bnw"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bnw_tsw"]->id);
  CHECK_EQ(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tnw_tsw"]->id);

  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::west()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::east()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::south()).id, domain["tsw_bne"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::north()).id, domain["tnw_bne"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bnw_tse"]->id);
  CHECK_EQ(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::top()).id, domain["tnw_tse"]->id);

  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tnw_bne"]->id);
  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["tnw_bsw"]->id);
  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bnw_tnw"]->id);
  CHECK_EQ(domain["tnw_bnw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tnw_tnw"]->id);

  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::west()).id, domain["tnw_bnw"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::east()).id, domain["tne_bnw"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::south()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bnw_tne"]->id);
  CHECK_EQ(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::top()).id, domain["tnw_tne"]->id);

  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tnw_tse"]->id);
  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Side<3>::south()).id, domain["tsw_tnw"]->id);
  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["tnw_tnw"]->id);
  CHECK_EQ(domain["tnw_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tnw_bsw"]->id);

  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::west()).id, domain["tnw_tsw"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::east()).id, domain["tne_tsw"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::south()).id, domain["tsw_tne"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::north()).id, domain["tnw_tne"]->id);
  CHECK_EQ(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tnw_bse"]->id);

  CHECK_EQ(domain["tnw_tnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tnw_tne"]->id);
  CHECK_EQ(domain["tnw_tnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["tnw_tsw"]->id);
  CHECK_EQ(domain["tnw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tnw_bnw"]->id);

  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Side<3>::west()).id, domain["tnw_tnw"]->id);
  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Side<3>::east()).id, domain["tne_tnw"]->id);
  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Side<3>::south()).id, domain["tnw_tse"]->id);
  CHECK_EQ(domain["tnw_tne"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tnw_bne"]->id);

  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::west()).id, domain["tnw_bse"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tne_bse"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::south()).id, domain["tse_bnw"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["tne_bnw"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bne_tsw"]->id);
  CHECK_EQ(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tne_tsw"]->id);

  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Side<3>::west()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Side<3>::south()).id, domain["tse_bne"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Side<3>::north()).id, domain["tne_bne"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bne_tse"]->id);
  CHECK_EQ(domain["tne_bse"]->getNormalNbrInfo(Side<3>::top()).id, domain["tne_tse"]->id);

  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::west()).id, domain["tnw_bne"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tne_bne"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["tne_bsw"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bne_tnw"]->id);
  CHECK_EQ(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::top()).id, domain["tne_tnw"]->id);

  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Side<3>::west()).id, domain["tne_bnw"]->id);
  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Side<3>::south()).id, domain["tne_bse"]->id);
  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["bne_tne"]->id);
  CHECK_EQ(domain["tne_bne"]->getNormalNbrInfo(Side<3>::top()).id, domain["tne_tne"]->id);

  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::west()).id, domain["tnw_tse"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tne_tse"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::south()).id, domain["tse_tnw"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::north()).id, domain["tne_tnw"]->id);
  CHECK_EQ(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tne_bsw"]->id);

  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Side<3>::west()).id, domain["tne_tsw"]->id);
  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Side<3>::south()).id, domain["tse_tne"]->id);
  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Side<3>::north()).id, domain["tne_tne"]->id);
  CHECK_EQ(domain["tne_tse"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tne_bse"]->id);

  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Side<3>::west()).id, domain["tnw_tne"]->id);
  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Side<3>::east()).id, domain["tne_tne"]->id);
  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Side<3>::south()).id, domain["tne_tsw"]->id);
  CHECK_EQ(domain["tne_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tne_bnw"]->id);

  CHECK_EQ(domain["tne_tne"]->getNormalNbrInfo(Side<3>::west()).id, domain["tne_tnw"]->id);
  CHECK_EQ(domain["tne_tne"]->getNormalNbrInfo(Side<3>::south()).id, domain["tne_tse"]->id);
  CHECK_EQ(domain["tne_tne"]->getNormalNbrInfo(Side<3>::bottom()).id, domain["tne_bne"]->id);
}
