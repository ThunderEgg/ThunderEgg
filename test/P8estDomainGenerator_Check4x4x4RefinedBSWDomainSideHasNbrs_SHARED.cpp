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
Check4x4x4RefinedBSWDomainSideHasNeighbors(const PatchVector& domain)
{
  CHECK(domain["bsw_bsw_bsw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["bsw_bsw_bsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_bsw_bsw"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["bsw_bsw_bsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_bsw_bsw"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bsw_bsw_bsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_bsw_bse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bsw_bsw_bse"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_bsw_bse"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["bsw_bsw_bse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_bsw_bse"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bsw_bsw_bse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_bsw_bnw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["bsw_bsw_bnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_bsw_bnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bsw_bsw_bnw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_bsw_bnw"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bsw_bsw_bnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_bsw_bne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bsw_bsw_bne"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_bsw_bne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bsw_bsw_bne"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_bsw_bne"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bsw_bsw_bne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_bsw_tsw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["bsw_bsw_tsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_bsw_tsw"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["bsw_bsw_tsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_bsw_tsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bsw_bsw_tsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_bsw_tse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bsw_bsw_tse"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_bsw_tse"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["bsw_bsw_tse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_bsw_tse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bsw_bsw_tse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_bsw_tnw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["bsw_bsw_tnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_bsw_tnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bsw_bsw_tnw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_bsw_tnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bsw_bsw_tnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_bsw_tne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bsw_bsw_tne"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_bsw_tne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bsw_bsw_tne"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_bsw_tne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bsw_bsw_tne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_bse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bsw_bse"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_bse"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["bsw_bse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_bse"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bsw_bse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_bnw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["bsw_bnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_bnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bsw_bnw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_bnw"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bsw_bnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_bne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bsw_bne"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_bne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bsw_bne"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_bne"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bsw_bne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_tsw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["bsw_tsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_tsw"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["bsw_tsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_tsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bsw_tsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_tse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bsw_tse"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_tse"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["bsw_tse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_tse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bsw_tse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_tnw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["bsw_tnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_tnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bsw_tnw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_tnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bsw_tnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bsw_tne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bsw_tne"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bsw_tne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bsw_tne"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bsw_tne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bsw_tne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bse_bsw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bse_bsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bse_bsw"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["bse_bsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bse_bsw"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bse_bsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bse_bse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bse_bse"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["bse_bse"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["bse_bse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bse_bse"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bse_bse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bse_bnw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bse_bnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bse_bnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bse_bnw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bse_bnw"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bse_bnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bse_bne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bse_bne"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["bse_bne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bse_bne"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bse_bne"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bse_bne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bse_tsw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bse_tsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bse_tsw"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["bse_tsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bse_tsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bse_tsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bse_tse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bse_tse"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["bse_tse"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["bse_tse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bse_tse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bse_tse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bse_tnw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bse_tnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bse_tnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bse_tnw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bse_tnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bse_tnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bse_tne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bse_tne"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["bse_tne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bse_tne"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bse_tne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bse_tne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bnw_bsw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["bnw_bsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bnw_bsw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bnw_bsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bnw_bsw"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bnw_bsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bnw_bse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bnw_bse"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bnw_bse"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bnw_bse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bnw_bse"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bnw_bse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bnw_bnw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["bnw_bnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bnw_bnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bnw_bnw"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["bnw_bnw"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bnw_bnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bnw_bne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bnw_bne"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bnw_bne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bnw_bne"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["bnw_bne"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bnw_bne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bnw_tsw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["bnw_tsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bnw_tsw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bnw_tsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bnw_tsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bnw_tsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bnw_tse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bnw_tse"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bnw_tse"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bnw_tse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bnw_tse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bnw_tse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bnw_tnw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["bnw_tnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bnw_tnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bnw_tnw"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["bnw_tnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bnw_tnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bnw_tne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bnw_tne"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bnw_tne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bnw_tne"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["bnw_tne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bnw_tne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bne_bsw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bne_bsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bne_bsw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bne_bsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bne_bsw"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bne_bsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bne_bse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bne_bse"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["bne_bse"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bne_bse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bne_bse"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bne_bse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bne_bnw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bne_bnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bne_bnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bne_bnw"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["bne_bnw"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bne_bnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bne_bne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bne_bne"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["bne_bne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bne_bne"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["bne_bne"]->hasNbr(Side<3>::bottom()) == false);
  CHECK(domain["bne_bne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bne_tsw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bne_tsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bne_tsw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bne_tsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bne_tsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bne_tsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bne_tse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bne_tse"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["bne_tse"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bne_tse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["bne_tse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bne_tse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bne_tnw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bne_tnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["bne_tnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bne_tnw"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["bne_tnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bne_tnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["bne_tne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["bne_tne"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["bne_tne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["bne_tne"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["bne_tne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["bne_tne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tsw_bsw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["tsw_bsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tsw_bsw"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["tsw_bsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tsw_bsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tsw_bsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tsw_bse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tsw_bse"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tsw_bse"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["tsw_bse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tsw_bse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tsw_bse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tsw_bnw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["tsw_bnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tsw_bnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tsw_bnw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tsw_bnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tsw_bnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tsw_bne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tsw_bne"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tsw_bne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tsw_bne"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tsw_bne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tsw_bne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tsw_tsw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["tsw_tsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tsw_tsw"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["tsw_tsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tsw_tsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tsw_tsw"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tsw_tse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tsw_tse"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tsw_tse"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["tsw_tse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tsw_tse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tsw_tse"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tsw_tnw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["tsw_tnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tsw_tnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tsw_tnw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tsw_tnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tsw_tnw"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tsw_tne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tsw_tne"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tsw_tne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tsw_tne"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tsw_tne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tsw_tne"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tse_bsw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tse_bsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tse_bsw"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["tse_bsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tse_bsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tse_bsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tse_bse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tse_bse"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["tse_bse"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["tse_bse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tse_bse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tse_bse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tse_bnw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tse_bnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tse_bnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tse_bnw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tse_bnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tse_bnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tse_bne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tse_bne"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["tse_bne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tse_bne"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tse_bne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tse_bne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tse_tsw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tse_tsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tse_tsw"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["tse_tsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tse_tsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tse_tsw"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tse_tse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tse_tse"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["tse_tse"]->hasNbr(Side<3>::south()) == false);
  CHECK(domain["tse_tse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tse_tse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tse_tse"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tse_tnw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tse_tnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tse_tnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tse_tnw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tse_tnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tse_tnw"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tse_tne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tse_tne"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["tse_tne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tse_tne"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tse_tne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tse_tne"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tnw_bsw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["tnw_bsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tnw_bsw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tnw_bsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tnw_bsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tnw_bsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tnw_bse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tnw_bse"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tnw_bse"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tnw_bse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tnw_bse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tnw_bse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tnw_bnw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["tnw_bnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tnw_bnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tnw_bnw"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["tnw_bnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tnw_bnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tnw_bne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tnw_bne"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tnw_bne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tnw_bne"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["tnw_bne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tnw_bne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tnw_tsw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["tnw_tsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tnw_tsw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tnw_tsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tnw_tsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tnw_tsw"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tnw_tse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tnw_tse"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tnw_tse"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tnw_tse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tnw_tse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tnw_tse"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tnw_tnw"]->hasNbr(Side<3>::west()) == false);
  CHECK(domain["tnw_tnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tnw_tnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tnw_tnw"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["tnw_tnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tnw_tnw"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tnw_tne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tnw_tne"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tnw_tne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tnw_tne"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["tnw_tne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tnw_tne"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tne_bsw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tne_bsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tne_bsw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tne_bsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tne_bsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tne_bsw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tne_bse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tne_bse"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["tne_bse"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tne_bse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tne_bse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tne_bse"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tne_bnw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tne_bnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tne_bnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tne_bnw"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["tne_bnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tne_bnw"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tne_bne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tne_bne"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["tne_bne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tne_bne"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["tne_bne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tne_bne"]->hasNbr(Side<3>::top()) == true);

  CHECK(domain["tne_tsw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tne_tsw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tne_tsw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tne_tsw"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tne_tsw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tne_tsw"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tne_tse"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tne_tse"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["tne_tse"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tne_tse"]->hasNbr(Side<3>::north()) == true);
  CHECK(domain["tne_tse"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tne_tse"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tne_tnw"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tne_tnw"]->hasNbr(Side<3>::east()) == true);
  CHECK(domain["tne_tnw"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tne_tnw"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["tne_tnw"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tne_tnw"]->hasNbr(Side<3>::top()) == false);

  CHECK(domain["tne_tne"]->hasNbr(Side<3>::west()) == true);
  CHECK(domain["tne_tne"]->hasNbr(Side<3>::east()) == false);
  CHECK(domain["tne_tne"]->hasNbr(Side<3>::south()) == true);
  CHECK(domain["tne_tne"]->hasNbr(Side<3>::north()) == false);
  CHECK(domain["tne_tne"]->hasNbr(Side<3>::bottom()) == true);
  CHECK(domain["tne_tne"]->hasNbr(Side<3>::top()) == false);
}
