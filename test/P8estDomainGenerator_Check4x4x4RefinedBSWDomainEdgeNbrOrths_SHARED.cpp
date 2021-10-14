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
Check4x4x4RefinedBSWDomainEdgeNeighborOrths(const PatchVector& domain)
{
  CHECK(domain["bsw_bsw_tnw"]->getCoarseNbrInfo(Edge::tn()).orth_on_coarse == Orthant<1>::lower());
  CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Edge::tn()).orth_on_coarse == Orthant<1>::upper());

  CHECK(domain["bsw_bsw_tse"]->getCoarseNbrInfo(Edge::te()).orth_on_coarse == Orthant<1>::lower());
  CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Edge::te()).orth_on_coarse == Orthant<1>::upper());

  CHECK(domain["bsw_bsw_bne"]->getCoarseNbrInfo(Edge::ne()).orth_on_coarse == Orthant<1>::lower());
  CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Edge::ne()).orth_on_coarse == Orthant<1>::upper());
}
