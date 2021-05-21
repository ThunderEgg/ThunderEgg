/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Gebneral Public Licenbse as published by
 *  the Free Software Foundation, either version 3 of the Licenbse, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be ubseful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Gebneral Public Licenbse for more details.
 *
 *  You should have received a copy of the GNU Gebneral Public Licenbse
 *  along with this program.  If not, bsee <https://www.gnu.org/licenbses/>.
 ***************************************************************************/

#include "P8estDomainGenerator_SHARED.h"
#include <catch2/catch_approx.hpp>
using namespace std;
using namespace ThunderEgg;

#include <catch2/catch_test_macros.hpp>

void Check4x4x4RefinedBSWDomainSideNeighborOrths(const PatchVector &domain)
{
	CHECK(domain["bsw_bsw_bse"]->getCoarseNbrInfo(Side<3>::east()).orth_on_coarse == Orthant<2>::sw());
	CHECK(domain["bsw_bsw_bne"]->getCoarseNbrInfo(Side<3>::east()).orth_on_coarse == Orthant<2>::se());
	CHECK(domain["bsw_bsw_tse"]->getCoarseNbrInfo(Side<3>::east()).orth_on_coarse == Orthant<2>::nw());
	CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Side<3>::east()).orth_on_coarse == Orthant<2>::ne());

	CHECK(domain["bsw_bsw_bnw"]->getCoarseNbrInfo(Side<3>::north()).orth_on_coarse == Orthant<2>::sw());
	CHECK(domain["bsw_bsw_bne"]->getCoarseNbrInfo(Side<3>::north()).orth_on_coarse == Orthant<2>::se());
	CHECK(domain["bsw_bsw_tnw"]->getCoarseNbrInfo(Side<3>::north()).orth_on_coarse == Orthant<2>::nw());
	CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Side<3>::north()).orth_on_coarse == Orthant<2>::ne());

	CHECK(domain["bsw_bsw_tsw"]->getCoarseNbrInfo(Side<3>::top()).orth_on_coarse == Orthant<2>::sw());
	CHECK(domain["bsw_bsw_tse"]->getCoarseNbrInfo(Side<3>::top()).orth_on_coarse == Orthant<2>::se());
	CHECK(domain["bsw_bsw_tnw"]->getCoarseNbrInfo(Side<3>::top()).orth_on_coarse == Orthant<2>::nw());
	CHECK(domain["bsw_bsw_tne"]->getCoarseNbrInfo(Side<3>::top()).orth_on_coarse == Orthant<2>::ne());
}
