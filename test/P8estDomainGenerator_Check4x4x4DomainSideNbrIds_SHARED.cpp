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

void Check4x4x4DomainSideNeighborIds(const PatchVector &domain)
{
	//side nbr id
	CHECK(domain["bsw_bsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bsw_bse"]->id);
	CHECK(domain["bsw_bsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["bsw_bnw"]->id);
	CHECK(domain["bsw_bsw"]->getNormalNbrInfo(Side<3>::top()).id == domain["bsw_tsw"]->id);

	CHECK(domain["bsw_bse"]->getNormalNbrInfo(Side<3>::west()).id == domain["bsw_bsw"]->id);
	CHECK(domain["bsw_bse"]->getNormalNbrInfo(Side<3>::east()).id == domain["bse_bsw"]->id);
	CHECK(domain["bsw_bse"]->getNormalNbrInfo(Side<3>::north()).id == domain["bsw_bne"]->id);
	CHECK(domain["bsw_bse"]->getNormalNbrInfo(Side<3>::top()).id == domain["bsw_tse"]->id);

	CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bsw_bne"]->id);
	CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["bsw_bsw"]->id);
	CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Side<3>::north()).id == domain["bnw_bsw"]->id);
	CHECK(domain["bsw_bnw"]->getNormalNbrInfo(Side<3>::top()).id == domain["bsw_tnw"]->id);

	CHECK(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::west()).id == domain["bsw_bnw"]->id);
	CHECK(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::east()).id == domain["bse_bnw"]->id);
	CHECK(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::south()).id == domain["bsw_bse"]->id);
	CHECK(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::north()).id == domain["bnw_bse"]->id);
	CHECK(domain["bsw_bne"]->getNormalNbrInfo(Side<3>::top()).id == domain["bsw_tne"]->id);

	CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bsw_tse"]->id);
	CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["bsw_tnw"]->id);
	CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bsw_bsw"]->id);
	CHECK(domain["bsw_tsw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tsw_bsw"]->id);

	CHECK(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::west()).id == domain["bsw_tsw"]->id);
	CHECK(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::east()).id == domain["bse_tsw"]->id);
	CHECK(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::north()).id == domain["bsw_tne"]->id);
	CHECK(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bsw_bse"]->id);
	CHECK(domain["bsw_tse"]->getNormalNbrInfo(Side<3>::top()).id == domain["tsw_bse"]->id);

	CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bsw_tne"]->id);
	CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["bsw_tsw"]->id);
	CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::north()).id == domain["bnw_tsw"]->id);
	CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bsw_bnw"]->id);
	CHECK(domain["bsw_tnw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tsw_bnw"]->id);

	CHECK(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::west()).id == domain["bsw_tnw"]->id);
	CHECK(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::east()).id == domain["bse_tnw"]->id);
	CHECK(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::south()).id == domain["bsw_tse"]->id);
	CHECK(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::north()).id == domain["bnw_tse"]->id);
	CHECK(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bsw_bne"]->id);
	CHECK(domain["bsw_tne"]->getNormalNbrInfo(Side<3>::top()).id == domain["tsw_bne"]->id);

	CHECK(domain["bse_bsw"]->getNormalNbrInfo(Side<3>::west()).id == domain["bsw_bse"]->id);
	CHECK(domain["bse_bsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bse_bse"]->id);
	CHECK(domain["bse_bsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["bse_bnw"]->id);
	CHECK(domain["bse_bsw"]->getNormalNbrInfo(Side<3>::top()).id == domain["bse_tsw"]->id);

	CHECK(domain["bse_bse"]->getNormalNbrInfo(Side<3>::west()).id == domain["bse_bsw"]->id);
	CHECK(domain["bse_bse"]->getNormalNbrInfo(Side<3>::north()).id == domain["bse_bne"]->id);
	CHECK(domain["bse_bse"]->getNormalNbrInfo(Side<3>::top()).id == domain["bse_tse"]->id);

	CHECK(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::west()).id == domain["bsw_bne"]->id);
	CHECK(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bse_bne"]->id);
	CHECK(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["bse_bsw"]->id);
	CHECK(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::north()).id == domain["bne_bsw"]->id);
	CHECK(domain["bse_bnw"]->getNormalNbrInfo(Side<3>::top()).id == domain["bse_tnw"]->id);

	CHECK(domain["bse_bne"]->getNormalNbrInfo(Side<3>::west()).id == domain["bse_bnw"]->id);
	CHECK(domain["bse_bne"]->getNormalNbrInfo(Side<3>::south()).id == domain["bse_bse"]->id);
	CHECK(domain["bse_bne"]->getNormalNbrInfo(Side<3>::north()).id == domain["bne_bse"]->id);
	CHECK(domain["bse_bne"]->getNormalNbrInfo(Side<3>::top()).id == domain["bse_tne"]->id);

	CHECK(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::west()).id == domain["bsw_tse"]->id);
	CHECK(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bse_tse"]->id);
	CHECK(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["bse_tnw"]->id);
	CHECK(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bse_bsw"]->id);
	CHECK(domain["bse_tsw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tse_bsw"]->id);

	CHECK(domain["bse_tse"]->getNormalNbrInfo(Side<3>::west()).id == domain["bse_tsw"]->id);
	CHECK(domain["bse_tse"]->getNormalNbrInfo(Side<3>::north()).id == domain["bse_tne"]->id);
	CHECK(domain["bse_tse"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bse_bse"]->id);
	CHECK(domain["bse_tse"]->getNormalNbrInfo(Side<3>::top()).id == domain["tse_bse"]->id);

	CHECK(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::west()).id == domain["bsw_tne"]->id);
	CHECK(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bse_tne"]->id);
	CHECK(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["bse_tsw"]->id);
	CHECK(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::north()).id == domain["bne_tsw"]->id);
	CHECK(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bse_bnw"]->id);
	CHECK(domain["bse_tnw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tse_bnw"]->id);

	CHECK(domain["bse_tne"]->getNormalNbrInfo(Side<3>::west()).id == domain["bse_tnw"]->id);
	CHECK(domain["bse_tne"]->getNormalNbrInfo(Side<3>::south()).id == domain["bse_tse"]->id);
	CHECK(domain["bse_tne"]->getNormalNbrInfo(Side<3>::north()).id == domain["bne_tse"]->id);
	CHECK(domain["bse_tne"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bse_bne"]->id);
	CHECK(domain["bse_tne"]->getNormalNbrInfo(Side<3>::top()).id == domain["tse_bne"]->id);

	CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bnw_bse"]->id);
	CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Side<3>::south()).id == domain["bsw_bnw"]->id);
	CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["bnw_bnw"]->id);
	CHECK(domain["bnw_bsw"]->getNormalNbrInfo(Side<3>::top()).id == domain["bnw_tsw"]->id);

	CHECK(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::west()).id == domain["bnw_bsw"]->id);
	CHECK(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::east()).id == domain["bne_bsw"]->id);
	CHECK(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::south()).id == domain["bsw_bne"]->id);
	CHECK(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::north()).id == domain["bnw_bne"]->id);
	CHECK(domain["bnw_bse"]->getNormalNbrInfo(Side<3>::top()).id == domain["bnw_tse"]->id);

	CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bnw_bne"]->id);
	CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["bnw_bsw"]->id);
	CHECK(domain["bnw_bnw"]->getNormalNbrInfo(Side<3>::top()).id == domain["bnw_tnw"]->id);

	CHECK(domain["bnw_bne"]->getNormalNbrInfo(Side<3>::west()).id == domain["bnw_bnw"]->id);
	CHECK(domain["bnw_bne"]->getNormalNbrInfo(Side<3>::east()).id == domain["bne_bnw"]->id);
	CHECK(domain["bnw_bne"]->getNormalNbrInfo(Side<3>::south()).id == domain["bnw_bse"]->id);
	CHECK(domain["bnw_bne"]->getNormalNbrInfo(Side<3>::top()).id == domain["bnw_tne"]->id);

	CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bnw_tse"]->id);
	CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::south()).id == domain["bsw_tnw"]->id);
	CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["bnw_tnw"]->id);
	CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bnw_bsw"]->id);
	CHECK(domain["bnw_tsw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tnw_bsw"]->id);

	CHECK(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::west()).id == domain["bnw_tsw"]->id);
	CHECK(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::east()).id == domain["bne_tsw"]->id);
	CHECK(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::south()).id == domain["bsw_tne"]->id);
	CHECK(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::north()).id == domain["bnw_tne"]->id);
	CHECK(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bnw_bse"]->id);
	CHECK(domain["bnw_tse"]->getNormalNbrInfo(Side<3>::top()).id == domain["tnw_bse"]->id);

	CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bnw_tne"]->id);
	CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["bnw_tsw"]->id);
	CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bnw_bnw"]->id);
	CHECK(domain["bnw_tnw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tnw_bnw"]->id);

	CHECK(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::west()).id == domain["bnw_tnw"]->id);
	CHECK(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::east()).id == domain["bne_tnw"]->id);
	CHECK(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::south()).id == domain["bnw_tse"]->id);
	CHECK(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bnw_bne"]->id);
	CHECK(domain["bnw_tne"]->getNormalNbrInfo(Side<3>::top()).id == domain["tnw_bne"]->id);

	CHECK(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::west()).id == domain["bnw_bse"]->id);
	CHECK(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bne_bse"]->id);
	CHECK(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::south()).id == domain["bse_bnw"]->id);
	CHECK(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["bne_bnw"]->id);
	CHECK(domain["bne_bsw"]->getNormalNbrInfo(Side<3>::top()).id == domain["bne_tsw"]->id);

	CHECK(domain["bne_bse"]->getNormalNbrInfo(Side<3>::west()).id == domain["bne_bsw"]->id);
	CHECK(domain["bne_bse"]->getNormalNbrInfo(Side<3>::south()).id == domain["bse_bne"]->id);
	CHECK(domain["bne_bse"]->getNormalNbrInfo(Side<3>::north()).id == domain["bne_bne"]->id);
	CHECK(domain["bne_bse"]->getNormalNbrInfo(Side<3>::top()).id == domain["bne_tse"]->id);

	CHECK(domain["bne_bnw"]->getNormalNbrInfo(Side<3>::west()).id == domain["bnw_bne"]->id);
	CHECK(domain["bne_bnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bne_bne"]->id);
	CHECK(domain["bne_bnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["bne_bsw"]->id);
	CHECK(domain["bne_bnw"]->getNormalNbrInfo(Side<3>::top()).id == domain["bne_tnw"]->id);

	CHECK(domain["bne_bne"]->getNormalNbrInfo(Side<3>::west()).id == domain["bne_bnw"]->id);
	CHECK(domain["bne_bne"]->getNormalNbrInfo(Side<3>::south()).id == domain["bne_bse"]->id);
	CHECK(domain["bne_bne"]->getNormalNbrInfo(Side<3>::top()).id == domain["bne_tne"]->id);

	CHECK(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::west()).id == domain["bnw_tse"]->id);
	CHECK(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bne_tse"]->id);
	CHECK(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::south()).id == domain["bse_tnw"]->id);
	CHECK(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["bne_tnw"]->id);
	CHECK(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bne_bsw"]->id);
	CHECK(domain["bne_tsw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tne_bsw"]->id);

	CHECK(domain["bne_tse"]->getNormalNbrInfo(Side<3>::west()).id == domain["bne_tsw"]->id);
	CHECK(domain["bne_tse"]->getNormalNbrInfo(Side<3>::south()).id == domain["bse_tne"]->id);
	CHECK(domain["bne_tse"]->getNormalNbrInfo(Side<3>::north()).id == domain["bne_tne"]->id);
	CHECK(domain["bne_tse"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bne_bse"]->id);
	CHECK(domain["bne_tse"]->getNormalNbrInfo(Side<3>::top()).id == domain["tne_bse"]->id);

	CHECK(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::west()).id == domain["bnw_tne"]->id);
	CHECK(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["bne_tne"]->id);
	CHECK(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["bne_tsw"]->id);
	CHECK(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bne_bnw"]->id);
	CHECK(domain["bne_tnw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tne_bnw"]->id);

	CHECK(domain["bne_tne"]->getNormalNbrInfo(Side<3>::west()).id == domain["bne_tnw"]->id);
	CHECK(domain["bne_tne"]->getNormalNbrInfo(Side<3>::south()).id == domain["bne_tse"]->id);
	CHECK(domain["bne_tne"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bne_bne"]->id);
	CHECK(domain["bne_tne"]->getNormalNbrInfo(Side<3>::top()).id == domain["tne_bne"]->id);

	CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tsw_bse"]->id);
	CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["tsw_bnw"]->id);
	CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bsw_tsw"]->id);
	CHECK(domain["tsw_bsw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tsw_tsw"]->id);

	CHECK(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::west()).id == domain["tsw_bsw"]->id);
	CHECK(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::east()).id == domain["tse_bsw"]->id);
	CHECK(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::north()).id == domain["tsw_bne"]->id);
	CHECK(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bsw_tse"]->id);
	CHECK(domain["tsw_bse"]->getNormalNbrInfo(Side<3>::top()).id == domain["tsw_tse"]->id);

	CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tsw_bne"]->id);
	CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["tsw_bsw"]->id);
	CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::north()).id == domain["tnw_bsw"]->id);
	CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bsw_tnw"]->id);
	CHECK(domain["tsw_bnw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tsw_tnw"]->id);

	CHECK(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::west()).id == domain["tsw_bnw"]->id);
	CHECK(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::east()).id == domain["tse_bnw"]->id);
	CHECK(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::south()).id == domain["tsw_bse"]->id);
	CHECK(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::north()).id == domain["tnw_bse"]->id);
	CHECK(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bsw_tne"]->id);
	CHECK(domain["tsw_bne"]->getNormalNbrInfo(Side<3>::top()).id == domain["tsw_tne"]->id);

	CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tsw_tse"]->id);
	CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["tsw_tnw"]->id);
	CHECK(domain["tsw_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tsw_bsw"]->id);

	CHECK(domain["tsw_tse"]->getNormalNbrInfo(Side<3>::west()).id == domain["tsw_tsw"]->id);
	CHECK(domain["tsw_tse"]->getNormalNbrInfo(Side<3>::east()).id == domain["tse_tsw"]->id);
	CHECK(domain["tsw_tse"]->getNormalNbrInfo(Side<3>::north()).id == domain["tsw_tne"]->id);
	CHECK(domain["tsw_tse"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tsw_bse"]->id);

	CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tsw_tne"]->id);
	CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["tsw_tsw"]->id);
	CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Side<3>::north()).id == domain["tnw_tsw"]->id);
	CHECK(domain["tsw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tsw_bnw"]->id);

	CHECK(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::west()).id == domain["tsw_tnw"]->id);
	CHECK(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::east()).id == domain["tse_tnw"]->id);
	CHECK(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::south()).id == domain["tsw_tse"]->id);
	CHECK(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::north()).id == domain["tnw_tse"]->id);
	CHECK(domain["tsw_tne"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tsw_bne"]->id);

	CHECK(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::west()).id == domain["tsw_bse"]->id);
	CHECK(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tse_bse"]->id);
	CHECK(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["tse_bnw"]->id);
	CHECK(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bse_tsw"]->id);
	CHECK(domain["tse_bsw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tse_tsw"]->id);

	CHECK(domain["tse_bse"]->getNormalNbrInfo(Side<3>::west()).id == domain["tse_bsw"]->id);
	CHECK(domain["tse_bse"]->getNormalNbrInfo(Side<3>::north()).id == domain["tse_bne"]->id);
	CHECK(domain["tse_bse"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bse_tse"]->id);
	CHECK(domain["tse_bse"]->getNormalNbrInfo(Side<3>::top()).id == domain["tse_tse"]->id);

	CHECK(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::west()).id == domain["tsw_bne"]->id);
	CHECK(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tse_bne"]->id);
	CHECK(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["tse_bsw"]->id);
	CHECK(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::north()).id == domain["tne_bsw"]->id);
	CHECK(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bse_tnw"]->id);
	CHECK(domain["tse_bnw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tse_tnw"]->id);

	CHECK(domain["tse_bne"]->getNormalNbrInfo(Side<3>::west()).id == domain["tse_bnw"]->id);
	CHECK(domain["tse_bne"]->getNormalNbrInfo(Side<3>::south()).id == domain["tse_bse"]->id);
	CHECK(domain["tse_bne"]->getNormalNbrInfo(Side<3>::north()).id == domain["tne_bse"]->id);
	CHECK(domain["tse_bne"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bse_tne"]->id);
	CHECK(domain["tse_bne"]->getNormalNbrInfo(Side<3>::top()).id == domain["tse_tne"]->id);

	CHECK(domain["tse_tsw"]->getNormalNbrInfo(Side<3>::west()).id == domain["tsw_tse"]->id);
	CHECK(domain["tse_tsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tse_tse"]->id);
	CHECK(domain["tse_tsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["tse_tnw"]->id);
	CHECK(domain["tse_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tse_bsw"]->id);

	CHECK(domain["tse_tse"]->getNormalNbrInfo(Side<3>::west()).id == domain["tse_tsw"]->id);
	CHECK(domain["tse_tse"]->getNormalNbrInfo(Side<3>::north()).id == domain["tse_tne"]->id);
	CHECK(domain["tse_tse"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tse_bse"]->id);

	CHECK(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::west()).id == domain["tsw_tne"]->id);
	CHECK(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tse_tne"]->id);
	CHECK(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["tse_tsw"]->id);
	CHECK(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::north()).id == domain["tne_tsw"]->id);
	CHECK(domain["tse_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tse_bnw"]->id);

	CHECK(domain["tse_tne"]->getNormalNbrInfo(Side<3>::west()).id == domain["tse_tnw"]->id);
	CHECK(domain["tse_tne"]->getNormalNbrInfo(Side<3>::south()).id == domain["tse_tse"]->id);
	CHECK(domain["tse_tne"]->getNormalNbrInfo(Side<3>::north()).id == domain["tne_tse"]->id);
	CHECK(domain["tse_tne"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tse_bne"]->id);

	CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tnw_bse"]->id);
	CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::south()).id == domain["tsw_bnw"]->id);
	CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["tnw_bnw"]->id);
	CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bnw_tsw"]->id);
	CHECK(domain["tnw_bsw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tnw_tsw"]->id);

	CHECK(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::west()).id == domain["tnw_bsw"]->id);
	CHECK(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::east()).id == domain["tne_bsw"]->id);
	CHECK(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::south()).id == domain["tsw_bne"]->id);
	CHECK(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::north()).id == domain["tnw_bne"]->id);
	CHECK(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bnw_tse"]->id);
	CHECK(domain["tnw_bse"]->getNormalNbrInfo(Side<3>::top()).id == domain["tnw_tse"]->id);

	CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tnw_bne"]->id);
	CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["tnw_bsw"]->id);
	CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bnw_tnw"]->id);
	CHECK(domain["tnw_bnw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tnw_tnw"]->id);

	CHECK(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::west()).id == domain["tnw_bnw"]->id);
	CHECK(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::east()).id == domain["tne_bnw"]->id);
	CHECK(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::south()).id == domain["tnw_bse"]->id);
	CHECK(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bnw_tne"]->id);
	CHECK(domain["tnw_bne"]->getNormalNbrInfo(Side<3>::top()).id == domain["tnw_tne"]->id);

	CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tnw_tse"]->id);
	CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Side<3>::south()).id == domain["tsw_tnw"]->id);
	CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["tnw_tnw"]->id);
	CHECK(domain["tnw_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tnw_bsw"]->id);

	CHECK(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::west()).id == domain["tnw_tsw"]->id);
	CHECK(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::east()).id == domain["tne_tsw"]->id);
	CHECK(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::south()).id == domain["tsw_tne"]->id);
	CHECK(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::north()).id == domain["tnw_tne"]->id);
	CHECK(domain["tnw_tse"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tnw_bse"]->id);

	CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tnw_tne"]->id);
	CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["tnw_tsw"]->id);
	CHECK(domain["tnw_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tnw_bnw"]->id);

	CHECK(domain["tnw_tne"]->getNormalNbrInfo(Side<3>::west()).id == domain["tnw_tnw"]->id);
	CHECK(domain["tnw_tne"]->getNormalNbrInfo(Side<3>::east()).id == domain["tne_tnw"]->id);
	CHECK(domain["tnw_tne"]->getNormalNbrInfo(Side<3>::south()).id == domain["tnw_tse"]->id);
	CHECK(domain["tnw_tne"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tnw_bne"]->id);

	CHECK(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::west()).id == domain["tnw_bse"]->id);
	CHECK(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tne_bse"]->id);
	CHECK(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::south()).id == domain["tse_bnw"]->id);
	CHECK(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["tne_bnw"]->id);
	CHECK(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bne_tsw"]->id);
	CHECK(domain["tne_bsw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tne_tsw"]->id);

	CHECK(domain["tne_bse"]->getNormalNbrInfo(Side<3>::west()).id == domain["tne_bsw"]->id);
	CHECK(domain["tne_bse"]->getNormalNbrInfo(Side<3>::south()).id == domain["tse_bne"]->id);
	CHECK(domain["tne_bse"]->getNormalNbrInfo(Side<3>::north()).id == domain["tne_bne"]->id);
	CHECK(domain["tne_bse"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bne_tse"]->id);
	CHECK(domain["tne_bse"]->getNormalNbrInfo(Side<3>::top()).id == domain["tne_tse"]->id);

	CHECK(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::west()).id == domain["tnw_bne"]->id);
	CHECK(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tne_bne"]->id);
	CHECK(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["tne_bsw"]->id);
	CHECK(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bne_tnw"]->id);
	CHECK(domain["tne_bnw"]->getNormalNbrInfo(Side<3>::top()).id == domain["tne_tnw"]->id);

	CHECK(domain["tne_bne"]->getNormalNbrInfo(Side<3>::west()).id == domain["tne_bnw"]->id);
	CHECK(domain["tne_bne"]->getNormalNbrInfo(Side<3>::south()).id == domain["tne_bse"]->id);
	CHECK(domain["tne_bne"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["bne_tne"]->id);
	CHECK(domain["tne_bne"]->getNormalNbrInfo(Side<3>::top()).id == domain["tne_tne"]->id);

	CHECK(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::west()).id == domain["tnw_tse"]->id);
	CHECK(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tne_tse"]->id);
	CHECK(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::south()).id == domain["tse_tnw"]->id);
	CHECK(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::north()).id == domain["tne_tnw"]->id);
	CHECK(domain["tne_tsw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tne_bsw"]->id);

	CHECK(domain["tne_tse"]->getNormalNbrInfo(Side<3>::west()).id == domain["tne_tsw"]->id);
	CHECK(domain["tne_tse"]->getNormalNbrInfo(Side<3>::south()).id == domain["tse_tne"]->id);
	CHECK(domain["tne_tse"]->getNormalNbrInfo(Side<3>::north()).id == domain["tne_tne"]->id);
	CHECK(domain["tne_tse"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tne_bse"]->id);

	CHECK(domain["tne_tnw"]->getNormalNbrInfo(Side<3>::west()).id == domain["tnw_tne"]->id);
	CHECK(domain["tne_tnw"]->getNormalNbrInfo(Side<3>::east()).id == domain["tne_tne"]->id);
	CHECK(domain["tne_tnw"]->getNormalNbrInfo(Side<3>::south()).id == domain["tne_tsw"]->id);
	CHECK(domain["tne_tnw"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tne_bnw"]->id);

	CHECK(domain["tne_tne"]->getNormalNbrInfo(Side<3>::west()).id == domain["tne_tnw"]->id);
	CHECK(domain["tne_tne"]->getNormalNbrInfo(Side<3>::south()).id == domain["tne_tse"]->id);
	CHECK(domain["tne_tne"]->getNormalNbrInfo(Side<3>::bottom()).id == domain["tne_bne"]->id);
}