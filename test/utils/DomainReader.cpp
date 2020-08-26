/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019-2020 ThunderEgg Developers. See AUTHORS.md file at the
 *  top-level directory.
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
#include "DomainReader.h"
template class DomainReader<2>;
template <> ThunderEgg::Orthant<2> DomainReader<2>::GetOrthant(nlohmann::json &orth_j)
{
	std::string            side_str = orth_j.get<std::string>();
	ThunderEgg::Orthant<2> orth;
	if (side_str == "SW") {
		orth = ThunderEgg::Orthant<2>::sw();
	} else if (side_str == "SE") {
		orth = ThunderEgg::Orthant<2>::se();
	} else if (side_str == "NW") {
		orth = ThunderEgg::Orthant<2>::nw();
	} else if (side_str == "NE") {
		orth = ThunderEgg::Orthant<2>::ne();
	} else {
		throw "parsing error";
	}
	return orth;
}
template <> ThunderEgg::Orthant<3> DomainReader<3>::GetOrthant(nlohmann::json &orth_j)
{
	std::string            side_str = orth_j.get<std::string>();
	ThunderEgg::Orthant<3> orth;
	if (side_str == "BSW") {
		orth = ThunderEgg::Orthant<3>::bsw();
	} else if (side_str == "BSE") {
		orth = ThunderEgg::Orthant<3>::bse();
	} else if (side_str == "BNW") {
		orth = ThunderEgg::Orthant<3>::bnw();
	} else if (side_str == "BNE") {
		orth = ThunderEgg::Orthant<3>::bne();
	} else if (side_str == "TSW") {
		orth = ThunderEgg::Orthant<3>::tsw();
	} else if (side_str == "TSE") {
		orth = ThunderEgg::Orthant<3>::tse();
	} else if (side_str == "TNW") {
		orth = ThunderEgg::Orthant<3>::tnw();
	} else if (side_str == "TNE") {
		orth = ThunderEgg::Orthant<3>::tne();
	} else {
		throw "parsing error";
	}
	return orth;
}