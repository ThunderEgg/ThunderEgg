/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
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

#include <ThunderEgg/Corner.h>
namespace ThunderEgg
{
void to_json(nlohmann::json &j, const Corner<1> &o)
{
	j = nullptr;
}
void to_json(nlohmann::json &j, const Corner<2> &o)
{
	if (o == Corner<2>::sw()) {
		j = "SW";
	} else if (o == Corner<2>::se()) {
		j = "SE";
	} else if (o == Corner<2>::nw()) {
		j = "NW";
	} else if (o == Corner<2>::ne()) {
		j = "NE";
	} else {
		j = nullptr;
	}
}
void to_json(nlohmann::json &j, const Corner<3> &o)
{
	if (o == Corner<3>::bsw()) {
		j = "BSW";
	} else if (o == Corner<3>::bse()) {
		j = "BSE";
	} else if (o == Corner<3>::bnw()) {
		j = "BNW";
	} else if (o == Corner<3>::bne()) {
		j = "BNE";
	} else if (o == Corner<3>::tsw()) {
		j = "TSW";
	} else if (o == Corner<3>::tse()) {
		j = "TSE";
	} else if (o == Corner<3>::tnw()) {
		j = "TNW";
	} else if (o == Corner<3>::tne()) {
		j = "TNE";
	} else {
		j = nullptr;
	}
}
void from_json(const nlohmann::json &j, Corner<1> &o)
{
	o = Corner<1>::null();
}
void from_json(const nlohmann::json &j, Corner<2> &o)
{
	if (j == "SW") {
		o = Corner<2>::sw();
	} else if (j == "SE") {
		o = Corner<2>::se();
	} else if (j == "NW") {
		o = Corner<2>::nw();
	} else if (j == "NE") {
		o = Corner<2>::ne();
	} else {
		o = Corner<2>::null();
	}
}
void from_json(const nlohmann::json &j, Corner<3> &o)
{
	if (j == "BSW") {
		o = Corner<3>::bsw();
	} else if (j == "BSE") {
		o = Corner<3>::bse();
	} else if (j == "BNW") {
		o = Corner<3>::bnw();
	} else if (j == "BNE") {
		o = Corner<3>::bne();
	} else if (j == "TSW") {
		o = Corner<3>::tsw();
	} else if (j == "TSE") {
		o = Corner<3>::tse();
	} else if (j == "TNW") {
		o = Corner<3>::tnw();
	} else if (j == "TNE") {
		o = Corner<3>::tne();
	} else {
		o = Corner<3>::null();
	}
}
} // namespace ThunderEgg