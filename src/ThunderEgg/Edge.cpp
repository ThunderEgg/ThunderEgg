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

#include <ThunderEgg/Edge.h>
namespace ThunderEgg
{
void to_json(nlohmann::json &j, const Edge<1> &o)
{
	j = nullptr;
}
void to_json(nlohmann::json &j, const Edge<2> &o)
{
	j = nullptr;
}
void to_json(nlohmann::json &j, const Edge<3> &o)
{
	if (o == Edge<3>::bs()) {
		j = "BS";
	} else if (o == Edge<3>::tn()) {
		j = "TN";
	} else if (o == Edge<3>::bn()) {
		j = "BN";
	} else if (o == Edge<3>::ts()) {
		j = "TS";
	} else if (o == Edge<3>::bw()) {
		j = "BW";
	} else if (o == Edge<3>::te()) {
		j = "TE";
	} else if (o == Edge<3>::be()) {
		j = "BE";
	} else if (o == Edge<3>::tw()) {
		j = "TW";
	} else if (o == Edge<3>::sw()) {
		j = "SW";
	} else if (o == Edge<3>::ne()) {
		j = "NE";
	} else if (o == Edge<3>::se()) {
		j = "SE";
	} else if (o == Edge<3>::nw()) {
		j = "NW";
	} else {
		j = nullptr;
	}
}
void from_json(const nlohmann::json &j, Edge<1> &o)
{
	o = Edge<1>::null();
}
void from_json(const nlohmann::json &j, Edge<2> &o)
{
	o = Edge<2>::null();
}
void from_json(const nlohmann::json &j, Edge<3> &o)
{
	if (j == "BS") {
		o = Edge<3>::bs();
	} else if (j == "TN") {
		o = Edge<3>::tn();
	} else if (j == "BN") {
		o = Edge<3>::bn();
	} else if (j == "TS") {
		o = Edge<3>::ts();
	} else if (j == "BW") {
		o = Edge<3>::bw();
	} else if (j == "TE") {
		o = Edge<3>::te();
	} else if (j == "BE") {
		o = Edge<3>::be();
	} else if (j == "TW") {
		o = Edge<3>::tw();
	} else if (j == "SW") {
		o = Edge<3>::sw();
	} else if (j == "NE") {
		o = Edge<3>::ne();
	} else if (j == "SE") {
		o = Edge<3>::se();
	} else if (j == "NW") {
		o = Edge<3>::nw();
	} else {
		o = Edge<3>::null();
	}
}
} // namespace ThunderEgg