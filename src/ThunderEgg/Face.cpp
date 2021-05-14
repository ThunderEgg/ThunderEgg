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

#include <ThunderEgg/Face.h>
namespace ThunderEgg
{
void to_json(nlohmann::json &j, const Side<1> &s)
{
	if (s == Side<1>::west()) {
		j = "WEST";
	} else if (s == Side<1>::east()) {
		j = "EAST";
	} else {
		j = nullptr;
	}
}
void to_json(nlohmann::json &j, const Side<2> &s)
{
	if (s == Side<2>::west()) {
		j = "WEST";
	} else if (s == Side<2>::east()) {
		j = "EAST";
	} else if (s == Side<2>::south()) {
		j = "SOUTH";
	} else if (s == Side<2>::north()) {
		j = "NORTH";
	} else {
		j = nullptr;
	}
}
void to_json(nlohmann::json &j, const Side<3> &s)
{
	if (s == Side<3>::west()) {
		j = "WEST";
	} else if (s == Side<3>::east()) {
		j = "EAST";
	} else if (s == Side<3>::south()) {
		j = "SOUTH";
	} else if (s == Side<3>::north()) {
		j = "NORTH";
	} else if (s == Side<3>::bottom()) {
		j = "BOTTOM";
	} else if (s == Side<3>::top()) {
		j = "TOP";
	} else {
		j = nullptr;
	}
}
void from_json(const nlohmann::json &j, Side<1> &s)
{
	if (j == "WEST") {
		s = Side<1>::west();
	} else if (j == "EAST") {
		s = Side<1>::east();
	} else {
		s = Side<1>::null();
	}
}
void from_json(const nlohmann::json &j, Side<2> &s)
{
	if (j == "WEST") {
		s = Side<2>::west();
	} else if (j == "EAST") {
		s = Side<2>::east();
	} else if (j == "SOUTH") {
		s = Side<2>::south();
	} else if (j == "NORTH") {
		s = Side<2>::north();
	} else {
		s = Side<2>::null();
	}
}
void from_json(const nlohmann::json &j, Side<3> &s)
{
	if (j == "WEST") {
		s = Side<3>::west();
	} else if (j == "EAST") {
		s = Side<3>::east();
	} else if (j == "SOUTH") {
		s = Side<3>::south();
	} else if (j == "NORTH") {
		s = Side<3>::north();
	} else if (j == "BOTTOM") {
		s = Side<3>::bottom();
	} else if (j == "TOP") {
		s = Side<3>::top();
	} else {
		s = Side<3>::null();
	}
}
std::ostream &operator<<(std::ostream &os, const Side<1> &s)
{
	if (s == Side<1>::east()) {
		os << "Side<1>::east()";
	} else if (s == Side<1>::west()) {
		os << "Side<1>::west()";
	} else if (s == Side<1>::null()) {
		os << "Side<1>::null()";
	} else {
		os << "Side<1> undefined value: " << s.getIndex();
	}
	return os;
}
std::ostream &operator<<(std::ostream &os, const Side<2> &s)
{
	if (s == Side<2>::east()) {
		os << "Side<2>::east()";
	} else if (s == Side<2>::west()) {
		os << "Side<2>::west()";
	} else if (s == Side<2>::south()) {
		os << "Side<2>::south()";
	} else if (s == Side<2>::north()) {
		os << "Side<2>::north()";
	} else if (s == Side<2>::null()) {
		os << "Side<2>::null()";
	} else {
		os << "Side<2> undefined value: " << s.getIndex();
	}
	return os;
}
std::ostream &operator<<(std::ostream &os, const Side<3> &s)
{
	if (s == Side<3>::east()) {
		os << "Side<3>::east()";
	} else if (s == Side<3>::west()) {
		os << "Side<3>::west()";
	} else if (s == Side<3>::south()) {
		os << "Side<3>::south()";
	} else if (s == Side<3>::north()) {
		os << "Side<3>::north()";
	} else if (s == Side<3>::bottom()) {
		os << "Side<3>::bottom()";
	} else if (s == Side<3>::top()) {
		os << "Side<3>::top()";
	} else if (s == Side<3>::null()) {
		os << "Side<3>::null()";
	} else {
		os << "Side<3> undefined value: " << s.getIndex();
	}
	return os;
}
void to_json(nlohmann::json &j, const Edge &o)
{
	if (o == Edge::bs()) {
		j = "BS";
	} else if (o == Edge::bn()) {
		j = "BN";
	} else if (o == Edge::ts()) {
		j = "TS";
	} else if (o == Edge::tn()) {
		j = "TN";
	} else if (o == Edge::bw()) {
		j = "BW";
	} else if (o == Edge::be()) {
		j = "BE";
	} else if (o == Edge::tw()) {
		j = "TW";
	} else if (o == Edge::te()) {
		j = "TE";
	} else if (o == Edge::sw()) {
		j = "SW";
	} else if (o == Edge::se()) {
		j = "SE";
	} else if (o == Edge::nw()) {
		j = "NW";
	} else if (o == Edge::ne()) {
		j = "NE";
	} else {
		j = nullptr;
	}
}
void from_json(const nlohmann::json &j, Edge &o)
{
	if (j == "BS") {
		o = Edge::bs();
	} else if (j == "BN") {
		o = Edge::bn();
	} else if (j == "TS") {
		o = Edge::ts();
	} else if (j == "TN") {
		o = Edge::tn();
	} else if (j == "BW") {
		o = Edge::bw();
	} else if (j == "BE") {
		o = Edge::be();
	} else if (j == "TW") {
		o = Edge::tw();
	} else if (j == "TE") {
		o = Edge::te();
	} else if (j == "SW") {
		o = Edge::sw();
	} else if (j == "SE") {
		o = Edge::se();
	} else if (j == "NW") {
		o = Edge::nw();
	} else if (j == "NE") {
		o = Edge::ne();
	} else {
		o = Edge::null();
	}
}
std::ostream &operator<<(std::ostream &os, const Edge &o)
{
	if (o == Edge::bs()) {
		os << "Edge::bs()";
	} else if (o == Edge::bn()) {
		os << "Edge::bn()";
	} else if (o == Edge::ts()) {
		os << "Edge::ts()";
	} else if (o == Edge::tn()) {
		os << "Edge::tn()";
	} else if (o == Edge::tw()) {
		os << "Edge::tw()";
	} else if (o == Edge::bw()) {
		os << "Edge::bw()";
	} else if (o == Edge::be()) {
		os << "Edge::be()";
	} else if (o == Edge::tw()) {
		os << "Edge::tw()";
	} else if (o == Edge::te()) {
		os << "Edge::te()";
	} else if (o == Edge::sw()) {
		os << "Edge::sw()";
	} else if (o == Edge::se()) {
		os << "Edge::se()";
	} else if (o == Edge::nw()) {
		os << "Edge::nw()";
	} else if (o == Edge::ne()) {
		os << "Edge::ne()";
	} else if (o == Edge::null()) {
		os << "Edge::null()";
	} else {
		os << "Edge invalid value: " << o.getIndex();
	}
	return os;
}
} // namespace ThunderEgg
