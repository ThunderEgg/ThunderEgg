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

#include <ThunderEgg/Side.h>
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
} // namespace ThunderEgg
