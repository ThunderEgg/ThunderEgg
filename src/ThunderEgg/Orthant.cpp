/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
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

#include <ThunderEgg/Orthant.h>
namespace ThunderEgg {
void
to_json(tpl::nlohmann::json& j, const Orthant<0>&)
{
  j = nullptr;
}
void
to_json(tpl::nlohmann::json& j, const Orthant<1>& o)
{
  if (o == Orthant<1>::lower()) {
    j = "LOWER";
  } else if (o == Orthant<1>::upper()) {
    j = "UPPER";
  } else {
    j = nullptr;
  }
}
void
to_json(tpl::nlohmann::json& j, const Orthant<2>& o)
{
  if (o == Orthant<2>::sw()) {
    j = "SW";
  } else if (o == Orthant<2>::se()) {
    j = "SE";
  } else if (o == Orthant<2>::nw()) {
    j = "NW";
  } else if (o == Orthant<2>::ne()) {
    j = "NE";
  } else {
    j = nullptr;
  }
}
void
to_json(tpl::nlohmann::json& j, const Orthant<3>& o)
{
  if (o == Orthant<3>::bsw()) {
    j = "BSW";
  } else if (o == Orthant<3>::bse()) {
    j = "BSE";
  } else if (o == Orthant<3>::bnw()) {
    j = "BNW";
  } else if (o == Orthant<3>::bne()) {
    j = "BNE";
  } else if (o == Orthant<3>::tsw()) {
    j = "TSW";
  } else if (o == Orthant<3>::tse()) {
    j = "TSE";
  } else if (o == Orthant<3>::tnw()) {
    j = "TNW";
  } else if (o == Orthant<3>::tne()) {
    j = "TNE";
  } else {
    j = nullptr;
  }
}
void
from_json(const tpl::nlohmann::json&, Orthant<0>& o)
{
  o = Orthant<0>::null();
}
void
from_json(const tpl::nlohmann::json& j, Orthant<1>& o)
{
  if (j == "LOWER") {
    o = Orthant<1>::lower();
  } else if (j == "UPPER") {
    o = Orthant<1>::upper();
  } else {
    o = Orthant<1>::null();
  }
}
void
from_json(const tpl::nlohmann::json& j, Orthant<2>& o)
{
  if (j == "SW") {
    o = Orthant<2>::sw();
  } else if (j == "SE") {
    o = Orthant<2>::se();
  } else if (j == "NW") {
    o = Orthant<2>::nw();
  } else if (j == "NE") {
    o = Orthant<2>::ne();
  } else {
    o = Orthant<2>::null();
  }
}
void
from_json(const tpl::nlohmann::json& j, Orthant<3>& o)
{
  if (j == "BSW") {
    o = Orthant<3>::bsw();
  } else if (j == "BSE") {
    o = Orthant<3>::bse();
  } else if (j == "BNW") {
    o = Orthant<3>::bnw();
  } else if (j == "BNE") {
    o = Orthant<3>::bne();
  } else if (j == "TSW") {
    o = Orthant<3>::tsw();
  } else if (j == "TSE") {
    o = Orthant<3>::tse();
  } else if (j == "TNW") {
    o = Orthant<3>::tnw();
  } else if (j == "TNE") {
    o = Orthant<3>::tne();
  } else {
    o = Orthant<3>::null();
  }
}
} // namespace ThunderEgg