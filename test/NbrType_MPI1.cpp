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
#include <ThunderEgg/NbrType.h>

#include <doctest.h>
#include <sstream>

using namespace ThunderEgg;
using namespace ThunderEgg::tpl;

TEST_CASE("Test ostream for NbrType")
{
  std::stringstream ss;
  ss << NbrType::Coarse;
  CHECK_EQ(ss.str(), "NbrType::Coarse");
  ss.str("");
  ss << NbrType::Fine;
  CHECK_EQ(ss.str(), "NbrType::Fine");
  ss.str("");
  ss << NbrType::Normal;
  CHECK_EQ(ss.str(), "NbrType::Normal");
}