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

#include <sstream>

#include <catch2/catch_test_macros.hpp>

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::tpl;

TEST_CASE("Orthant<0> unsigned char constructor works", "[Orthant]")
{
  Orthant<0> o(13);
  CHECK(o.getIndex() == 13);
}
TEST_CASE("Orthant<1> unsigned char constructor works", "[Orthant]")
{
  Orthant<1> o(13);
  CHECK(o.getIndex() == 13);
}
TEST_CASE("Orthant<2> unsigned char constructor works", "[Orthant]")
{
  Orthant<2> o(13);
  CHECK(o.getIndex() == 13);
}
TEST_CASE("Orthant<3> unsigned char constructor works", "[Orthant]")
{
  Orthant<3> o(13);
  CHECK(o.getIndex() == 13);
}
TEST_CASE("Orthant<0> Default constructor works", "[Orthant]")
{
  Orthant<0> o;
  CHECK(o == Orthant<0>::null());
}
TEST_CASE("Orthant<1> Default constructor works", "[Orthant]")
{
  Orthant<1> o;
  CHECK(o == Orthant<1>::null());
}
TEST_CASE("Orthant<2> Default constructor works", "[Orthant]")
{
  Orthant<2> o;
  CHECK(o == Orthant<2>::null());
}
TEST_CASE("Orthant<3> Default constructor works", "[Orthant]")
{
  Orthant<3> o;
  CHECK(o == Orthant<3>::null());
}
TEST_CASE("Orthant<0> named constructors give expected index values", "[Orthant]")
{
  CHECK(Orthant<0>::null().getIndex() == 1);
}
TEST_CASE("Orthant<1> named constructors give expected index values", "[Orthant]")
{
  CHECK(Orthant<1>::lower().getIndex() == 0);
  CHECK(Orthant<1>::upper().getIndex() == 1);
  CHECK(Orthant<1>::null().getIndex() == 2);
}
TEST_CASE("Orthant<2> named constructors give expected index values", "[Orthant]")
{
  CHECK(Orthant<2>::sw().getIndex() == 0);
  CHECK(Orthant<2>::se().getIndex() == 1);
  CHECK(Orthant<2>::nw().getIndex() == 2);
  CHECK(Orthant<2>::ne().getIndex() == 3);
  CHECK(Orthant<2>::null().getIndex() == 4);
}
TEST_CASE("Orthant<3> named constructors give expected index values", "[Orthant]")
{
  CHECK(Orthant<3>::bsw().getIndex() == 0);
  CHECK(Orthant<3>::bse().getIndex() == 1);
  CHECK(Orthant<3>::bnw().getIndex() == 2);
  CHECK(Orthant<3>::bne().getIndex() == 3);
  CHECK(Orthant<3>::tsw().getIndex() == 4);
  CHECK(Orthant<3>::tse().getIndex() == 5);
  CHECK(Orthant<3>::tnw().getIndex() == 6);
  CHECK(Orthant<3>::tne().getIndex() == 7);
  CHECK(Orthant<3>::null().getIndex() == 8);
}
TEST_CASE("Orthant<1> getNbrOnSide is as expected", "[Orthant]")
{
  CHECK(Orthant<1>::lower().getNbrOnSide(Side<1>::west()) == Orthant<1>::upper());
  CHECK(Orthant<1>::lower().getNbrOnSide(Side<1>::east()) == Orthant<1>::upper());

  CHECK(Orthant<1>::upper().getNbrOnSide(Side<1>::west()) == Orthant<1>::lower());
  CHECK(Orthant<1>::upper().getNbrOnSide(Side<1>::east()) == Orthant<1>::lower());
}
TEST_CASE("Orthant<2> getNbrOnSide is as expected", "[Orthant]")
{
  CHECK(Orthant<2>::sw().getNbrOnSide(Side<2>::west()) == Orthant<2>::se());
  CHECK(Orthant<2>::sw().getNbrOnSide(Side<2>::east()) == Orthant<2>::se());
  CHECK(Orthant<2>::sw().getNbrOnSide(Side<2>::south()) == Orthant<2>::nw());
  CHECK(Orthant<2>::sw().getNbrOnSide(Side<2>::north()) == Orthant<2>::nw());

  CHECK(Orthant<2>::se().getNbrOnSide(Side<2>::west()) == Orthant<2>::sw());
  CHECK(Orthant<2>::se().getNbrOnSide(Side<2>::east()) == Orthant<2>::sw());
  CHECK(Orthant<2>::se().getNbrOnSide(Side<2>::south()) == Orthant<2>::ne());
  CHECK(Orthant<2>::se().getNbrOnSide(Side<2>::north()) == Orthant<2>::ne());

  CHECK(Orthant<2>::nw().getNbrOnSide(Side<2>::west()) == Orthant<2>::ne());
  CHECK(Orthant<2>::nw().getNbrOnSide(Side<2>::east()) == Orthant<2>::ne());
  CHECK(Orthant<2>::nw().getNbrOnSide(Side<2>::south()) == Orthant<2>::sw());
  CHECK(Orthant<2>::nw().getNbrOnSide(Side<2>::north()) == Orthant<2>::sw());

  CHECK(Orthant<2>::ne().getNbrOnSide(Side<2>::west()) == Orthant<2>::nw());
  CHECK(Orthant<2>::ne().getNbrOnSide(Side<2>::east()) == Orthant<2>::nw());
  CHECK(Orthant<2>::ne().getNbrOnSide(Side<2>::south()) == Orthant<2>::se());
  CHECK(Orthant<2>::ne().getNbrOnSide(Side<2>::north()) == Orthant<2>::se());
}
TEST_CASE("Orthant<3> getNbrOnSide is as expected", "[Orthant]")
{
  CHECK(Orthant<3>::bsw().getNbrOnSide(Side<3>::west()) == Orthant<3>::bse());
  CHECK(Orthant<3>::bsw().getNbrOnSide(Side<3>::east()) == Orthant<3>::bse());
  CHECK(Orthant<3>::bsw().getNbrOnSide(Side<3>::south()) == Orthant<3>::bnw());
  CHECK(Orthant<3>::bsw().getNbrOnSide(Side<3>::north()) == Orthant<3>::bnw());
  CHECK(Orthant<3>::bsw().getNbrOnSide(Side<3>::bottom()) == Orthant<3>::tsw());
  CHECK(Orthant<3>::bsw().getNbrOnSide(Side<3>::top()) == Orthant<3>::tsw());

  CHECK(Orthant<3>::bse().getNbrOnSide(Side<3>::west()) == Orthant<3>::bsw());
  CHECK(Orthant<3>::bse().getNbrOnSide(Side<3>::east()) == Orthant<3>::bsw());
  CHECK(Orthant<3>::bse().getNbrOnSide(Side<3>::south()) == Orthant<3>::bne());
  CHECK(Orthant<3>::bse().getNbrOnSide(Side<3>::north()) == Orthant<3>::bne());
  CHECK(Orthant<3>::bse().getNbrOnSide(Side<3>::bottom()) == Orthant<3>::tse());
  CHECK(Orthant<3>::bse().getNbrOnSide(Side<3>::top()) == Orthant<3>::tse());

  CHECK(Orthant<3>::bnw().getNbrOnSide(Side<3>::west()) == Orthant<3>::bne());
  CHECK(Orthant<3>::bnw().getNbrOnSide(Side<3>::east()) == Orthant<3>::bne());
  CHECK(Orthant<3>::bnw().getNbrOnSide(Side<3>::south()) == Orthant<3>::bsw());
  CHECK(Orthant<3>::bnw().getNbrOnSide(Side<3>::north()) == Orthant<3>::bsw());
  CHECK(Orthant<3>::bnw().getNbrOnSide(Side<3>::bottom()) == Orthant<3>::tnw());
  CHECK(Orthant<3>::bnw().getNbrOnSide(Side<3>::top()) == Orthant<3>::tnw());

  CHECK(Orthant<3>::bne().getNbrOnSide(Side<3>::west()) == Orthant<3>::bnw());
  CHECK(Orthant<3>::bne().getNbrOnSide(Side<3>::east()) == Orthant<3>::bnw());
  CHECK(Orthant<3>::bne().getNbrOnSide(Side<3>::south()) == Orthant<3>::bse());
  CHECK(Orthant<3>::bne().getNbrOnSide(Side<3>::north()) == Orthant<3>::bse());
  CHECK(Orthant<3>::bne().getNbrOnSide(Side<3>::bottom()) == Orthant<3>::tne());
  CHECK(Orthant<3>::bne().getNbrOnSide(Side<3>::top()) == Orthant<3>::tne());

  CHECK(Orthant<3>::tsw().getNbrOnSide(Side<3>::west()) == Orthant<3>::tse());
  CHECK(Orthant<3>::tsw().getNbrOnSide(Side<3>::east()) == Orthant<3>::tse());
  CHECK(Orthant<3>::tsw().getNbrOnSide(Side<3>::south()) == Orthant<3>::tnw());
  CHECK(Orthant<3>::tsw().getNbrOnSide(Side<3>::north()) == Orthant<3>::tnw());
  CHECK(Orthant<3>::tsw().getNbrOnSide(Side<3>::bottom()) == Orthant<3>::bsw());
  CHECK(Orthant<3>::tsw().getNbrOnSide(Side<3>::top()) == Orthant<3>::bsw());

  CHECK(Orthant<3>::tse().getNbrOnSide(Side<3>::east()) == Orthant<3>::tsw());
  CHECK(Orthant<3>::tse().getNbrOnSide(Side<3>::west()) == Orthant<3>::tsw());
  CHECK(Orthant<3>::tse().getNbrOnSide(Side<3>::south()) == Orthant<3>::tne());
  CHECK(Orthant<3>::tse().getNbrOnSide(Side<3>::north()) == Orthant<3>::tne());
  CHECK(Orthant<3>::tse().getNbrOnSide(Side<3>::bottom()) == Orthant<3>::bse());
  CHECK(Orthant<3>::tse().getNbrOnSide(Side<3>::top()) == Orthant<3>::bse());

  CHECK(Orthant<3>::tnw().getNbrOnSide(Side<3>::west()) == Orthant<3>::tne());
  CHECK(Orthant<3>::tnw().getNbrOnSide(Side<3>::east()) == Orthant<3>::tne());
  CHECK(Orthant<3>::tnw().getNbrOnSide(Side<3>::south()) == Orthant<3>::tsw());
  CHECK(Orthant<3>::tnw().getNbrOnSide(Side<3>::north()) == Orthant<3>::tsw());
  CHECK(Orthant<3>::tnw().getNbrOnSide(Side<3>::bottom()) == Orthant<3>::bnw());
  CHECK(Orthant<3>::tnw().getNbrOnSide(Side<3>::top()) == Orthant<3>::bnw());

  CHECK(Orthant<3>::tne().getNbrOnSide(Side<3>::east()) == Orthant<3>::tnw());
  CHECK(Orthant<3>::tne().getNbrOnSide(Side<3>::west()) == Orthant<3>::tnw());
  CHECK(Orthant<3>::tne().getNbrOnSide(Side<3>::south()) == Orthant<3>::tse());
  CHECK(Orthant<3>::tne().getNbrOnSide(Side<3>::north()) == Orthant<3>::tse());
  CHECK(Orthant<3>::tne().getNbrOnSide(Side<3>::bottom()) == Orthant<3>::bne());
  CHECK(Orthant<3>::tne().getNbrOnSide(Side<3>::top()) == Orthant<3>::bne());
}
TEST_CASE("Orthant<1> getInteriorSides is as expected", "[Orthant]")
{
  {
    auto array = Orthant<1>::lower().getInteriorSides();
    CHECK(array[0] == Side<1>::east());
  }
  {
    auto array = Orthant<1>::upper().getInteriorSides();
    CHECK(array[0] == Side<1>::west());
  }
}
TEST_CASE("Orthant<2> getInteriorSides is as expected", "[Orthant]")
{
  {
    auto array = Orthant<2>::sw().getInteriorSides();
    CHECK(array[0] == Side<2>::east());
    CHECK(array[1] == Side<2>::north());
  }
  {
    auto array = Orthant<2>::se().getInteriorSides();
    CHECK(array[0] == Side<2>::west());
    CHECK(array[1] == Side<2>::north());
  }
  {
    auto array = Orthant<2>::nw().getInteriorSides();
    CHECK(array[0] == Side<2>::east());
    CHECK(array[1] == Side<2>::south());
  }
  {
    auto array = Orthant<2>::ne().getInteriorSides();
    CHECK(array[0] == Side<2>::west());
    CHECK(array[1] == Side<2>::south());
  }
}
TEST_CASE("Orthant<3> getInteriorSides is as expected", "[Orthant]")
{
  {
    auto array = Orthant<3>::bsw().getInteriorSides();
    CHECK(array[0] == Side<3>::east());
    CHECK(array[1] == Side<3>::north());
    CHECK(array[2] == Side<3>::top());
  }
  {
    auto array = Orthant<3>::bse().getInteriorSides();
    CHECK(array[0] == Side<3>::west());
    CHECK(array[1] == Side<3>::north());
    CHECK(array[2] == Side<3>::top());
  }
  {
    auto array = Orthant<3>::bnw().getInteriorSides();
    CHECK(array[0] == Side<3>::east());
    CHECK(array[1] == Side<3>::south());
    CHECK(array[2] == Side<3>::top());
  }
  {
    auto array = Orthant<3>::bne().getInteriorSides();
    CHECK(array[0] == Side<3>::west());
    CHECK(array[1] == Side<3>::south());
    CHECK(array[2] == Side<3>::top());
  }
  {
    auto array = Orthant<3>::tsw().getInteriorSides();
    CHECK(array[0] == Side<3>::east());
    CHECK(array[1] == Side<3>::north());
    CHECK(array[2] == Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::tse().getInteriorSides();
    CHECK(array[0] == Side<3>::west());
    CHECK(array[1] == Side<3>::north());
    CHECK(array[2] == Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::tnw().getInteriorSides();
    CHECK(array[0] == Side<3>::east());
    CHECK(array[1] == Side<3>::south());
    CHECK(array[2] == Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::tne().getInteriorSides();
    CHECK(array[0] == Side<3>::west());
    CHECK(array[1] == Side<3>::south());
    CHECK(array[2] == Side<3>::bottom());
  }
}
TEST_CASE("Orthant<1> getExteriorSides is as expected", "[Orthant]")
{
  {
    auto array = Orthant<1>::lower().getExteriorSides();
    CHECK(array[0] == Side<1>::west());
  }
  {
    auto array = Orthant<1>::upper().getExteriorSides();
    CHECK(array[0] == Side<1>::east());
  }
}
TEST_CASE("Orthant<2> getExteriorSides is as expected", "[Orthant]")
{
  {
    auto array = Orthant<2>::sw().getExteriorSides();
    CHECK(array[0] == Side<2>::west());
    CHECK(array[1] == Side<2>::south());
  }
  {
    auto array = Orthant<2>::se().getExteriorSides();
    CHECK(array[0] == Side<2>::east());
    CHECK(array[1] == Side<2>::south());
  }
  {
    auto array = Orthant<2>::nw().getExteriorSides();
    CHECK(array[0] == Side<2>::west());
    CHECK(array[1] == Side<2>::north());
  }
  {
    auto array = Orthant<2>::ne().getExteriorSides();
    CHECK(array[0] == Side<2>::east());
    CHECK(array[1] == Side<2>::north());
  }
}
TEST_CASE("Orthant<3> getExteriorSides is as expected", "[Orthant]")
{
  {
    auto array = Orthant<3>::bsw().getExteriorSides();
    CHECK(array[0] == Side<3>::west());
    CHECK(array[1] == Side<3>::south());
    CHECK(array[2] == Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::bse().getExteriorSides();
    CHECK(array[0] == Side<3>::east());
    CHECK(array[1] == Side<3>::south());
    CHECK(array[2] == Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::bnw().getExteriorSides();
    CHECK(array[0] == Side<3>::west());
    CHECK(array[1] == Side<3>::north());
    CHECK(array[2] == Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::bne().getExteriorSides();
    CHECK(array[0] == Side<3>::east());
    CHECK(array[1] == Side<3>::north());
    CHECK(array[2] == Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::tsw().getExteriorSides();
    CHECK(array[0] == Side<3>::west());
    CHECK(array[1] == Side<3>::south());
    CHECK(array[2] == Side<3>::top());
  }
  {
    auto array = Orthant<3>::tse().getExteriorSides();
    CHECK(array[0] == Side<3>::east());
    CHECK(array[1] == Side<3>::south());
    CHECK(array[2] == Side<3>::top());
  }
  {
    auto array = Orthant<3>::tnw().getExteriorSides();
    CHECK(array[0] == Side<3>::west());
    CHECK(array[1] == Side<3>::north());
    CHECK(array[2] == Side<3>::top());
  }
  {
    auto array = Orthant<3>::tne().getExteriorSides();
    CHECK(array[0] == Side<3>::east());
    CHECK(array[1] == Side<3>::north());
    CHECK(array[2] == Side<3>::top());
  }
}
TEST_CASE("Orthant<1> isOnSide is as expected", "[Orthant]")
{
  CHECK(Orthant<1>::lower().isOnSide(Side<1>::west()));
  CHECK_FALSE(Orthant<1>::lower().isOnSide(Side<1>::east()));

  CHECK_FALSE(Orthant<1>::upper().isOnSide(Side<1>::west()));
  CHECK(Orthant<1>::upper().isOnSide(Side<1>::east()));
}
TEST_CASE("Orthant<2> isOnSide is as expected", "[Orthant]")
{
  CHECK(Orthant<2>::sw().isOnSide(Side<2>::west()));
  CHECK_FALSE(Orthant<2>::sw().isOnSide(Side<2>::east()));
  CHECK(Orthant<2>::sw().isOnSide(Side<2>::south()));
  CHECK_FALSE(Orthant<2>::sw().isOnSide(Side<2>::north()));

  CHECK_FALSE(Orthant<2>::se().isOnSide(Side<2>::west()));
  CHECK(Orthant<2>::se().isOnSide(Side<2>::east()));
  CHECK(Orthant<2>::se().isOnSide(Side<2>::south()));
  CHECK_FALSE(Orthant<2>::se().isOnSide(Side<2>::north()));

  CHECK(Orthant<2>::nw().isOnSide(Side<2>::west()));
  CHECK_FALSE(Orthant<2>::nw().isOnSide(Side<2>::east()));
  CHECK_FALSE(Orthant<2>::nw().isOnSide(Side<2>::south()));
  CHECK(Orthant<2>::nw().isOnSide(Side<2>::north()));

  CHECK_FALSE(Orthant<2>::ne().isOnSide(Side<2>::west()));
  CHECK(Orthant<2>::ne().isOnSide(Side<2>::east()));
  CHECK_FALSE(Orthant<2>::ne().isOnSide(Side<2>::south()));
  CHECK(Orthant<2>::ne().isOnSide(Side<2>::north()));
}
TEST_CASE("Orthant<3> isOnSide is as expected", "[Orthant]")
{
  CHECK(Orthant<3>::bsw().isOnSide(Side<3>::west()));
  CHECK_FALSE(Orthant<3>::bsw().isOnSide(Side<3>::east()));
  CHECK(Orthant<3>::bsw().isOnSide(Side<3>::south()));
  CHECK_FALSE(Orthant<3>::bsw().isOnSide(Side<3>::north()));
  CHECK(Orthant<3>::bsw().isOnSide(Side<3>::bottom()));
  CHECK_FALSE(Orthant<3>::bsw().isOnSide(Side<3>::top()));

  CHECK_FALSE(Orthant<3>::bse().isOnSide(Side<3>::west()));
  CHECK(Orthant<3>::bse().isOnSide(Side<3>::east()));
  CHECK(Orthant<3>::bse().isOnSide(Side<3>::south()));
  CHECK_FALSE(Orthant<3>::bse().isOnSide(Side<3>::north()));
  CHECK(Orthant<3>::bse().isOnSide(Side<3>::bottom()));
  CHECK_FALSE(Orthant<3>::bse().isOnSide(Side<3>::top()));

  CHECK(Orthant<3>::bnw().isOnSide(Side<3>::west()));
  CHECK_FALSE(Orthant<3>::bnw().isOnSide(Side<3>::east()));
  CHECK_FALSE(Orthant<3>::bnw().isOnSide(Side<3>::south()));
  CHECK(Orthant<3>::bnw().isOnSide(Side<3>::north()));
  CHECK(Orthant<3>::bnw().isOnSide(Side<3>::bottom()));
  CHECK_FALSE(Orthant<3>::bnw().isOnSide(Side<3>::top()));

  CHECK_FALSE(Orthant<3>::bne().isOnSide(Side<3>::west()));
  CHECK(Orthant<3>::bne().isOnSide(Side<3>::east()));
  CHECK_FALSE(Orthant<3>::bne().isOnSide(Side<3>::south()));
  CHECK(Orthant<3>::bne().isOnSide(Side<3>::north()));
  CHECK(Orthant<3>::bne().isOnSide(Side<3>::bottom()));
  CHECK_FALSE(Orthant<3>::bne().isOnSide(Side<3>::top()));

  CHECK(Orthant<3>::tsw().isOnSide(Side<3>::west()));
  CHECK_FALSE(Orthant<3>::tsw().isOnSide(Side<3>::east()));
  CHECK(Orthant<3>::tsw().isOnSide(Side<3>::south()));
  CHECK_FALSE(Orthant<3>::tsw().isOnSide(Side<3>::north()));
  CHECK_FALSE(Orthant<3>::tsw().isOnSide(Side<3>::bottom()));
  CHECK(Orthant<3>::tsw().isOnSide(Side<3>::top()));

  CHECK_FALSE(Orthant<3>::tse().isOnSide(Side<3>::west()));
  CHECK(Orthant<3>::tse().isOnSide(Side<3>::east()));
  CHECK(Orthant<3>::tse().isOnSide(Side<3>::south()));
  CHECK_FALSE(Orthant<3>::tse().isOnSide(Side<3>::north()));
  CHECK_FALSE(Orthant<3>::tse().isOnSide(Side<3>::bottom()));
  CHECK(Orthant<3>::tse().isOnSide(Side<3>::top()));

  CHECK(Orthant<3>::tnw().isOnSide(Side<3>::west()));
  CHECK_FALSE(Orthant<3>::tnw().isOnSide(Side<3>::east()));
  CHECK_FALSE(Orthant<3>::tnw().isOnSide(Side<3>::south()));
  CHECK(Orthant<3>::tnw().isOnSide(Side<3>::north()));
  CHECK_FALSE(Orthant<3>::tnw().isOnSide(Side<3>::bottom()));
  CHECK(Orthant<3>::tnw().isOnSide(Side<3>::top()));

  CHECK_FALSE(Orthant<3>::tne().isOnSide(Side<3>::west()));
  CHECK(Orthant<3>::tne().isOnSide(Side<3>::east()));
  CHECK_FALSE(Orthant<3>::tne().isOnSide(Side<3>::south()));
  CHECK(Orthant<3>::tne().isOnSide(Side<3>::north()));
  CHECK_FALSE(Orthant<3>::tne().isOnSide(Side<3>::bottom()));
  CHECK(Orthant<3>::tne().isOnSide(Side<3>::top()));
}
TEST_CASE("Orthant<1> getValuesOnSide is as expected", "[Orthant]")
{
  SECTION("Side<1>::west()")
  {
    std::array<Orthant<1>, 1> values = Orthant<1>::getValuesOnSide(Side<1>::west());
    CHECK(values[0] == Orthant<1>::lower());
  }
  SECTION("Side<1>::east()")
  {
    std::array<Orthant<1>, 1> values = Orthant<1>::getValuesOnSide(Side<1>::east());
    CHECK(values[0] == Orthant<1>::upper());
  }
}
TEST_CASE("Orthant<2> getValuesOnSide is as expected", "[Orthant]")
{
  {
    std::array<Orthant<2>, 2> values = Orthant<2>::getValuesOnSide(Side<2>::west());
    CHECK(values[0] == Orthant<2>::sw());
    CHECK(values[1] == Orthant<2>::nw());
  }

  {
    std::array<Orthant<2>, 2> values = Orthant<2>::getValuesOnSide(Side<2>::east());
    CHECK(values[0] == Orthant<2>::se());
    CHECK(values[1] == Orthant<2>::ne());
  }

  {
    std::array<Orthant<2>, 2> values = Orthant<2>::getValuesOnSide(Side<2>::south());
    CHECK(values[0] == Orthant<2>::sw());
    CHECK(values[1] == Orthant<2>::se());
  }

  {
    std::array<Orthant<2>, 2> values = Orthant<2>::getValuesOnSide(Side<2>::north());
    CHECK(values[0] == Orthant<2>::nw());
    CHECK(values[1] == Orthant<2>::ne());
  }
}
TEST_CASE("Orthant<3> getValuesOnSide is as expected", "[Orthant]")
{
  {
    std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::west());
    CHECK(values[0] == Orthant<3>::bsw());
    CHECK(values[1] == Orthant<3>::bnw());
    CHECK(values[2] == Orthant<3>::tsw());
    CHECK(values[3] == Orthant<3>::tnw());
  }

  {
    std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::east());
    CHECK(values[0] == Orthant<3>::bse());
    CHECK(values[1] == Orthant<3>::bne());
    CHECK(values[2] == Orthant<3>::tse());
    CHECK(values[3] == Orthant<3>::tne());
  }

  {
    std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::south());
    CHECK(values[0] == Orthant<3>::bsw());
    CHECK(values[1] == Orthant<3>::bse());
    CHECK(values[2] == Orthant<3>::tsw());
    CHECK(values[3] == Orthant<3>::tse());
  }

  {
    std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::north());
    CHECK(values[0] == Orthant<3>::bnw());
    CHECK(values[1] == Orthant<3>::bne());
    CHECK(values[2] == Orthant<3>::tnw());
    CHECK(values[3] == Orthant<3>::tne());
  }

  {
    std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::bottom());
    CHECK(values[0] == Orthant<3>::bsw());
    CHECK(values[1] == Orthant<3>::bse());
    CHECK(values[2] == Orthant<3>::bnw());
    CHECK(values[3] == Orthant<3>::bne());
  }

  {
    std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::top());
    CHECK(values[0] == Orthant<3>::tsw());
    CHECK(values[1] == Orthant<3>::tse());
    CHECK(values[2] == Orthant<3>::tnw());
    CHECK(values[3] == Orthant<3>::tne());
  }
}
TEST_CASE("Orthant<1> collapseOnAxis is as expected", "[Orthant]")
{
  CHECK(Orthant<1>::lower().collapseOnAxis(0) == Orthant<0>(0));
  CHECK(Orthant<1>::upper().collapseOnAxis(0) == Orthant<0>(0));
}
TEST_CASE("Orthant<2> collapseOnAxis is as expected", "[Orthant]")
{
  CHECK(Orthant<2>::sw().collapseOnAxis(0) == Orthant<1>::lower());
  CHECK(Orthant<2>::sw().collapseOnAxis(1) == Orthant<1>::lower());

  CHECK(Orthant<2>::se().collapseOnAxis(0) == Orthant<1>::lower());
  CHECK(Orthant<2>::se().collapseOnAxis(1) == Orthant<1>::upper());

  CHECK(Orthant<2>::nw().collapseOnAxis(0) == Orthant<1>::upper());
  CHECK(Orthant<2>::nw().collapseOnAxis(1) == Orthant<1>::lower());

  CHECK(Orthant<2>::ne().collapseOnAxis(0) == Orthant<1>::upper());
  CHECK(Orthant<2>::ne().collapseOnAxis(1) == Orthant<1>::upper());
}
TEST_CASE("Orthant<3> collapseOnAxis is as expected", "[Orthant]")
{
  CHECK(Orthant<3>::bsw().collapseOnAxis(0) == Orthant<2>::sw());
  CHECK(Orthant<3>::bsw().collapseOnAxis(1) == Orthant<2>::sw());
  CHECK(Orthant<3>::bsw().collapseOnAxis(2) == Orthant<2>::sw());

  CHECK(Orthant<3>::bse().collapseOnAxis(0) == Orthant<2>::sw());
  CHECK(Orthant<3>::bse().collapseOnAxis(1) == Orthant<2>::se());
  CHECK(Orthant<3>::bse().collapseOnAxis(2) == Orthant<2>::se());

  CHECK(Orthant<3>::bnw().collapseOnAxis(0) == Orthant<2>::se());
  CHECK(Orthant<3>::bnw().collapseOnAxis(1) == Orthant<2>::sw());
  CHECK(Orthant<3>::bnw().collapseOnAxis(2) == Orthant<2>::nw());

  CHECK(Orthant<3>::bne().collapseOnAxis(0) == Orthant<2>::se());
  CHECK(Orthant<3>::bne().collapseOnAxis(1) == Orthant<2>::se());
  CHECK(Orthant<3>::bne().collapseOnAxis(2) == Orthant<2>::ne());

  CHECK(Orthant<3>::tsw().collapseOnAxis(0) == Orthant<2>::nw());
  CHECK(Orthant<3>::tsw().collapseOnAxis(1) == Orthant<2>::nw());
  CHECK(Orthant<3>::tsw().collapseOnAxis(2) == Orthant<2>::sw());

  CHECK(Orthant<3>::tse().collapseOnAxis(0) == Orthant<2>::nw());
  CHECK(Orthant<3>::tse().collapseOnAxis(1) == Orthant<2>::ne());
  CHECK(Orthant<3>::tse().collapseOnAxis(2) == Orthant<2>::se());

  CHECK(Orthant<3>::tnw().collapseOnAxis(0) == Orthant<2>::ne());
  CHECK(Orthant<3>::tnw().collapseOnAxis(1) == Orthant<2>::nw());
  CHECK(Orthant<3>::tnw().collapseOnAxis(2) == Orthant<2>::nw());

  CHECK(Orthant<3>::tne().collapseOnAxis(0) == Orthant<2>::ne());
  CHECK(Orthant<3>::tne().collapseOnAxis(1) == Orthant<2>::ne());
  CHECK(Orthant<3>::tne().collapseOnAxis(2) == Orthant<2>::ne());
}
TEST_CASE("Orthant<0> ==", "[Orthant]")
{
  CHECK(Orthant<0>::null() == Orthant<0>::null());
}
TEST_CASE("Orthant<1> ==", "[Orthant]")
{
  CHECK(Orthant<1>::lower() == Orthant<1>::lower());
  CHECK_FALSE(Orthant<1>::lower() == Orthant<1>::upper());
  CHECK_FALSE(Orthant<1>::lower() == Orthant<1>::null());

  CHECK_FALSE(Orthant<1>::upper() == Orthant<1>::lower());
  CHECK(Orthant<1>::upper() == Orthant<1>::upper());
  CHECK_FALSE(Orthant<1>::upper() == Orthant<1>::null());

  CHECK_FALSE(Orthant<1>::null() == Orthant<1>::lower());
  CHECK_FALSE(Orthant<1>::null() == Orthant<1>::upper());
  CHECK(Orthant<1>::null() == Orthant<1>::null());
}
TEST_CASE("Orthant<2> ==", "[Orthant]")
{
  CHECK(Orthant<2>::sw() == Orthant<2>::sw());
  CHECK_FALSE(Orthant<2>::sw() == Orthant<2>::se());
  CHECK_FALSE(Orthant<2>::sw() == Orthant<2>::nw());
  CHECK_FALSE(Orthant<2>::sw() == Orthant<2>::ne());
  CHECK_FALSE(Orthant<2>::sw() == Orthant<2>::null());

  CHECK_FALSE(Orthant<2>::se() == Orthant<2>::sw());
  CHECK(Orthant<2>::se() == Orthant<2>::se());
  CHECK_FALSE(Orthant<2>::se() == Orthant<2>::nw());
  CHECK_FALSE(Orthant<2>::se() == Orthant<2>::ne());
  CHECK_FALSE(Orthant<2>::se() == Orthant<2>::null());

  CHECK_FALSE(Orthant<2>::nw() == Orthant<2>::sw());
  CHECK_FALSE(Orthant<2>::nw() == Orthant<2>::se());
  CHECK(Orthant<2>::nw() == Orthant<2>::nw());
  CHECK_FALSE(Orthant<2>::nw() == Orthant<2>::ne());
  CHECK_FALSE(Orthant<2>::nw() == Orthant<2>::null());

  CHECK_FALSE(Orthant<2>::ne() == Orthant<2>::sw());
  CHECK_FALSE(Orthant<2>::ne() == Orthant<2>::se());
  CHECK_FALSE(Orthant<2>::ne() == Orthant<2>::nw());
  CHECK(Orthant<2>::ne() == Orthant<2>::ne());
  CHECK_FALSE(Orthant<2>::ne() == Orthant<2>::null());

  CHECK_FALSE(Orthant<2>::null() == Orthant<2>::sw());
  CHECK_FALSE(Orthant<2>::null() == Orthant<2>::se());
  CHECK_FALSE(Orthant<2>::null() == Orthant<2>::nw());
  CHECK_FALSE(Orthant<2>::null() == Orthant<2>::ne());
  CHECK(Orthant<2>::null() == Orthant<2>::null());
}
TEST_CASE("Orthant<3> ==", "[Orthant]")
{
  CHECK(Orthant<3>::bsw() == Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::bsw() == Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::bsw() == Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::bsw() == Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::bsw() == Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::bsw() == Orthant<3>::tse());
  CHECK_FALSE(Orthant<3>::bsw() == Orthant<3>::tnw());
  CHECK_FALSE(Orthant<3>::bsw() == Orthant<3>::tne());
  CHECK_FALSE(Orthant<3>::bsw() == Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::bse() == Orthant<3>::bsw());
  CHECK(Orthant<3>::bse() == Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::bse() == Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::bse() == Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::bse() == Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::bse() == Orthant<3>::tse());
  CHECK_FALSE(Orthant<3>::bse() == Orthant<3>::tnw());
  CHECK_FALSE(Orthant<3>::bse() == Orthant<3>::tne());
  CHECK_FALSE(Orthant<3>::bse() == Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::bnw() == Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::bnw() == Orthant<3>::bse());
  CHECK(Orthant<3>::bnw() == Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::bnw() == Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::bnw() == Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::bnw() == Orthant<3>::tse());
  CHECK_FALSE(Orthant<3>::bnw() == Orthant<3>::tnw());
  CHECK_FALSE(Orthant<3>::bnw() == Orthant<3>::tne());
  CHECK_FALSE(Orthant<3>::bnw() == Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::bne() == Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::bne() == Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::bne() == Orthant<3>::bnw());
  CHECK(Orthant<3>::bne() == Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::bne() == Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::bne() == Orthant<3>::tse());
  CHECK_FALSE(Orthant<3>::bne() == Orthant<3>::tnw());
  CHECK_FALSE(Orthant<3>::bne() == Orthant<3>::tne());
  CHECK_FALSE(Orthant<3>::bne() == Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::tsw() == Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::tsw() == Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::tsw() == Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::tsw() == Orthant<3>::bne());
  CHECK(Orthant<3>::tsw() == Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::tsw() == Orthant<3>::tse());
  CHECK_FALSE(Orthant<3>::tsw() == Orthant<3>::tnw());
  CHECK_FALSE(Orthant<3>::tsw() == Orthant<3>::tne());
  CHECK_FALSE(Orthant<3>::tsw() == Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::tse() == Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::tse() == Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::tse() == Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::tse() == Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::tse() == Orthant<3>::tsw());
  CHECK(Orthant<3>::tse() == Orthant<3>::tse());
  CHECK_FALSE(Orthant<3>::tse() == Orthant<3>::tnw());
  CHECK_FALSE(Orthant<3>::tse() == Orthant<3>::tne());
  CHECK_FALSE(Orthant<3>::tse() == Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::tnw() == Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::tnw() == Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::tnw() == Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::tnw() == Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::tnw() == Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::tnw() == Orthant<3>::tse());
  CHECK(Orthant<3>::tnw() == Orthant<3>::tnw());
  CHECK_FALSE(Orthant<3>::tnw() == Orthant<3>::tne());
  CHECK_FALSE(Orthant<3>::tnw() == Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::tne() == Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::tne() == Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::tne() == Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::tne() == Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::tne() == Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::tne() == Orthant<3>::tse());
  CHECK_FALSE(Orthant<3>::tne() == Orthant<3>::tnw());
  CHECK(Orthant<3>::tne() == Orthant<3>::tne());
  CHECK_FALSE(Orthant<3>::tne() == Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::null() == Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::null() == Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::null() == Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::null() == Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::null() == Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::null() == Orthant<3>::tse());
  CHECK_FALSE(Orthant<3>::null() == Orthant<3>::tnw());
  CHECK_FALSE(Orthant<3>::null() == Orthant<3>::tne());
  CHECK(Orthant<3>::null() == Orthant<3>::null());
}
TEST_CASE("Orthant<1> !=", "[Orthant]")
{
  CHECK_FALSE(Orthant<1>::lower() != Orthant<1>::lower());
  CHECK(Orthant<1>::lower() != Orthant<1>::upper());
  CHECK(Orthant<1>::lower() != Orthant<1>::null());

  CHECK(Orthant<1>::upper() != Orthant<1>::lower());
  CHECK_FALSE(Orthant<1>::upper() != Orthant<1>::upper());
  CHECK(Orthant<1>::upper() != Orthant<1>::null());

  CHECK(Orthant<1>::null() != Orthant<1>::lower());
  CHECK(Orthant<1>::null() != Orthant<1>::upper());
  CHECK_FALSE(Orthant<1>::null() != Orthant<1>::null());
}
TEST_CASE("Orthant<2> !=", "[Orthant]")
{
  CHECK_FALSE(Orthant<2>::sw() != Orthant<2>::sw());
  CHECK(Orthant<2>::sw() != Orthant<2>::se());
  CHECK(Orthant<2>::sw() != Orthant<2>::nw());
  CHECK(Orthant<2>::sw() != Orthant<2>::ne());
  CHECK(Orthant<2>::sw() != Orthant<2>::null());

  CHECK(Orthant<2>::se() != Orthant<2>::sw());
  CHECK_FALSE(Orthant<2>::se() != Orthant<2>::se());
  CHECK(Orthant<2>::se() != Orthant<2>::nw());
  CHECK(Orthant<2>::se() != Orthant<2>::ne());
  CHECK(Orthant<2>::se() != Orthant<2>::null());

  CHECK(Orthant<2>::nw() != Orthant<2>::sw());
  CHECK(Orthant<2>::nw() != Orthant<2>::se());
  CHECK_FALSE(Orthant<2>::nw() != Orthant<2>::nw());
  CHECK(Orthant<2>::nw() != Orthant<2>::ne());
  CHECK(Orthant<2>::nw() != Orthant<2>::null());

  CHECK(Orthant<2>::ne() != Orthant<2>::sw());
  CHECK(Orthant<2>::ne() != Orthant<2>::se());
  CHECK(Orthant<2>::ne() != Orthant<2>::nw());
  CHECK_FALSE(Orthant<2>::ne() != Orthant<2>::ne());
  CHECK(Orthant<2>::ne() != Orthant<2>::null());

  CHECK(Orthant<2>::null() != Orthant<2>::sw());
  CHECK(Orthant<2>::null() != Orthant<2>::se());
  CHECK(Orthant<2>::null() != Orthant<2>::nw());
  CHECK(Orthant<2>::null() != Orthant<2>::ne());
  CHECK_FALSE(Orthant<2>::null() != Orthant<2>::null());
}
TEST_CASE("Orthant<3> !=", "[Orthant]")
{
  CHECK_FALSE(Orthant<3>::bsw() != Orthant<3>::bsw());
  CHECK(Orthant<3>::bsw() != Orthant<3>::bse());
  CHECK(Orthant<3>::bsw() != Orthant<3>::bnw());
  CHECK(Orthant<3>::bsw() != Orthant<3>::bne());
  CHECK(Orthant<3>::bsw() != Orthant<3>::tsw());
  CHECK(Orthant<3>::bsw() != Orthant<3>::tse());
  CHECK(Orthant<3>::bsw() != Orthant<3>::tnw());
  CHECK(Orthant<3>::bsw() != Orthant<3>::tne());
  CHECK(Orthant<3>::bsw() != Orthant<3>::null());

  CHECK(Orthant<3>::bse() != Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::bse() != Orthant<3>::bse());
  CHECK(Orthant<3>::bse() != Orthant<3>::bnw());
  CHECK(Orthant<3>::bse() != Orthant<3>::bne());
  CHECK(Orthant<3>::bse() != Orthant<3>::tsw());
  CHECK(Orthant<3>::bse() != Orthant<3>::tse());
  CHECK(Orthant<3>::bse() != Orthant<3>::tnw());
  CHECK(Orthant<3>::bse() != Orthant<3>::tne());
  CHECK(Orthant<3>::bse() != Orthant<3>::null());

  CHECK(Orthant<3>::bnw() != Orthant<3>::bsw());
  CHECK(Orthant<3>::bnw() != Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::bnw() != Orthant<3>::bnw());
  CHECK(Orthant<3>::bnw() != Orthant<3>::bne());
  CHECK(Orthant<3>::bnw() != Orthant<3>::tsw());
  CHECK(Orthant<3>::bnw() != Orthant<3>::tse());
  CHECK(Orthant<3>::bnw() != Orthant<3>::tnw());
  CHECK(Orthant<3>::bnw() != Orthant<3>::tne());
  CHECK(Orthant<3>::bnw() != Orthant<3>::null());

  CHECK(Orthant<3>::bne() != Orthant<3>::bsw());
  CHECK(Orthant<3>::bne() != Orthant<3>::bse());
  CHECK(Orthant<3>::bne() != Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::bne() != Orthant<3>::bne());
  CHECK(Orthant<3>::bne() != Orthant<3>::tsw());
  CHECK(Orthant<3>::bne() != Orthant<3>::tse());
  CHECK(Orthant<3>::bne() != Orthant<3>::tnw());
  CHECK(Orthant<3>::bne() != Orthant<3>::tne());
  CHECK(Orthant<3>::bne() != Orthant<3>::null());

  CHECK(Orthant<3>::tsw() != Orthant<3>::bsw());
  CHECK(Orthant<3>::tsw() != Orthant<3>::bse());
  CHECK(Orthant<3>::tsw() != Orthant<3>::bnw());
  CHECK(Orthant<3>::tsw() != Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::tsw() != Orthant<3>::tsw());
  CHECK(Orthant<3>::tsw() != Orthant<3>::tse());
  CHECK(Orthant<3>::tsw() != Orthant<3>::tnw());
  CHECK(Orthant<3>::tsw() != Orthant<3>::tne());
  CHECK(Orthant<3>::tsw() != Orthant<3>::null());

  CHECK(Orthant<3>::tse() != Orthant<3>::bsw());
  CHECK(Orthant<3>::tse() != Orthant<3>::bse());
  CHECK(Orthant<3>::tse() != Orthant<3>::bnw());
  CHECK(Orthant<3>::tse() != Orthant<3>::bne());
  CHECK(Orthant<3>::tse() != Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::tse() != Orthant<3>::tse());
  CHECK(Orthant<3>::tse() != Orthant<3>::tnw());
  CHECK(Orthant<3>::tse() != Orthant<3>::tne());
  CHECK(Orthant<3>::tse() != Orthant<3>::null());

  CHECK(Orthant<3>::tnw() != Orthant<3>::bsw());
  CHECK(Orthant<3>::tnw() != Orthant<3>::bse());
  CHECK(Orthant<3>::tnw() != Orthant<3>::bnw());
  CHECK(Orthant<3>::tnw() != Orthant<3>::bne());
  CHECK(Orthant<3>::tnw() != Orthant<3>::tsw());
  CHECK(Orthant<3>::tnw() != Orthant<3>::tse());
  CHECK_FALSE(Orthant<3>::tnw() != Orthant<3>::tnw());
  CHECK(Orthant<3>::tnw() != Orthant<3>::tne());
  CHECK(Orthant<3>::tnw() != Orthant<3>::null());

  CHECK(Orthant<3>::tne() != Orthant<3>::bsw());
  CHECK(Orthant<3>::tne() != Orthant<3>::bse());
  CHECK(Orthant<3>::tne() != Orthant<3>::bnw());
  CHECK(Orthant<3>::tne() != Orthant<3>::bne());
  CHECK(Orthant<3>::tne() != Orthant<3>::tsw());
  CHECK(Orthant<3>::tne() != Orthant<3>::tse());
  CHECK(Orthant<3>::tne() != Orthant<3>::tnw());
  CHECK_FALSE(Orthant<3>::tne() != Orthant<3>::tne());
  CHECK(Orthant<3>::tne() != Orthant<3>::null());

  CHECK(Orthant<3>::null() != Orthant<3>::bsw());
  CHECK(Orthant<3>::null() != Orthant<3>::bse());
  CHECK(Orthant<3>::null() != Orthant<3>::bnw());
  CHECK(Orthant<3>::null() != Orthant<3>::bne());
  CHECK(Orthant<3>::null() != Orthant<3>::tsw());
  CHECK(Orthant<3>::null() != Orthant<3>::tse());
  CHECK(Orthant<3>::null() != Orthant<3>::tnw());
  CHECK(Orthant<3>::null() != Orthant<3>::tne());
  CHECK_FALSE(Orthant<3>::null() != Orthant<3>::null());
}
TEST_CASE("Orthant<1> <", "[Orthant]")
{
  CHECK_FALSE(Orthant<1>::lower() < Orthant<1>::lower());
  CHECK(Orthant<1>::lower() < Orthant<1>::upper());
  CHECK(Orthant<1>::lower() < Orthant<1>::null());

  CHECK_FALSE(Orthant<1>::upper() < Orthant<1>::lower());
  CHECK_FALSE(Orthant<1>::upper() < Orthant<1>::upper());
  CHECK(Orthant<1>::upper() < Orthant<1>::null());

  CHECK_FALSE(Orthant<1>::null() < Orthant<1>::lower());
  CHECK_FALSE(Orthant<1>::null() < Orthant<1>::upper());
  CHECK_FALSE(Orthant<1>::null() < Orthant<1>::null());
}
TEST_CASE("Orthant<2> <", "[Orthant]")
{
  CHECK_FALSE(Orthant<2>::sw() < Orthant<2>::sw());
  CHECK(Orthant<2>::sw() < Orthant<2>::se());
  CHECK(Orthant<2>::sw() < Orthant<2>::nw());
  CHECK(Orthant<2>::sw() < Orthant<2>::ne());
  CHECK(Orthant<2>::sw() < Orthant<2>::null());

  CHECK_FALSE(Orthant<2>::se() < Orthant<2>::sw());
  CHECK_FALSE(Orthant<2>::se() < Orthant<2>::se());
  CHECK(Orthant<2>::se() < Orthant<2>::nw());
  CHECK(Orthant<2>::se() < Orthant<2>::ne());
  CHECK(Orthant<2>::se() < Orthant<2>::null());

  CHECK_FALSE(Orthant<2>::nw() < Orthant<2>::sw());
  CHECK_FALSE(Orthant<2>::nw() < Orthant<2>::se());
  CHECK_FALSE(Orthant<2>::nw() < Orthant<2>::nw());
  CHECK(Orthant<2>::nw() < Orthant<2>::ne());
  CHECK(Orthant<2>::nw() < Orthant<2>::null());

  CHECK_FALSE(Orthant<2>::ne() < Orthant<2>::sw());
  CHECK_FALSE(Orthant<2>::ne() < Orthant<2>::se());
  CHECK_FALSE(Orthant<2>::ne() < Orthant<2>::nw());
  CHECK_FALSE(Orthant<2>::ne() < Orthant<2>::ne());
  CHECK(Orthant<2>::ne() < Orthant<2>::null());

  CHECK_FALSE(Orthant<2>::null() < Orthant<2>::sw());
  CHECK_FALSE(Orthant<2>::null() < Orthant<2>::se());
  CHECK_FALSE(Orthant<2>::null() < Orthant<2>::nw());
  CHECK_FALSE(Orthant<2>::null() < Orthant<2>::ne());
  CHECK_FALSE(Orthant<2>::null() < Orthant<2>::null());
}
TEST_CASE("Orthant<3> <", "[Orthant]")
{
  CHECK_FALSE(Orthant<3>::bsw() < Orthant<3>::bsw());
  CHECK(Orthant<3>::bsw() < Orthant<3>::bse());
  CHECK(Orthant<3>::bsw() < Orthant<3>::bnw());
  CHECK(Orthant<3>::bsw() < Orthant<3>::bne());
  CHECK(Orthant<3>::bsw() < Orthant<3>::tsw());
  CHECK(Orthant<3>::bsw() < Orthant<3>::tse());
  CHECK(Orthant<3>::bsw() < Orthant<3>::tnw());
  CHECK(Orthant<3>::bsw() < Orthant<3>::tne());
  CHECK(Orthant<3>::bsw() < Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::bse() < Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::bse() < Orthant<3>::bse());
  CHECK(Orthant<3>::bse() < Orthant<3>::bnw());
  CHECK(Orthant<3>::bse() < Orthant<3>::bne());
  CHECK(Orthant<3>::bse() < Orthant<3>::tsw());
  CHECK(Orthant<3>::bse() < Orthant<3>::tse());
  CHECK(Orthant<3>::bse() < Orthant<3>::tnw());
  CHECK(Orthant<3>::bse() < Orthant<3>::tne());
  CHECK(Orthant<3>::bse() < Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::bnw() < Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::bnw() < Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::bnw() < Orthant<3>::bnw());
  CHECK(Orthant<3>::bnw() < Orthant<3>::bne());
  CHECK(Orthant<3>::bnw() < Orthant<3>::tsw());
  CHECK(Orthant<3>::bnw() < Orthant<3>::tse());
  CHECK(Orthant<3>::bnw() < Orthant<3>::tnw());
  CHECK(Orthant<3>::bnw() < Orthant<3>::tne());
  CHECK(Orthant<3>::bnw() < Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::bne() < Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::bne() < Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::bne() < Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::bne() < Orthant<3>::bne());
  CHECK(Orthant<3>::bne() < Orthant<3>::tsw());
  CHECK(Orthant<3>::bne() < Orthant<3>::tse());
  CHECK(Orthant<3>::bne() < Orthant<3>::tnw());
  CHECK(Orthant<3>::bne() < Orthant<3>::tne());
  CHECK(Orthant<3>::bne() < Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::tsw() < Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::tsw() < Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::tsw() < Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::tsw() < Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::tsw() < Orthant<3>::tsw());
  CHECK(Orthant<3>::tsw() < Orthant<3>::tse());
  CHECK(Orthant<3>::tsw() < Orthant<3>::tnw());
  CHECK(Orthant<3>::tsw() < Orthant<3>::tne());
  CHECK(Orthant<3>::tsw() < Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::tse() < Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::tse() < Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::tse() < Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::tse() < Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::tse() < Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::tse() < Orthant<3>::tse());
  CHECK(Orthant<3>::tse() < Orthant<3>::tnw());
  CHECK(Orthant<3>::tse() < Orthant<3>::tne());
  CHECK(Orthant<3>::tse() < Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::tnw() < Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::tnw() < Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::tnw() < Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::tnw() < Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::tnw() < Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::tnw() < Orthant<3>::tse());
  CHECK_FALSE(Orthant<3>::tnw() < Orthant<3>::tnw());
  CHECK(Orthant<3>::tnw() < Orthant<3>::tne());
  CHECK(Orthant<3>::tnw() < Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::tne() < Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::tne() < Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::tne() < Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::tne() < Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::tne() < Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::tne() < Orthant<3>::tse());
  CHECK_FALSE(Orthant<3>::tne() < Orthant<3>::tnw());
  CHECK_FALSE(Orthant<3>::tne() < Orthant<3>::tne());
  CHECK(Orthant<3>::tne() < Orthant<3>::null());

  CHECK_FALSE(Orthant<3>::null() < Orthant<3>::bsw());
  CHECK_FALSE(Orthant<3>::null() < Orthant<3>::bse());
  CHECK_FALSE(Orthant<3>::null() < Orthant<3>::bnw());
  CHECK_FALSE(Orthant<3>::null() < Orthant<3>::bne());
  CHECK_FALSE(Orthant<3>::null() < Orthant<3>::tsw());
  CHECK_FALSE(Orthant<3>::null() < Orthant<3>::tse());
  CHECK_FALSE(Orthant<3>::null() < Orthant<3>::tnw());
  CHECK_FALSE(Orthant<3>::null() < Orthant<3>::tne());
  CHECK_FALSE(Orthant<3>::null() < Orthant<3>::null());
}
TEST_CASE("Test ostream for Orthant<0>", "[Orthant]")
{
  stringstream ss;
  ss << Orthant<0>::null();
  CHECK(ss.str() == "Orthant<0>::null()");
  ss.str("");
  ss << Orthant<0>(13);
  CHECK(ss.str() == "Orthant<0> invalid value: 13");
}
TEST_CASE("Test ostream for Orthant<1>", "[Orthant]")
{
  stringstream ss;
  ss << Orthant<1>::lower();
  CHECK(ss.str() == "Orthant<1>::lower()");
  ss.str("");
  ss << Orthant<1>::upper();
  CHECK(ss.str() == "Orthant<1>::upper()");
  ss.str("");
  ss << Orthant<1>::null();
  CHECK(ss.str() == "Orthant<1>::null()");
  ss.str("");
  ss << Orthant<1>(13);
  CHECK(ss.str() == "Orthant<1> invalid value: 13");
}
TEST_CASE("Test ostream for Orthant<2>", "[Orthant]")
{
  stringstream ss;
  ss << Orthant<2>::sw();
  CHECK(ss.str() == "Orthant<2>::sw()");
  ss.str("");
  ss << Orthant<2>::se();
  CHECK(ss.str() == "Orthant<2>::se()");
  ss.str("");
  ss << Orthant<2>::nw();
  CHECK(ss.str() == "Orthant<2>::nw()");
  ss.str("");
  ss << Orthant<2>::ne();
  CHECK(ss.str() == "Orthant<2>::ne()");
  ss.str("");
  ss << Orthant<2>::null();
  CHECK(ss.str() == "Orthant<2>::null()");
  ss.str("");
  ss << Orthant<2>(13);
  CHECK(ss.str() == "Orthant<2> invalid value: 13");
}
TEST_CASE("Test ostream for Orthant<3>", "[Orthant]")
{
  stringstream ss;
  ss << Orthant<3>::bsw();
  CHECK(ss.str() == "Orthant<3>::bsw()");
  ss.str("");
  ss << Orthant<3>::bse();
  CHECK(ss.str() == "Orthant<3>::bse()");
  ss.str("");
  ss << Orthant<3>::bnw();
  CHECK(ss.str() == "Orthant<3>::bnw()");
  ss.str("");
  ss << Orthant<3>::bne();
  CHECK(ss.str() == "Orthant<3>::bne()");
  ss.str("");
  ss << Orthant<3>::tsw();
  CHECK(ss.str() == "Orthant<3>::tsw()");
  ss.str("");
  ss << Orthant<3>::tse();
  CHECK(ss.str() == "Orthant<3>::tse()");
  ss.str("");
  ss << Orthant<3>::tnw();
  CHECK(ss.str() == "Orthant<3>::tnw()");
  ss.str("");
  ss << Orthant<3>::tne();
  CHECK(ss.str() == "Orthant<3>::tne()");
  ss.str("");
  ss << Orthant<3>::null();
  CHECK(ss.str() == "Orthant<3>::null()");
  ss.str("");
  ss << Orthant<3>(13);
  CHECK(ss.str() == "Orthant<3> invalid value: 13");
}
TEST_CASE("Test iterator for Orthant<0>", "[Orthant]")
{
  auto iter = Orthant<0>::getValues().begin();
  CHECK(*iter == Orthant<0>(0));
  ++iter;
  CHECK(*iter == Orthant<0>::null());
  CHECK(iter == Orthant<0>::getValues().end());
}
TEST_CASE("Test iterator for Orthant<1>", "[Orthant]")
{
  auto iter = Orthant<1>::getValues().begin();
  CHECK(iter == Orthant<1>::getValues().begin());
  CHECK(iter != Orthant<1>::getValues().end());
  CHECK(*iter == Orthant<1>::lower());
  ++iter;
  CHECK(iter->getIndex() == 1);
  CHECK(*iter == Orthant<1>::upper());
  ++iter;
  CHECK(*iter == Orthant<1>::null());
  CHECK(iter == Orthant<1>::getValues().end());
}
TEST_CASE("Test iterator for Orthant<2>", "[Orthant]")
{
  auto iter = Orthant<2>::getValues().begin();
  CHECK(iter == Orthant<2>::getValues().begin());
  CHECK(iter != Orthant<2>::getValues().end());
  CHECK(*iter == Orthant<2>::sw());
  ++iter;
  CHECK(iter->getIndex() == 1);
  CHECK(*iter == Orthant<2>::se());
  ++iter;
  CHECK(*iter == Orthant<2>::nw());
  ++iter;
  CHECK(*iter == Orthant<2>::ne());
  ++iter;
  CHECK(*iter == Orthant<2>::null());
  CHECK(iter == Orthant<2>::getValues().end());
}
TEST_CASE("Test iterator for Orthant<3>", "[Orthant]")
{
  auto iter = Orthant<3>::getValues().begin();
  CHECK(iter == Orthant<3>::getValues().begin());
  CHECK(iter != Orthant<3>::getValues().end());
  CHECK(*iter == Orthant<3>::bsw());
  ++iter;
  CHECK(iter->getIndex() == 1);
  CHECK(*iter == Orthant<3>::bse());
  ++iter;
  CHECK(*iter == Orthant<3>::bnw());
  ++iter;
  CHECK(*iter == Orthant<3>::bne());
  ++iter;
  CHECK(*iter == Orthant<3>::tsw());
  ++iter;
  CHECK(*iter == Orthant<3>::tse());
  ++iter;
  CHECK(*iter == Orthant<3>::tnw());
  ++iter;
  CHECK(*iter == Orthant<3>::tne());
  ++iter;
  CHECK(*iter == Orthant<3>::null());
  CHECK(iter == Orthant<3>::getValues().end());
}
TEST_CASE("Test from_json for Orthant<0>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  CHECK(j["null"].get<Orthant<0>>() == Orthant<0>::null());
}
TEST_CASE("Test from_json for Orthant<1>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["lower"] = "LOWER";
  j["upper"] = "UPPER";
  CHECK(j["null"].get<Orthant<1>>() == Orthant<1>::null());
  CHECK(j["lower"].get<Orthant<1>>() == Orthant<1>::lower());
  CHECK(j["upper"].get<Orthant<1>>() == Orthant<1>::upper());
}
TEST_CASE("Test from_json for Orthant<2>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["sw"] = "SW";
  j["se"] = "SE";
  j["nw"] = "NW";
  j["ne"] = "NE";
  CHECK(j["null"].get<Orthant<2>>() == Orthant<2>::null());
  CHECK(j["sw"].get<Orthant<2>>() == Orthant<2>::sw());
  CHECK(j["se"].get<Orthant<2>>() == Orthant<2>::se());
  CHECK(j["nw"].get<Orthant<2>>() == Orthant<2>::nw());
  CHECK(j["ne"].get<Orthant<2>>() == Orthant<2>::ne());
}
TEST_CASE("Test from_json for Orthant<3>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = nullptr;
  j["bsw"] = "BSW";
  j["bse"] = "BSE";
  j["bnw"] = "BNW";
  j["bne"] = "BNE";
  j["tsw"] = "TSW";
  j["tse"] = "TSE";
  j["tnw"] = "TNW";
  j["tne"] = "TNE";
  CHECK(j["null"].get<Orthant<3>>() == Orthant<3>::null());
  CHECK(j["bsw"].get<Orthant<3>>() == Orthant<3>::bsw());
  CHECK(j["bse"].get<Orthant<3>>() == Orthant<3>::bse());
  CHECK(j["bnw"].get<Orthant<3>>() == Orthant<3>::bnw());
  CHECK(j["bne"].get<Orthant<3>>() == Orthant<3>::bne());
  CHECK(j["tsw"].get<Orthant<3>>() == Orthant<3>::tsw());
  CHECK(j["tse"].get<Orthant<3>>() == Orthant<3>::tse());
  CHECK(j["tnw"].get<Orthant<3>>() == Orthant<3>::tnw());
  CHECK(j["tne"].get<Orthant<3>>() == Orthant<3>::tne());
}
TEST_CASE("Test to_json for Orthant<0>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = Orthant<0>::null();
  CHECK(j["null"] == nullptr);
}
TEST_CASE("Test to_json for Orthant<1>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = Orthant<1>::null();
  j["lower"] = Orthant<1>::lower();
  j["upper"] = Orthant<1>::upper();
  CHECK(j["null"] == nullptr);
  CHECK(j["lower"] == "LOWER");
  CHECK(j["upper"] == "UPPER");
}
TEST_CASE("Test to_json for Orthant<2>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = Orthant<2>::null();
  j["sw"] = Orthant<2>::sw();
  j["se"] = Orthant<2>::se();
  j["nw"] = Orthant<2>::nw();
  j["ne"] = Orthant<2>::ne();
  CHECK(j["null"] == nullptr);
  CHECK(j["sw"] == "SW");
  CHECK(j["se"] == "SE");
  CHECK(j["nw"] == "NW");
  CHECK(j["ne"] == "NE");
}
TEST_CASE("Test to_json for Orthant<3>", "[Orthant]")
{
  nlohmann::json j;
  j["null"] = Orthant<3>::null();
  j["bsw"] = Orthant<3>::bsw();
  j["bse"] = Orthant<3>::bse();
  j["bnw"] = Orthant<3>::bnw();
  j["bne"] = Orthant<3>::bne();
  j["tsw"] = Orthant<3>::tsw();
  j["tse"] = Orthant<3>::tse();
  j["tnw"] = Orthant<3>::tnw();
  j["tne"] = Orthant<3>::tne();
  CHECK(j["null"] == nullptr);
  CHECK(j["bsw"] == "BSW");
  CHECK(j["bse"] == "BSE");
  CHECK(j["bnw"] == "BNW");
  CHECK(j["bne"] == "BNE");
  CHECK(j["tsw"] == "TSW");
  CHECK(j["tse"] == "TSE");
  CHECK(j["tnw"] == "TNW");
  CHECK(j["tne"] == "TNE");
}