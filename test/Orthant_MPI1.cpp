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

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Orthant<0> unsigned char constructor works")
{
  Orthant<0> o(13);
  CHECK_EQ(o.getIndex(), 13);
}
TEST_CASE("Orthant<1> unsigned char constructor works")
{
  Orthant<1> o(13);
  CHECK_EQ(o.getIndex(), 13);
}
TEST_CASE("Orthant<2> unsigned char constructor works")
{
  Orthant<2> o(13);
  CHECK_EQ(o.getIndex(), 13);
}
TEST_CASE("Orthant<3> unsigned char constructor works")
{
  Orthant<3> o(13);
  CHECK_EQ(o.getIndex(), 13);
}
TEST_CASE("Orthant<0> Default constructor works")
{
  Orthant<0> o;
  CHECK_EQ(o, Orthant<0>::null());
}
TEST_CASE("Orthant<1> Default constructor works")
{
  Orthant<1> o;
  CHECK_EQ(o, Orthant<1>::null());
}
TEST_CASE("Orthant<2> Default constructor works")
{
  Orthant<2> o;
  CHECK_EQ(o, Orthant<2>::null());
}
TEST_CASE("Orthant<3> Default constructor works")
{
  Orthant<3> o;
  CHECK_EQ(o, Orthant<3>::null());
}
TEST_CASE("Orthant<0> named constructors give expected index values")
{
  CHECK_EQ(Orthant<0>::null().getIndex(), 1);
}
TEST_CASE("Orthant<1> named constructors give expected index values")
{
  CHECK_EQ(Orthant<1>::lower().getIndex(), 0);
  CHECK_EQ(Orthant<1>::upper().getIndex(), 1);
  CHECK_EQ(Orthant<1>::null().getIndex(), 2);
}
TEST_CASE("Orthant<2> named constructors give expected index values")
{
  CHECK_EQ(Orthant<2>::sw().getIndex(), 0);
  CHECK_EQ(Orthant<2>::se().getIndex(), 1);
  CHECK_EQ(Orthant<2>::nw().getIndex(), 2);
  CHECK_EQ(Orthant<2>::ne().getIndex(), 3);
  CHECK_EQ(Orthant<2>::null().getIndex(), 4);
}
TEST_CASE("Orthant<3> named constructors give expected index values")
{
  CHECK_EQ(Orthant<3>::bsw().getIndex(), 0);
  CHECK_EQ(Orthant<3>::bse().getIndex(), 1);
  CHECK_EQ(Orthant<3>::bnw().getIndex(), 2);
  CHECK_EQ(Orthant<3>::bne().getIndex(), 3);
  CHECK_EQ(Orthant<3>::tsw().getIndex(), 4);
  CHECK_EQ(Orthant<3>::tse().getIndex(), 5);
  CHECK_EQ(Orthant<3>::tnw().getIndex(), 6);
  CHECK_EQ(Orthant<3>::tne().getIndex(), 7);
  CHECK_EQ(Orthant<3>::null().getIndex(), 8);
}
TEST_CASE("Orthant<1> getNbrOnSide is as expected")
{
  CHECK_EQ(Orthant<1>::lower().getNbrOnSide(Side<1>::west()), Orthant<1>::upper());
  CHECK_EQ(Orthant<1>::lower().getNbrOnSide(Side<1>::east()), Orthant<1>::upper());

  CHECK_EQ(Orthant<1>::upper().getNbrOnSide(Side<1>::west()), Orthant<1>::lower());
  CHECK_EQ(Orthant<1>::upper().getNbrOnSide(Side<1>::east()), Orthant<1>::lower());
}
TEST_CASE("Orthant<2> getNbrOnSide is as expected")
{
  CHECK_EQ(Orthant<2>::sw().getNbrOnSide(Side<2>::west()), Orthant<2>::se());
  CHECK_EQ(Orthant<2>::sw().getNbrOnSide(Side<2>::east()), Orthant<2>::se());
  CHECK_EQ(Orthant<2>::sw().getNbrOnSide(Side<2>::south()), Orthant<2>::nw());
  CHECK_EQ(Orthant<2>::sw().getNbrOnSide(Side<2>::north()), Orthant<2>::nw());

  CHECK_EQ(Orthant<2>::se().getNbrOnSide(Side<2>::west()), Orthant<2>::sw());
  CHECK_EQ(Orthant<2>::se().getNbrOnSide(Side<2>::east()), Orthant<2>::sw());
  CHECK_EQ(Orthant<2>::se().getNbrOnSide(Side<2>::south()), Orthant<2>::ne());
  CHECK_EQ(Orthant<2>::se().getNbrOnSide(Side<2>::north()), Orthant<2>::ne());

  CHECK_EQ(Orthant<2>::nw().getNbrOnSide(Side<2>::west()), Orthant<2>::ne());
  CHECK_EQ(Orthant<2>::nw().getNbrOnSide(Side<2>::east()), Orthant<2>::ne());
  CHECK_EQ(Orthant<2>::nw().getNbrOnSide(Side<2>::south()), Orthant<2>::sw());
  CHECK_EQ(Orthant<2>::nw().getNbrOnSide(Side<2>::north()), Orthant<2>::sw());

  CHECK_EQ(Orthant<2>::ne().getNbrOnSide(Side<2>::west()), Orthant<2>::nw());
  CHECK_EQ(Orthant<2>::ne().getNbrOnSide(Side<2>::east()), Orthant<2>::nw());
  CHECK_EQ(Orthant<2>::ne().getNbrOnSide(Side<2>::south()), Orthant<2>::se());
  CHECK_EQ(Orthant<2>::ne().getNbrOnSide(Side<2>::north()), Orthant<2>::se());
}
TEST_CASE("Orthant<3> getNbrOnSide is as expected")
{
  CHECK_EQ(Orthant<3>::bsw().getNbrOnSide(Side<3>::west()), Orthant<3>::bse());
  CHECK_EQ(Orthant<3>::bsw().getNbrOnSide(Side<3>::east()), Orthant<3>::bse());
  CHECK_EQ(Orthant<3>::bsw().getNbrOnSide(Side<3>::south()), Orthant<3>::bnw());
  CHECK_EQ(Orthant<3>::bsw().getNbrOnSide(Side<3>::north()), Orthant<3>::bnw());
  CHECK_EQ(Orthant<3>::bsw().getNbrOnSide(Side<3>::bottom()), Orthant<3>::tsw());
  CHECK_EQ(Orthant<3>::bsw().getNbrOnSide(Side<3>::top()), Orthant<3>::tsw());

  CHECK_EQ(Orthant<3>::bse().getNbrOnSide(Side<3>::west()), Orthant<3>::bsw());
  CHECK_EQ(Orthant<3>::bse().getNbrOnSide(Side<3>::east()), Orthant<3>::bsw());
  CHECK_EQ(Orthant<3>::bse().getNbrOnSide(Side<3>::south()), Orthant<3>::bne());
  CHECK_EQ(Orthant<3>::bse().getNbrOnSide(Side<3>::north()), Orthant<3>::bne());
  CHECK_EQ(Orthant<3>::bse().getNbrOnSide(Side<3>::bottom()), Orthant<3>::tse());
  CHECK_EQ(Orthant<3>::bse().getNbrOnSide(Side<3>::top()), Orthant<3>::tse());

  CHECK_EQ(Orthant<3>::bnw().getNbrOnSide(Side<3>::west()), Orthant<3>::bne());
  CHECK_EQ(Orthant<3>::bnw().getNbrOnSide(Side<3>::east()), Orthant<3>::bne());
  CHECK_EQ(Orthant<3>::bnw().getNbrOnSide(Side<3>::south()), Orthant<3>::bsw());
  CHECK_EQ(Orthant<3>::bnw().getNbrOnSide(Side<3>::north()), Orthant<3>::bsw());
  CHECK_EQ(Orthant<3>::bnw().getNbrOnSide(Side<3>::bottom()), Orthant<3>::tnw());
  CHECK_EQ(Orthant<3>::bnw().getNbrOnSide(Side<3>::top()), Orthant<3>::tnw());

  CHECK_EQ(Orthant<3>::bne().getNbrOnSide(Side<3>::west()), Orthant<3>::bnw());
  CHECK_EQ(Orthant<3>::bne().getNbrOnSide(Side<3>::east()), Orthant<3>::bnw());
  CHECK_EQ(Orthant<3>::bne().getNbrOnSide(Side<3>::south()), Orthant<3>::bse());
  CHECK_EQ(Orthant<3>::bne().getNbrOnSide(Side<3>::north()), Orthant<3>::bse());
  CHECK_EQ(Orthant<3>::bne().getNbrOnSide(Side<3>::bottom()), Orthant<3>::tne());
  CHECK_EQ(Orthant<3>::bne().getNbrOnSide(Side<3>::top()), Orthant<3>::tne());

  CHECK_EQ(Orthant<3>::tsw().getNbrOnSide(Side<3>::west()), Orthant<3>::tse());
  CHECK_EQ(Orthant<3>::tsw().getNbrOnSide(Side<3>::east()), Orthant<3>::tse());
  CHECK_EQ(Orthant<3>::tsw().getNbrOnSide(Side<3>::south()), Orthant<3>::tnw());
  CHECK_EQ(Orthant<3>::tsw().getNbrOnSide(Side<3>::north()), Orthant<3>::tnw());
  CHECK_EQ(Orthant<3>::tsw().getNbrOnSide(Side<3>::bottom()), Orthant<3>::bsw());
  CHECK_EQ(Orthant<3>::tsw().getNbrOnSide(Side<3>::top()), Orthant<3>::bsw());

  CHECK_EQ(Orthant<3>::tse().getNbrOnSide(Side<3>::east()), Orthant<3>::tsw());
  CHECK_EQ(Orthant<3>::tse().getNbrOnSide(Side<3>::west()), Orthant<3>::tsw());
  CHECK_EQ(Orthant<3>::tse().getNbrOnSide(Side<3>::south()), Orthant<3>::tne());
  CHECK_EQ(Orthant<3>::tse().getNbrOnSide(Side<3>::north()), Orthant<3>::tne());
  CHECK_EQ(Orthant<3>::tse().getNbrOnSide(Side<3>::bottom()), Orthant<3>::bse());
  CHECK_EQ(Orthant<3>::tse().getNbrOnSide(Side<3>::top()), Orthant<3>::bse());

  CHECK_EQ(Orthant<3>::tnw().getNbrOnSide(Side<3>::west()), Orthant<3>::tne());
  CHECK_EQ(Orthant<3>::tnw().getNbrOnSide(Side<3>::east()), Orthant<3>::tne());
  CHECK_EQ(Orthant<3>::tnw().getNbrOnSide(Side<3>::south()), Orthant<3>::tsw());
  CHECK_EQ(Orthant<3>::tnw().getNbrOnSide(Side<3>::north()), Orthant<3>::tsw());
  CHECK_EQ(Orthant<3>::tnw().getNbrOnSide(Side<3>::bottom()), Orthant<3>::bnw());
  CHECK_EQ(Orthant<3>::tnw().getNbrOnSide(Side<3>::top()), Orthant<3>::bnw());

  CHECK_EQ(Orthant<3>::tne().getNbrOnSide(Side<3>::east()), Orthant<3>::tnw());
  CHECK_EQ(Orthant<3>::tne().getNbrOnSide(Side<3>::west()), Orthant<3>::tnw());
  CHECK_EQ(Orthant<3>::tne().getNbrOnSide(Side<3>::south()), Orthant<3>::tse());
  CHECK_EQ(Orthant<3>::tne().getNbrOnSide(Side<3>::north()), Orthant<3>::tse());
  CHECK_EQ(Orthant<3>::tne().getNbrOnSide(Side<3>::bottom()), Orthant<3>::bne());
  CHECK_EQ(Orthant<3>::tne().getNbrOnSide(Side<3>::top()), Orthant<3>::bne());
}
TEST_CASE("Orthant<1> getInteriorSides is as expected")
{
  {
    auto array = Orthant<1>::lower().getInteriorSides();
    CHECK_EQ(array[0], Side<1>::east());
  }
  {
    auto array = Orthant<1>::upper().getInteriorSides();
    CHECK_EQ(array[0], Side<1>::west());
  }
}
TEST_CASE("Orthant<2> getInteriorSides is as expected")
{
  {
    auto array = Orthant<2>::sw().getInteriorSides();
    CHECK_EQ(array[0], Side<2>::east());
    CHECK_EQ(array[1], Side<2>::north());
  }
  {
    auto array = Orthant<2>::se().getInteriorSides();
    CHECK_EQ(array[0], Side<2>::west());
    CHECK_EQ(array[1], Side<2>::north());
  }
  {
    auto array = Orthant<2>::nw().getInteriorSides();
    CHECK_EQ(array[0], Side<2>::east());
    CHECK_EQ(array[1], Side<2>::south());
  }
  {
    auto array = Orthant<2>::ne().getInteriorSides();
    CHECK_EQ(array[0], Side<2>::west());
    CHECK_EQ(array[1], Side<2>::south());
  }
}
TEST_CASE("Orthant<3> getInteriorSides is as expected")
{
  {
    auto array = Orthant<3>::bsw().getInteriorSides();
    CHECK_EQ(array[0], Side<3>::east());
    CHECK_EQ(array[1], Side<3>::north());
    CHECK_EQ(array[2], Side<3>::top());
  }
  {
    auto array = Orthant<3>::bse().getInteriorSides();
    CHECK_EQ(array[0], Side<3>::west());
    CHECK_EQ(array[1], Side<3>::north());
    CHECK_EQ(array[2], Side<3>::top());
  }
  {
    auto array = Orthant<3>::bnw().getInteriorSides();
    CHECK_EQ(array[0], Side<3>::east());
    CHECK_EQ(array[1], Side<3>::south());
    CHECK_EQ(array[2], Side<3>::top());
  }
  {
    auto array = Orthant<3>::bne().getInteriorSides();
    CHECK_EQ(array[0], Side<3>::west());
    CHECK_EQ(array[1], Side<3>::south());
    CHECK_EQ(array[2], Side<3>::top());
  }
  {
    auto array = Orthant<3>::tsw().getInteriorSides();
    CHECK_EQ(array[0], Side<3>::east());
    CHECK_EQ(array[1], Side<3>::north());
    CHECK_EQ(array[2], Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::tse().getInteriorSides();
    CHECK_EQ(array[0], Side<3>::west());
    CHECK_EQ(array[1], Side<3>::north());
    CHECK_EQ(array[2], Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::tnw().getInteriorSides();
    CHECK_EQ(array[0], Side<3>::east());
    CHECK_EQ(array[1], Side<3>::south());
    CHECK_EQ(array[2], Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::tne().getInteriorSides();
    CHECK_EQ(array[0], Side<3>::west());
    CHECK_EQ(array[1], Side<3>::south());
    CHECK_EQ(array[2], Side<3>::bottom());
  }
}
TEST_CASE("Orthant<1> getExteriorSides is as expected")
{
  {
    auto array = Orthant<1>::lower().getExteriorSides();
    CHECK_EQ(array[0], Side<1>::west());
  }
  {
    auto array = Orthant<1>::upper().getExteriorSides();
    CHECK_EQ(array[0], Side<1>::east());
  }
}
TEST_CASE("Orthant<2> getExteriorSides is as expected")
{
  {
    auto array = Orthant<2>::sw().getExteriorSides();
    CHECK_EQ(array[0], Side<2>::west());
    CHECK_EQ(array[1], Side<2>::south());
  }
  {
    auto array = Orthant<2>::se().getExteriorSides();
    CHECK_EQ(array[0], Side<2>::east());
    CHECK_EQ(array[1], Side<2>::south());
  }
  {
    auto array = Orthant<2>::nw().getExteriorSides();
    CHECK_EQ(array[0], Side<2>::west());
    CHECK_EQ(array[1], Side<2>::north());
  }
  {
    auto array = Orthant<2>::ne().getExteriorSides();
    CHECK_EQ(array[0], Side<2>::east());
    CHECK_EQ(array[1], Side<2>::north());
  }
}
TEST_CASE("Orthant<3> getExteriorSides is as expected")
{
  {
    auto array = Orthant<3>::bsw().getExteriorSides();
    CHECK_EQ(array[0], Side<3>::west());
    CHECK_EQ(array[1], Side<3>::south());
    CHECK_EQ(array[2], Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::bse().getExteriorSides();
    CHECK_EQ(array[0], Side<3>::east());
    CHECK_EQ(array[1], Side<3>::south());
    CHECK_EQ(array[2], Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::bnw().getExteriorSides();
    CHECK_EQ(array[0], Side<3>::west());
    CHECK_EQ(array[1], Side<3>::north());
    CHECK_EQ(array[2], Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::bne().getExteriorSides();
    CHECK_EQ(array[0], Side<3>::east());
    CHECK_EQ(array[1], Side<3>::north());
    CHECK_EQ(array[2], Side<3>::bottom());
  }
  {
    auto array = Orthant<3>::tsw().getExteriorSides();
    CHECK_EQ(array[0], Side<3>::west());
    CHECK_EQ(array[1], Side<3>::south());
    CHECK_EQ(array[2], Side<3>::top());
  }
  {
    auto array = Orthant<3>::tse().getExteriorSides();
    CHECK_EQ(array[0], Side<3>::east());
    CHECK_EQ(array[1], Side<3>::south());
    CHECK_EQ(array[2], Side<3>::top());
  }
  {
    auto array = Orthant<3>::tnw().getExteriorSides();
    CHECK_EQ(array[0], Side<3>::west());
    CHECK_EQ(array[1], Side<3>::north());
    CHECK_EQ(array[2], Side<3>::top());
  }
  {
    auto array = Orthant<3>::tne().getExteriorSides();
    CHECK_EQ(array[0], Side<3>::east());
    CHECK_EQ(array[1], Side<3>::north());
    CHECK_EQ(array[2], Side<3>::top());
  }
}
TEST_CASE("Orthant<1> isOnSide is as expected")
{
  CHECK_UNARY(Orthant<1>::lower().isOnSide(Side<1>::west()));
  CHECK_UNARY_FALSE(Orthant<1>::lower().isOnSide(Side<1>::east()));

  CHECK_UNARY_FALSE(Orthant<1>::upper().isOnSide(Side<1>::west()));
  CHECK_UNARY(Orthant<1>::upper().isOnSide(Side<1>::east()));
}
TEST_CASE("Orthant<2> isOnSide is as expected")
{
  CHECK_UNARY(Orthant<2>::sw().isOnSide(Side<2>::west()));
  CHECK_UNARY_FALSE(Orthant<2>::sw().isOnSide(Side<2>::east()));
  CHECK_UNARY(Orthant<2>::sw().isOnSide(Side<2>::south()));
  CHECK_UNARY_FALSE(Orthant<2>::sw().isOnSide(Side<2>::north()));

  CHECK_UNARY_FALSE(Orthant<2>::se().isOnSide(Side<2>::west()));
  CHECK_UNARY(Orthant<2>::se().isOnSide(Side<2>::east()));
  CHECK_UNARY(Orthant<2>::se().isOnSide(Side<2>::south()));
  CHECK_UNARY_FALSE(Orthant<2>::se().isOnSide(Side<2>::north()));

  CHECK_UNARY(Orthant<2>::nw().isOnSide(Side<2>::west()));
  CHECK_UNARY_FALSE(Orthant<2>::nw().isOnSide(Side<2>::east()));
  CHECK_UNARY_FALSE(Orthant<2>::nw().isOnSide(Side<2>::south()));
  CHECK_UNARY(Orthant<2>::nw().isOnSide(Side<2>::north()));

  CHECK_UNARY_FALSE(Orthant<2>::ne().isOnSide(Side<2>::west()));
  CHECK_UNARY(Orthant<2>::ne().isOnSide(Side<2>::east()));
  CHECK_UNARY_FALSE(Orthant<2>::ne().isOnSide(Side<2>::south()));
  CHECK_UNARY(Orthant<2>::ne().isOnSide(Side<2>::north()));
}
TEST_CASE("Orthant<3> isOnSide is as expected")
{
  CHECK_UNARY(Orthant<3>::bsw().isOnSide(Side<3>::west()));
  CHECK_UNARY_FALSE(Orthant<3>::bsw().isOnSide(Side<3>::east()));
  CHECK_UNARY(Orthant<3>::bsw().isOnSide(Side<3>::south()));
  CHECK_UNARY_FALSE(Orthant<3>::bsw().isOnSide(Side<3>::north()));
  CHECK_UNARY(Orthant<3>::bsw().isOnSide(Side<3>::bottom()));
  CHECK_UNARY_FALSE(Orthant<3>::bsw().isOnSide(Side<3>::top()));

  CHECK_UNARY_FALSE(Orthant<3>::bse().isOnSide(Side<3>::west()));
  CHECK_UNARY(Orthant<3>::bse().isOnSide(Side<3>::east()));
  CHECK_UNARY(Orthant<3>::bse().isOnSide(Side<3>::south()));
  CHECK_UNARY_FALSE(Orthant<3>::bse().isOnSide(Side<3>::north()));
  CHECK_UNARY(Orthant<3>::bse().isOnSide(Side<3>::bottom()));
  CHECK_UNARY_FALSE(Orthant<3>::bse().isOnSide(Side<3>::top()));

  CHECK_UNARY(Orthant<3>::bnw().isOnSide(Side<3>::west()));
  CHECK_UNARY_FALSE(Orthant<3>::bnw().isOnSide(Side<3>::east()));
  CHECK_UNARY_FALSE(Orthant<3>::bnw().isOnSide(Side<3>::south()));
  CHECK_UNARY(Orthant<3>::bnw().isOnSide(Side<3>::north()));
  CHECK_UNARY(Orthant<3>::bnw().isOnSide(Side<3>::bottom()));
  CHECK_UNARY_FALSE(Orthant<3>::bnw().isOnSide(Side<3>::top()));

  CHECK_UNARY_FALSE(Orthant<3>::bne().isOnSide(Side<3>::west()));
  CHECK_UNARY(Orthant<3>::bne().isOnSide(Side<3>::east()));
  CHECK_UNARY_FALSE(Orthant<3>::bne().isOnSide(Side<3>::south()));
  CHECK_UNARY(Orthant<3>::bne().isOnSide(Side<3>::north()));
  CHECK_UNARY(Orthant<3>::bne().isOnSide(Side<3>::bottom()));
  CHECK_UNARY_FALSE(Orthant<3>::bne().isOnSide(Side<3>::top()));

  CHECK_UNARY(Orthant<3>::tsw().isOnSide(Side<3>::west()));
  CHECK_UNARY_FALSE(Orthant<3>::tsw().isOnSide(Side<3>::east()));
  CHECK_UNARY(Orthant<3>::tsw().isOnSide(Side<3>::south()));
  CHECK_UNARY_FALSE(Orthant<3>::tsw().isOnSide(Side<3>::north()));
  CHECK_UNARY_FALSE(Orthant<3>::tsw().isOnSide(Side<3>::bottom()));
  CHECK_UNARY(Orthant<3>::tsw().isOnSide(Side<3>::top()));

  CHECK_UNARY_FALSE(Orthant<3>::tse().isOnSide(Side<3>::west()));
  CHECK_UNARY(Orthant<3>::tse().isOnSide(Side<3>::east()));
  CHECK_UNARY(Orthant<3>::tse().isOnSide(Side<3>::south()));
  CHECK_UNARY_FALSE(Orthant<3>::tse().isOnSide(Side<3>::north()));
  CHECK_UNARY_FALSE(Orthant<3>::tse().isOnSide(Side<3>::bottom()));
  CHECK_UNARY(Orthant<3>::tse().isOnSide(Side<3>::top()));

  CHECK_UNARY(Orthant<3>::tnw().isOnSide(Side<3>::west()));
  CHECK_UNARY_FALSE(Orthant<3>::tnw().isOnSide(Side<3>::east()));
  CHECK_UNARY_FALSE(Orthant<3>::tnw().isOnSide(Side<3>::south()));
  CHECK_UNARY(Orthant<3>::tnw().isOnSide(Side<3>::north()));
  CHECK_UNARY_FALSE(Orthant<3>::tnw().isOnSide(Side<3>::bottom()));
  CHECK_UNARY(Orthant<3>::tnw().isOnSide(Side<3>::top()));

  CHECK_UNARY_FALSE(Orthant<3>::tne().isOnSide(Side<3>::west()));
  CHECK_UNARY(Orthant<3>::tne().isOnSide(Side<3>::east()));
  CHECK_UNARY_FALSE(Orthant<3>::tne().isOnSide(Side<3>::south()));
  CHECK_UNARY(Orthant<3>::tne().isOnSide(Side<3>::north()));
  CHECK_UNARY_FALSE(Orthant<3>::tne().isOnSide(Side<3>::bottom()));
  CHECK_UNARY(Orthant<3>::tne().isOnSide(Side<3>::top()));
}
TEST_CASE("Orthant<1> getValuesOnSide is as expected")
{
  {
    std::array<Orthant<1>, 1> values = Orthant<1>::getValuesOnSide(Side<1>::west());
    CHECK_EQ(values[0], Orthant<1>::lower());
  }
  {
    std::array<Orthant<1>, 1> values = Orthant<1>::getValuesOnSide(Side<1>::east());
    CHECK_EQ(values[0], Orthant<1>::upper());
  }
}
TEST_CASE("Orthant<2> getValuesOnSide is as expected")
{
  {
    std::array<Orthant<2>, 2> values = Orthant<2>::getValuesOnSide(Side<2>::west());
    CHECK_EQ(values[0], Orthant<2>::sw());
    CHECK_EQ(values[1], Orthant<2>::nw());
  }

  {
    std::array<Orthant<2>, 2> values = Orthant<2>::getValuesOnSide(Side<2>::east());
    CHECK_EQ(values[0], Orthant<2>::se());
    CHECK_EQ(values[1], Orthant<2>::ne());
  }

  {
    std::array<Orthant<2>, 2> values = Orthant<2>::getValuesOnSide(Side<2>::south());
    CHECK_EQ(values[0], Orthant<2>::sw());
    CHECK_EQ(values[1], Orthant<2>::se());
  }

  {
    std::array<Orthant<2>, 2> values = Orthant<2>::getValuesOnSide(Side<2>::north());
    CHECK_EQ(values[0], Orthant<2>::nw());
    CHECK_EQ(values[1], Orthant<2>::ne());
  }
}
TEST_CASE("Orthant<3> getValuesOnSide is as expected")
{
  {
    std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::west());
    CHECK_EQ(values[0], Orthant<3>::bsw());
    CHECK_EQ(values[1], Orthant<3>::bnw());
    CHECK_EQ(values[2], Orthant<3>::tsw());
    CHECK_EQ(values[3], Orthant<3>::tnw());
  }

  {
    std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::east());
    CHECK_EQ(values[0], Orthant<3>::bse());
    CHECK_EQ(values[1], Orthant<3>::bne());
    CHECK_EQ(values[2], Orthant<3>::tse());
    CHECK_EQ(values[3], Orthant<3>::tne());
  }

  {
    std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::south());
    CHECK_EQ(values[0], Orthant<3>::bsw());
    CHECK_EQ(values[1], Orthant<3>::bse());
    CHECK_EQ(values[2], Orthant<3>::tsw());
    CHECK_EQ(values[3], Orthant<3>::tse());
  }

  {
    std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::north());
    CHECK_EQ(values[0], Orthant<3>::bnw());
    CHECK_EQ(values[1], Orthant<3>::bne());
    CHECK_EQ(values[2], Orthant<3>::tnw());
    CHECK_EQ(values[3], Orthant<3>::tne());
  }

  {
    std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::bottom());
    CHECK_EQ(values[0], Orthant<3>::bsw());
    CHECK_EQ(values[1], Orthant<3>::bse());
    CHECK_EQ(values[2], Orthant<3>::bnw());
    CHECK_EQ(values[3], Orthant<3>::bne());
  }

  {
    std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::top());
    CHECK_EQ(values[0], Orthant<3>::tsw());
    CHECK_EQ(values[1], Orthant<3>::tse());
    CHECK_EQ(values[2], Orthant<3>::tnw());
    CHECK_EQ(values[3], Orthant<3>::tne());
  }
}
TEST_CASE("Orthant<1> collapseOnAxis is as expected")
{
  CHECK_EQ(Orthant<1>::lower().collapseOnAxis(0), Orthant<0>(0));
  CHECK_EQ(Orthant<1>::upper().collapseOnAxis(0), Orthant<0>(0));
}
TEST_CASE("Orthant<2> collapseOnAxis is as expected")
{
  CHECK_EQ(Orthant<2>::sw().collapseOnAxis(0), Orthant<1>::lower());
  CHECK_EQ(Orthant<2>::sw().collapseOnAxis(1), Orthant<1>::lower());

  CHECK_EQ(Orthant<2>::se().collapseOnAxis(0), Orthant<1>::lower());
  CHECK_EQ(Orthant<2>::se().collapseOnAxis(1), Orthant<1>::upper());

  CHECK_EQ(Orthant<2>::nw().collapseOnAxis(0), Orthant<1>::upper());
  CHECK_EQ(Orthant<2>::nw().collapseOnAxis(1), Orthant<1>::lower());

  CHECK_EQ(Orthant<2>::ne().collapseOnAxis(0), Orthant<1>::upper());
  CHECK_EQ(Orthant<2>::ne().collapseOnAxis(1), Orthant<1>::upper());
}
TEST_CASE("Orthant<3> collapseOnAxis is as expected")
{
  CHECK_EQ(Orthant<3>::bsw().collapseOnAxis(0), Orthant<2>::sw());
  CHECK_EQ(Orthant<3>::bsw().collapseOnAxis(1), Orthant<2>::sw());
  CHECK_EQ(Orthant<3>::bsw().collapseOnAxis(2), Orthant<2>::sw());

  CHECK_EQ(Orthant<3>::bse().collapseOnAxis(0), Orthant<2>::sw());
  CHECK_EQ(Orthant<3>::bse().collapseOnAxis(1), Orthant<2>::se());
  CHECK_EQ(Orthant<3>::bse().collapseOnAxis(2), Orthant<2>::se());

  CHECK_EQ(Orthant<3>::bnw().collapseOnAxis(0), Orthant<2>::se());
  CHECK_EQ(Orthant<3>::bnw().collapseOnAxis(1), Orthant<2>::sw());
  CHECK_EQ(Orthant<3>::bnw().collapseOnAxis(2), Orthant<2>::nw());

  CHECK_EQ(Orthant<3>::bne().collapseOnAxis(0), Orthant<2>::se());
  CHECK_EQ(Orthant<3>::bne().collapseOnAxis(1), Orthant<2>::se());
  CHECK_EQ(Orthant<3>::bne().collapseOnAxis(2), Orthant<2>::ne());

  CHECK_EQ(Orthant<3>::tsw().collapseOnAxis(0), Orthant<2>::nw());
  CHECK_EQ(Orthant<3>::tsw().collapseOnAxis(1), Orthant<2>::nw());
  CHECK_EQ(Orthant<3>::tsw().collapseOnAxis(2), Orthant<2>::sw());

  CHECK_EQ(Orthant<3>::tse().collapseOnAxis(0), Orthant<2>::nw());
  CHECK_EQ(Orthant<3>::tse().collapseOnAxis(1), Orthant<2>::ne());
  CHECK_EQ(Orthant<3>::tse().collapseOnAxis(2), Orthant<2>::se());

  CHECK_EQ(Orthant<3>::tnw().collapseOnAxis(0), Orthant<2>::ne());
  CHECK_EQ(Orthant<3>::tnw().collapseOnAxis(1), Orthant<2>::nw());
  CHECK_EQ(Orthant<3>::tnw().collapseOnAxis(2), Orthant<2>::nw());

  CHECK_EQ(Orthant<3>::tne().collapseOnAxis(0), Orthant<2>::ne());
  CHECK_EQ(Orthant<3>::tne().collapseOnAxis(1), Orthant<2>::ne());
  CHECK_EQ(Orthant<3>::tne().collapseOnAxis(2), Orthant<2>::ne());
}
TEST_CASE("Orthant<0> ==")
{
  CHECK_EQ(Orthant<0>::null(), Orthant<0>::null());
}
TEST_CASE("Orthant<1> ==")
{
  CHECK_EQ(Orthant<1>::lower(), Orthant<1>::lower());
  CHECK_UNARY_FALSE(Orthant<1>::lower() == Orthant<1>::upper());
  CHECK_UNARY_FALSE(Orthant<1>::lower() == Orthant<1>::null());

  CHECK_UNARY_FALSE(Orthant<1>::upper() == Orthant<1>::lower());
  CHECK_EQ(Orthant<1>::upper(), Orthant<1>::upper());
  CHECK_UNARY_FALSE(Orthant<1>::upper() == Orthant<1>::null());

  CHECK_UNARY_FALSE(Orthant<1>::null() == Orthant<1>::lower());
  CHECK_UNARY_FALSE(Orthant<1>::null() == Orthant<1>::upper());
  CHECK_EQ(Orthant<1>::null(), Orthant<1>::null());
}
TEST_CASE("Orthant<2> ==")
{
  CHECK_EQ(Orthant<2>::sw(), Orthant<2>::sw());
  CHECK_UNARY_FALSE(Orthant<2>::sw() == Orthant<2>::se());
  CHECK_UNARY_FALSE(Orthant<2>::sw() == Orthant<2>::nw());
  CHECK_UNARY_FALSE(Orthant<2>::sw() == Orthant<2>::ne());
  CHECK_UNARY_FALSE(Orthant<2>::sw() == Orthant<2>::null());

  CHECK_UNARY_FALSE(Orthant<2>::se() == Orthant<2>::sw());
  CHECK_EQ(Orthant<2>::se(), Orthant<2>::se());
  CHECK_UNARY_FALSE(Orthant<2>::se() == Orthant<2>::nw());
  CHECK_UNARY_FALSE(Orthant<2>::se() == Orthant<2>::ne());
  CHECK_UNARY_FALSE(Orthant<2>::se() == Orthant<2>::null());

  CHECK_UNARY_FALSE(Orthant<2>::nw() == Orthant<2>::sw());
  CHECK_UNARY_FALSE(Orthant<2>::nw() == Orthant<2>::se());
  CHECK_EQ(Orthant<2>::nw(), Orthant<2>::nw());
  CHECK_UNARY_FALSE(Orthant<2>::nw() == Orthant<2>::ne());
  CHECK_UNARY_FALSE(Orthant<2>::nw() == Orthant<2>::null());

  CHECK_UNARY_FALSE(Orthant<2>::ne() == Orthant<2>::sw());
  CHECK_UNARY_FALSE(Orthant<2>::ne() == Orthant<2>::se());
  CHECK_UNARY_FALSE(Orthant<2>::ne() == Orthant<2>::nw());
  CHECK_EQ(Orthant<2>::ne(), Orthant<2>::ne());
  CHECK_UNARY_FALSE(Orthant<2>::ne() == Orthant<2>::null());

  CHECK_UNARY_FALSE(Orthant<2>::null() == Orthant<2>::sw());
  CHECK_UNARY_FALSE(Orthant<2>::null() == Orthant<2>::se());
  CHECK_UNARY_FALSE(Orthant<2>::null() == Orthant<2>::nw());
  CHECK_UNARY_FALSE(Orthant<2>::null() == Orthant<2>::ne());
  CHECK_EQ(Orthant<2>::null(), Orthant<2>::null());
}
TEST_CASE("Orthant<3> ==")
{
  CHECK_EQ(Orthant<3>::bsw(), Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::bsw() == Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::bsw() == Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::bsw() == Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::bsw() == Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::bsw() == Orthant<3>::tse());
  CHECK_UNARY_FALSE(Orthant<3>::bsw() == Orthant<3>::tnw());
  CHECK_UNARY_FALSE(Orthant<3>::bsw() == Orthant<3>::tne());
  CHECK_UNARY_FALSE(Orthant<3>::bsw() == Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::bse() == Orthant<3>::bsw());
  CHECK_EQ(Orthant<3>::bse(), Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::bse() == Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::bse() == Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::bse() == Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::bse() == Orthant<3>::tse());
  CHECK_UNARY_FALSE(Orthant<3>::bse() == Orthant<3>::tnw());
  CHECK_UNARY_FALSE(Orthant<3>::bse() == Orthant<3>::tne());
  CHECK_UNARY_FALSE(Orthant<3>::bse() == Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::bnw() == Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::bnw() == Orthant<3>::bse());
  CHECK_EQ(Orthant<3>::bnw(), Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::bnw() == Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::bnw() == Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::bnw() == Orthant<3>::tse());
  CHECK_UNARY_FALSE(Orthant<3>::bnw() == Orthant<3>::tnw());
  CHECK_UNARY_FALSE(Orthant<3>::bnw() == Orthant<3>::tne());
  CHECK_UNARY_FALSE(Orthant<3>::bnw() == Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::bne() == Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::bne() == Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::bne() == Orthant<3>::bnw());
  CHECK_EQ(Orthant<3>::bne(), Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::bne() == Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::bne() == Orthant<3>::tse());
  CHECK_UNARY_FALSE(Orthant<3>::bne() == Orthant<3>::tnw());
  CHECK_UNARY_FALSE(Orthant<3>::bne() == Orthant<3>::tne());
  CHECK_UNARY_FALSE(Orthant<3>::bne() == Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::tsw() == Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::tsw() == Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::tsw() == Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::tsw() == Orthant<3>::bne());
  CHECK_EQ(Orthant<3>::tsw(), Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::tsw() == Orthant<3>::tse());
  CHECK_UNARY_FALSE(Orthant<3>::tsw() == Orthant<3>::tnw());
  CHECK_UNARY_FALSE(Orthant<3>::tsw() == Orthant<3>::tne());
  CHECK_UNARY_FALSE(Orthant<3>::tsw() == Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::tse() == Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::tse() == Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::tse() == Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::tse() == Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::tse() == Orthant<3>::tsw());
  CHECK_EQ(Orthant<3>::tse(), Orthant<3>::tse());
  CHECK_UNARY_FALSE(Orthant<3>::tse() == Orthant<3>::tnw());
  CHECK_UNARY_FALSE(Orthant<3>::tse() == Orthant<3>::tne());
  CHECK_UNARY_FALSE(Orthant<3>::tse() == Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::tnw() == Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() == Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() == Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() == Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() == Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() == Orthant<3>::tse());
  CHECK_EQ(Orthant<3>::tnw(), Orthant<3>::tnw());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() == Orthant<3>::tne());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() == Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::tne() == Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::tne() == Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::tne() == Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::tne() == Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::tne() == Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::tne() == Orthant<3>::tse());
  CHECK_UNARY_FALSE(Orthant<3>::tne() == Orthant<3>::tnw());
  CHECK_EQ(Orthant<3>::tne(), Orthant<3>::tne());
  CHECK_UNARY_FALSE(Orthant<3>::tne() == Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::null() == Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::null() == Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::null() == Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::null() == Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::null() == Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::null() == Orthant<3>::tse());
  CHECK_UNARY_FALSE(Orthant<3>::null() == Orthant<3>::tnw());
  CHECK_UNARY_FALSE(Orthant<3>::null() == Orthant<3>::tne());
  CHECK_EQ(Orthant<3>::null(), Orthant<3>::null());
}
TEST_CASE("Orthant<1> !=")
{
  CHECK_UNARY_FALSE(Orthant<1>::lower() != Orthant<1>::lower());
  CHECK_NE(Orthant<1>::lower(), Orthant<1>::upper());
  CHECK_NE(Orthant<1>::lower(), Orthant<1>::null());

  CHECK_NE(Orthant<1>::upper(), Orthant<1>::lower());
  CHECK_UNARY_FALSE(Orthant<1>::upper() != Orthant<1>::upper());
  CHECK_NE(Orthant<1>::upper(), Orthant<1>::null());

  CHECK_NE(Orthant<1>::null(), Orthant<1>::lower());
  CHECK_NE(Orthant<1>::null(), Orthant<1>::upper());
  CHECK_UNARY_FALSE(Orthant<1>::null() != Orthant<1>::null());
}
TEST_CASE("Orthant<2> !=")
{
  CHECK_UNARY_FALSE(Orthant<2>::sw() != Orthant<2>::sw());
  CHECK_NE(Orthant<2>::sw(), Orthant<2>::se());
  CHECK_NE(Orthant<2>::sw(), Orthant<2>::nw());
  CHECK_NE(Orthant<2>::sw(), Orthant<2>::ne());
  CHECK_NE(Orthant<2>::sw(), Orthant<2>::null());

  CHECK_NE(Orthant<2>::se(), Orthant<2>::sw());
  CHECK_UNARY_FALSE(Orthant<2>::se() != Orthant<2>::se());
  CHECK_NE(Orthant<2>::se(), Orthant<2>::nw());
  CHECK_NE(Orthant<2>::se(), Orthant<2>::ne());
  CHECK_NE(Orthant<2>::se(), Orthant<2>::null());

  CHECK_NE(Orthant<2>::nw(), Orthant<2>::sw());
  CHECK_NE(Orthant<2>::nw(), Orthant<2>::se());
  CHECK_UNARY_FALSE(Orthant<2>::nw() != Orthant<2>::nw());
  CHECK_NE(Orthant<2>::nw(), Orthant<2>::ne());
  CHECK_NE(Orthant<2>::nw(), Orthant<2>::null());

  CHECK_NE(Orthant<2>::ne(), Orthant<2>::sw());
  CHECK_NE(Orthant<2>::ne(), Orthant<2>::se());
  CHECK_NE(Orthant<2>::ne(), Orthant<2>::nw());
  CHECK_UNARY_FALSE(Orthant<2>::ne() != Orthant<2>::ne());
  CHECK_NE(Orthant<2>::ne(), Orthant<2>::null());

  CHECK_NE(Orthant<2>::null(), Orthant<2>::sw());
  CHECK_NE(Orthant<2>::null(), Orthant<2>::se());
  CHECK_NE(Orthant<2>::null(), Orthant<2>::nw());
  CHECK_NE(Orthant<2>::null(), Orthant<2>::ne());
  CHECK_UNARY_FALSE(Orthant<2>::null() != Orthant<2>::null());
}
TEST_CASE("Orthant<3> !=")
{
  CHECK_UNARY_FALSE(Orthant<3>::bsw() != Orthant<3>::bsw());
  CHECK_NE(Orthant<3>::bsw(), Orthant<3>::bse());
  CHECK_NE(Orthant<3>::bsw(), Orthant<3>::bnw());
  CHECK_NE(Orthant<3>::bsw(), Orthant<3>::bne());
  CHECK_NE(Orthant<3>::bsw(), Orthant<3>::tsw());
  CHECK_NE(Orthant<3>::bsw(), Orthant<3>::tse());
  CHECK_NE(Orthant<3>::bsw(), Orthant<3>::tnw());
  CHECK_NE(Orthant<3>::bsw(), Orthant<3>::tne());
  CHECK_NE(Orthant<3>::bsw(), Orthant<3>::null());

  CHECK_NE(Orthant<3>::bse(), Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::bse() != Orthant<3>::bse());
  CHECK_NE(Orthant<3>::bse(), Orthant<3>::bnw());
  CHECK_NE(Orthant<3>::bse(), Orthant<3>::bne());
  CHECK_NE(Orthant<3>::bse(), Orthant<3>::tsw());
  CHECK_NE(Orthant<3>::bse(), Orthant<3>::tse());
  CHECK_NE(Orthant<3>::bse(), Orthant<3>::tnw());
  CHECK_NE(Orthant<3>::bse(), Orthant<3>::tne());
  CHECK_NE(Orthant<3>::bse(), Orthant<3>::null());

  CHECK_NE(Orthant<3>::bnw(), Orthant<3>::bsw());
  CHECK_NE(Orthant<3>::bnw(), Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::bnw() != Orthant<3>::bnw());
  CHECK_NE(Orthant<3>::bnw(), Orthant<3>::bne());
  CHECK_NE(Orthant<3>::bnw(), Orthant<3>::tsw());
  CHECK_NE(Orthant<3>::bnw(), Orthant<3>::tse());
  CHECK_NE(Orthant<3>::bnw(), Orthant<3>::tnw());
  CHECK_NE(Orthant<3>::bnw(), Orthant<3>::tne());
  CHECK_NE(Orthant<3>::bnw(), Orthant<3>::null());

  CHECK_NE(Orthant<3>::bne(), Orthant<3>::bsw());
  CHECK_NE(Orthant<3>::bne(), Orthant<3>::bse());
  CHECK_NE(Orthant<3>::bne(), Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::bne() != Orthant<3>::bne());
  CHECK_NE(Orthant<3>::bne(), Orthant<3>::tsw());
  CHECK_NE(Orthant<3>::bne(), Orthant<3>::tse());
  CHECK_NE(Orthant<3>::bne(), Orthant<3>::tnw());
  CHECK_NE(Orthant<3>::bne(), Orthant<3>::tne());
  CHECK_NE(Orthant<3>::bne(), Orthant<3>::null());

  CHECK_NE(Orthant<3>::tsw(), Orthant<3>::bsw());
  CHECK_NE(Orthant<3>::tsw(), Orthant<3>::bse());
  CHECK_NE(Orthant<3>::tsw(), Orthant<3>::bnw());
  CHECK_NE(Orthant<3>::tsw(), Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::tsw() != Orthant<3>::tsw());
  CHECK_NE(Orthant<3>::tsw(), Orthant<3>::tse());
  CHECK_NE(Orthant<3>::tsw(), Orthant<3>::tnw());
  CHECK_NE(Orthant<3>::tsw(), Orthant<3>::tne());
  CHECK_NE(Orthant<3>::tsw(), Orthant<3>::null());

  CHECK_NE(Orthant<3>::tse(), Orthant<3>::bsw());
  CHECK_NE(Orthant<3>::tse(), Orthant<3>::bse());
  CHECK_NE(Orthant<3>::tse(), Orthant<3>::bnw());
  CHECK_NE(Orthant<3>::tse(), Orthant<3>::bne());
  CHECK_NE(Orthant<3>::tse(), Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::tse() != Orthant<3>::tse());
  CHECK_NE(Orthant<3>::tse(), Orthant<3>::tnw());
  CHECK_NE(Orthant<3>::tse(), Orthant<3>::tne());
  CHECK_NE(Orthant<3>::tse(), Orthant<3>::null());

  CHECK_NE(Orthant<3>::tnw(), Orthant<3>::bsw());
  CHECK_NE(Orthant<3>::tnw(), Orthant<3>::bse());
  CHECK_NE(Orthant<3>::tnw(), Orthant<3>::bnw());
  CHECK_NE(Orthant<3>::tnw(), Orthant<3>::bne());
  CHECK_NE(Orthant<3>::tnw(), Orthant<3>::tsw());
  CHECK_NE(Orthant<3>::tnw(), Orthant<3>::tse());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() != Orthant<3>::tnw());
  CHECK_NE(Orthant<3>::tnw(), Orthant<3>::tne());
  CHECK_NE(Orthant<3>::tnw(), Orthant<3>::null());

  CHECK_NE(Orthant<3>::tne(), Orthant<3>::bsw());
  CHECK_NE(Orthant<3>::tne(), Orthant<3>::bse());
  CHECK_NE(Orthant<3>::tne(), Orthant<3>::bnw());
  CHECK_NE(Orthant<3>::tne(), Orthant<3>::bne());
  CHECK_NE(Orthant<3>::tne(), Orthant<3>::tsw());
  CHECK_NE(Orthant<3>::tne(), Orthant<3>::tse());
  CHECK_NE(Orthant<3>::tne(), Orthant<3>::tnw());
  CHECK_UNARY_FALSE(Orthant<3>::tne() != Orthant<3>::tne());
  CHECK_NE(Orthant<3>::tne(), Orthant<3>::null());

  CHECK_NE(Orthant<3>::null(), Orthant<3>::bsw());
  CHECK_NE(Orthant<3>::null(), Orthant<3>::bse());
  CHECK_NE(Orthant<3>::null(), Orthant<3>::bnw());
  CHECK_NE(Orthant<3>::null(), Orthant<3>::bne());
  CHECK_NE(Orthant<3>::null(), Orthant<3>::tsw());
  CHECK_NE(Orthant<3>::null(), Orthant<3>::tse());
  CHECK_NE(Orthant<3>::null(), Orthant<3>::tnw());
  CHECK_NE(Orthant<3>::null(), Orthant<3>::tne());
  CHECK_UNARY_FALSE(Orthant<3>::null() != Orthant<3>::null());
}
TEST_CASE("Orthant<1> <")
{
  CHECK_UNARY_FALSE(Orthant<1>::lower() < Orthant<1>::lower());
  CHECK_LT(Orthant<1>::lower(), Orthant<1>::upper());
  CHECK_LT(Orthant<1>::lower(), Orthant<1>::null());

  CHECK_UNARY_FALSE(Orthant<1>::upper() < Orthant<1>::lower());
  CHECK_UNARY_FALSE(Orthant<1>::upper() < Orthant<1>::upper());
  CHECK_LT(Orthant<1>::upper(), Orthant<1>::null());

  CHECK_UNARY_FALSE(Orthant<1>::null() < Orthant<1>::lower());
  CHECK_UNARY_FALSE(Orthant<1>::null() < Orthant<1>::upper());
  CHECK_UNARY_FALSE(Orthant<1>::null() < Orthant<1>::null());
}
TEST_CASE("Orthant<2> <")
{
  CHECK_UNARY_FALSE(Orthant<2>::sw() < Orthant<2>::sw());
  CHECK_LT(Orthant<2>::sw(), Orthant<2>::se());
  CHECK_LT(Orthant<2>::sw(), Orthant<2>::nw());
  CHECK_LT(Orthant<2>::sw(), Orthant<2>::ne());
  CHECK_LT(Orthant<2>::sw(), Orthant<2>::null());

  CHECK_UNARY_FALSE(Orthant<2>::se() < Orthant<2>::sw());
  CHECK_UNARY_FALSE(Orthant<2>::se() < Orthant<2>::se());
  CHECK_LT(Orthant<2>::se(), Orthant<2>::nw());
  CHECK_LT(Orthant<2>::se(), Orthant<2>::ne());
  CHECK_LT(Orthant<2>::se(), Orthant<2>::null());

  CHECK_UNARY_FALSE(Orthant<2>::nw() < Orthant<2>::sw());
  CHECK_UNARY_FALSE(Orthant<2>::nw() < Orthant<2>::se());
  CHECK_UNARY_FALSE(Orthant<2>::nw() < Orthant<2>::nw());
  CHECK_LT(Orthant<2>::nw(), Orthant<2>::ne());
  CHECK_LT(Orthant<2>::nw(), Orthant<2>::null());

  CHECK_UNARY_FALSE(Orthant<2>::ne() < Orthant<2>::sw());
  CHECK_UNARY_FALSE(Orthant<2>::ne() < Orthant<2>::se());
  CHECK_UNARY_FALSE(Orthant<2>::ne() < Orthant<2>::nw());
  CHECK_UNARY_FALSE(Orthant<2>::ne() < Orthant<2>::ne());
  CHECK_LT(Orthant<2>::ne(), Orthant<2>::null());

  CHECK_UNARY_FALSE(Orthant<2>::null() < Orthant<2>::sw());
  CHECK_UNARY_FALSE(Orthant<2>::null() < Orthant<2>::se());
  CHECK_UNARY_FALSE(Orthant<2>::null() < Orthant<2>::nw());
  CHECK_UNARY_FALSE(Orthant<2>::null() < Orthant<2>::ne());
  CHECK_UNARY_FALSE(Orthant<2>::null() < Orthant<2>::null());
}
TEST_CASE("Orthant<3> <")
{
  CHECK_UNARY_FALSE(Orthant<3>::bsw() < Orthant<3>::bsw());
  CHECK_LT(Orthant<3>::bsw(), Orthant<3>::bse());
  CHECK_LT(Orthant<3>::bsw(), Orthant<3>::bnw());
  CHECK_LT(Orthant<3>::bsw(), Orthant<3>::bne());
  CHECK_LT(Orthant<3>::bsw(), Orthant<3>::tsw());
  CHECK_LT(Orthant<3>::bsw(), Orthant<3>::tse());
  CHECK_LT(Orthant<3>::bsw(), Orthant<3>::tnw());
  CHECK_LT(Orthant<3>::bsw(), Orthant<3>::tne());
  CHECK_LT(Orthant<3>::bsw(), Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::bse() < Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::bse() < Orthant<3>::bse());
  CHECK_LT(Orthant<3>::bse(), Orthant<3>::bnw());
  CHECK_LT(Orthant<3>::bse(), Orthant<3>::bne());
  CHECK_LT(Orthant<3>::bse(), Orthant<3>::tsw());
  CHECK_LT(Orthant<3>::bse(), Orthant<3>::tse());
  CHECK_LT(Orthant<3>::bse(), Orthant<3>::tnw());
  CHECK_LT(Orthant<3>::bse(), Orthant<3>::tne());
  CHECK_LT(Orthant<3>::bse(), Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::bnw() < Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::bnw() < Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::bnw() < Orthant<3>::bnw());
  CHECK_LT(Orthant<3>::bnw(), Orthant<3>::bne());
  CHECK_LT(Orthant<3>::bnw(), Orthant<3>::tsw());
  CHECK_LT(Orthant<3>::bnw(), Orthant<3>::tse());
  CHECK_LT(Orthant<3>::bnw(), Orthant<3>::tnw());
  CHECK_LT(Orthant<3>::bnw(), Orthant<3>::tne());
  CHECK_LT(Orthant<3>::bnw(), Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::bne() < Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::bne() < Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::bne() < Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::bne() < Orthant<3>::bne());
  CHECK_LT(Orthant<3>::bne(), Orthant<3>::tsw());
  CHECK_LT(Orthant<3>::bne(), Orthant<3>::tse());
  CHECK_LT(Orthant<3>::bne(), Orthant<3>::tnw());
  CHECK_LT(Orthant<3>::bne(), Orthant<3>::tne());
  CHECK_LT(Orthant<3>::bne(), Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::tsw() < Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::tsw() < Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::tsw() < Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::tsw() < Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::tsw() < Orthant<3>::tsw());
  CHECK_LT(Orthant<3>::tsw(), Orthant<3>::tse());
  CHECK_LT(Orthant<3>::tsw(), Orthant<3>::tnw());
  CHECK_LT(Orthant<3>::tsw(), Orthant<3>::tne());
  CHECK_LT(Orthant<3>::tsw(), Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::tse() < Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::tse() < Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::tse() < Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::tse() < Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::tse() < Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::tse() < Orthant<3>::tse());
  CHECK_LT(Orthant<3>::tse(), Orthant<3>::tnw());
  CHECK_LT(Orthant<3>::tse(), Orthant<3>::tne());
  CHECK_LT(Orthant<3>::tse(), Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::tnw() < Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() < Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() < Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() < Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() < Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() < Orthant<3>::tse());
  CHECK_UNARY_FALSE(Orthant<3>::tnw() < Orthant<3>::tnw());
  CHECK_LT(Orthant<3>::tnw(), Orthant<3>::tne());
  CHECK_LT(Orthant<3>::tnw(), Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::tne() < Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::tne() < Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::tne() < Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::tne() < Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::tne() < Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::tne() < Orthant<3>::tse());
  CHECK_UNARY_FALSE(Orthant<3>::tne() < Orthant<3>::tnw());
  CHECK_UNARY_FALSE(Orthant<3>::tne() < Orthant<3>::tne());
  CHECK_LT(Orthant<3>::tne(), Orthant<3>::null());

  CHECK_UNARY_FALSE(Orthant<3>::null() < Orthant<3>::bsw());
  CHECK_UNARY_FALSE(Orthant<3>::null() < Orthant<3>::bse());
  CHECK_UNARY_FALSE(Orthant<3>::null() < Orthant<3>::bnw());
  CHECK_UNARY_FALSE(Orthant<3>::null() < Orthant<3>::bne());
  CHECK_UNARY_FALSE(Orthant<3>::null() < Orthant<3>::tsw());
  CHECK_UNARY_FALSE(Orthant<3>::null() < Orthant<3>::tse());
  CHECK_UNARY_FALSE(Orthant<3>::null() < Orthant<3>::tnw());
  CHECK_UNARY_FALSE(Orthant<3>::null() < Orthant<3>::tne());
  CHECK_UNARY_FALSE(Orthant<3>::null() < Orthant<3>::null());
}
TEST_CASE("Test ostream for Orthant<0>")
{
  stringstream ss;
  ss << Orthant<0>::null();
  CHECK_EQ(ss.str(), "Orthant<0>::null()");
  ss.str("");
  ss << Orthant<0>(13);
  CHECK_EQ(ss.str(), "Orthant<0> invalid value: 13");
}
TEST_CASE("Test ostream for Orthant<1>")
{
  stringstream ss;
  ss << Orthant<1>::lower();
  CHECK_EQ(ss.str(), "Orthant<1>::lower()");
  ss.str("");
  ss << Orthant<1>::upper();
  CHECK_EQ(ss.str(), "Orthant<1>::upper()");
  ss.str("");
  ss << Orthant<1>::null();
  CHECK_EQ(ss.str(), "Orthant<1>::null()");
  ss.str("");
  ss << Orthant<1>(13);
  CHECK_EQ(ss.str(), "Orthant<1> invalid value: 13");
}
TEST_CASE("Test ostream for Orthant<2>")
{
  stringstream ss;
  ss << Orthant<2>::sw();
  CHECK_EQ(ss.str(), "Orthant<2>::sw()");
  ss.str("");
  ss << Orthant<2>::se();
  CHECK_EQ(ss.str(), "Orthant<2>::se()");
  ss.str("");
  ss << Orthant<2>::nw();
  CHECK_EQ(ss.str(), "Orthant<2>::nw()");
  ss.str("");
  ss << Orthant<2>::ne();
  CHECK_EQ(ss.str(), "Orthant<2>::ne()");
  ss.str("");
  ss << Orthant<2>::null();
  CHECK_EQ(ss.str(), "Orthant<2>::null()");
  ss.str("");
  ss << Orthant<2>(13);
  CHECK_EQ(ss.str(), "Orthant<2> invalid value: 13");
}
TEST_CASE("Test ostream for Orthant<3>")
{
  stringstream ss;
  ss << Orthant<3>::bsw();
  CHECK_EQ(ss.str(), "Orthant<3>::bsw()");
  ss.str("");
  ss << Orthant<3>::bse();
  CHECK_EQ(ss.str(), "Orthant<3>::bse()");
  ss.str("");
  ss << Orthant<3>::bnw();
  CHECK_EQ(ss.str(), "Orthant<3>::bnw()");
  ss.str("");
  ss << Orthant<3>::bne();
  CHECK_EQ(ss.str(), "Orthant<3>::bne()");
  ss.str("");
  ss << Orthant<3>::tsw();
  CHECK_EQ(ss.str(), "Orthant<3>::tsw()");
  ss.str("");
  ss << Orthant<3>::tse();
  CHECK_EQ(ss.str(), "Orthant<3>::tse()");
  ss.str("");
  ss << Orthant<3>::tnw();
  CHECK_EQ(ss.str(), "Orthant<3>::tnw()");
  ss.str("");
  ss << Orthant<3>::tne();
  CHECK_EQ(ss.str(), "Orthant<3>::tne()");
  ss.str("");
  ss << Orthant<3>::null();
  CHECK_EQ(ss.str(), "Orthant<3>::null()");
  ss.str("");
  ss << Orthant<3>(13);
  CHECK_EQ(ss.str(), "Orthant<3> invalid value: 13");
}
TEST_CASE("Test iterator for Orthant<0>")
{
  auto iter = Orthant<0>::getValues().begin();
  CHECK_EQ(*iter, Orthant<0>(0));
  ++iter;
  CHECK_EQ(*iter, Orthant<0>::null());
  CHECK_EQ(iter, Orthant<0>::getValues().end());
}
TEST_CASE("Test iterator for Orthant<1>")
{
  auto iter = Orthant<1>::getValues().begin();
  CHECK_EQ(iter, Orthant<1>::getValues().begin());
  CHECK_NE(iter, Orthant<1>::getValues().end());
  CHECK_EQ(*iter, Orthant<1>::lower());
  ++iter;
  CHECK_EQ(iter->getIndex(), 1);
  CHECK_EQ(*iter, Orthant<1>::upper());
  ++iter;
  CHECK_EQ(*iter, Orthant<1>::null());
  CHECK_EQ(iter, Orthant<1>::getValues().end());
}
TEST_CASE("Test iterator for Orthant<2>")
{
  auto iter = Orthant<2>::getValues().begin();
  CHECK_EQ(iter, Orthant<2>::getValues().begin());
  CHECK_NE(iter, Orthant<2>::getValues().end());
  CHECK_EQ(*iter, Orthant<2>::sw());
  ++iter;
  CHECK_EQ(iter->getIndex(), 1);
  CHECK_EQ(*iter, Orthant<2>::se());
  ++iter;
  CHECK_EQ(*iter, Orthant<2>::nw());
  ++iter;
  CHECK_EQ(*iter, Orthant<2>::ne());
  ++iter;
  CHECK_EQ(*iter, Orthant<2>::null());
  CHECK_EQ(iter, Orthant<2>::getValues().end());
}
TEST_CASE("Test iterator for Orthant<3>")
{
  auto iter = Orthant<3>::getValues().begin();
  CHECK_EQ(iter, Orthant<3>::getValues().begin());
  CHECK_NE(iter, Orthant<3>::getValues().end());
  CHECK_EQ(*iter, Orthant<3>::bsw());
  ++iter;
  CHECK_EQ(iter->getIndex(), 1);
  CHECK_EQ(*iter, Orthant<3>::bse());
  ++iter;
  CHECK_EQ(*iter, Orthant<3>::bnw());
  ++iter;
  CHECK_EQ(*iter, Orthant<3>::bne());
  ++iter;
  CHECK_EQ(*iter, Orthant<3>::tsw());
  ++iter;
  CHECK_EQ(*iter, Orthant<3>::tse());
  ++iter;
  CHECK_EQ(*iter, Orthant<3>::tnw());
  ++iter;
  CHECK_EQ(*iter, Orthant<3>::tne());
  ++iter;
  CHECK_EQ(*iter, Orthant<3>::null());
  CHECK_EQ(iter, Orthant<3>::getValues().end());
}
