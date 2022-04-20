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
#include <ThunderEgg/Face.h>
#include <sstream>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::tpl;
TEST_CASE("Test num_sides for Side<1>")
{
  size_t sides = Side<1>::number_of;
  CHECK_EQ(sides, 2);
}
TEST_CASE("Test num_sides for Side<2>")
{
  size_t sides = Side<2>::number_of;
  CHECK_EQ(sides, 4);
}
TEST_CASE("Test num_sides for Side<3>")
{
  size_t sides = Side<3>::number_of;
  CHECK_EQ(sides, 6);
}
TEST_CASE("Test dimensionality for Side<1>")
{
  size_t dimensionality = Side<1>::dimensionality;
  CHECK_EQ(dimensionality, 0);
}
TEST_CASE("Test dimensionality for Side<2>")
{
  size_t dimensionality = Side<2>::dimensionality;
  CHECK_EQ(dimensionality, 1);
}
TEST_CASE("Test dimensionality for Side<3>")
{
  size_t dimensionality = Side<3>::dimensionality;
  CHECK_EQ(dimensionality, 2);
}
TEST_CASE("Test that named constructors for Side<1> give expected index values")
{
  CHECK_EQ(Side<1>::west().getIndex(), 0);
  CHECK_EQ(Side<1>::east().getIndex(), 1);
  CHECK_EQ(Side<1>::null().getIndex(), 2);
}
TEST_CASE("Test that named constructors for Side<2> give expected index values")
{
  CHECK_EQ(Side<2>::west().getIndex(), 0);
  CHECK_EQ(Side<2>::east().getIndex(), 1);
  CHECK_EQ(Side<2>::south().getIndex(), 2);
  CHECK_EQ(Side<2>::north().getIndex(), 3);
  CHECK_EQ(Side<2>::null().getIndex(), 4);
}
TEST_CASE("Test that named constructors for Side<3> give expected index values")
{
  CHECK_EQ(Side<3>::west().getIndex(), 0);
  CHECK_EQ(Side<3>::east().getIndex(), 1);
  CHECK_EQ(Side<3>::south().getIndex(), 2);
  CHECK_EQ(Side<3>::north().getIndex(), 3);
  CHECK_EQ(Side<3>::bottom().getIndex(), 4);
  CHECK_EQ(Side<3>::top().getIndex(), 5);
  CHECK_EQ(Side<3>::null().getIndex(), 6);
}
TEST_CASE("Test that named constructors for Side<1> give expected axis index values")
{
  CHECK_EQ(Side<1>::west().getAxisIndex(), 0);
  CHECK_EQ(Side<1>::east().getAxisIndex(), 0);
  CHECK_EQ(Side<1>::null().getAxisIndex(), 1);
}
TEST_CASE("Test that named constructors for Side<2> give expected axis index values")
{
  CHECK_EQ(Side<2>::west().getAxisIndex(), 0);
  CHECK_EQ(Side<2>::east().getAxisIndex(), 0);
  CHECK_EQ(Side<2>::south().getAxisIndex(), 1);
  CHECK_EQ(Side<2>::north().getAxisIndex(), 1);
  CHECK_EQ(Side<2>::null().getAxisIndex(), 2);
}
TEST_CASE("Test that named constructors for Side<3> give expected axis index values")
{
  CHECK_EQ(Side<3>::west().getAxisIndex(), 0);
  CHECK_EQ(Side<3>::east().getAxisIndex(), 0);
  CHECK_EQ(Side<3>::south().getAxisIndex(), 1);
  CHECK_EQ(Side<3>::north().getAxisIndex(), 1);
  CHECK_EQ(Side<3>::bottom().getAxisIndex(), 2);
  CHECK_EQ(Side<3>::top().getAxisIndex(), 2);
  CHECK_EQ(Side<3>::null().getAxisIndex(), 3);
}
TEST_CASE("Test opposite for Side<1>")
{
  CHECK_EQ(Side<1>::west().opposite(), Side<1>::east());
  CHECK_EQ(Side<1>::east().opposite(), Side<1>::west());
}
TEST_CASE("Test opposite for Side<2>")
{
  CHECK_EQ(Side<2>::west().opposite(), Side<2>::east());
  CHECK_EQ(Side<2>::east().opposite(), Side<2>::west());
  CHECK_EQ(Side<2>::south().opposite(), Side<2>::north());
  CHECK_EQ(Side<2>::north().opposite(), Side<2>::south());
}
TEST_CASE("Test opposite for Side<3>")
{
  CHECK_EQ(Side<3>::west().opposite(), Side<3>::east());
  CHECK_EQ(Side<3>::east().opposite(), Side<3>::west());
  CHECK_EQ(Side<3>::south().opposite(), Side<3>::north());
  CHECK_EQ(Side<3>::north().opposite(), Side<3>::south());
  CHECK_EQ(Side<3>::bottom().opposite(), Side<3>::top());
  CHECK_EQ(Side<3>::top().opposite(), Side<3>::bottom());
}
TEST_CASE("Test getSides for Side<1>")
{
  CHECK_EQ(Side<1>::west().getSides()[0], Side<1>::west());
  CHECK_EQ(Side<1>::east().getSides()[0], Side<1>::east());
}
TEST_CASE("Test getSides for Side<2>")
{
  CHECK_EQ(Side<2>::west().getSides()[0], Side<2>::west());
  CHECK_EQ(Side<2>::east().getSides()[0], Side<2>::east());
  CHECK_EQ(Side<2>::south().getSides()[0], Side<2>::south());
  CHECK_EQ(Side<2>::north().getSides()[0], Side<2>::north());
}
TEST_CASE("Test getSides for Side<3>")
{
  CHECK_EQ(Side<3>::west().getSides()[0], Side<3>::west());
  CHECK_EQ(Side<3>::east().getSides()[0], Side<3>::east());
  CHECK_EQ(Side<3>::south().getSides()[0], Side<3>::south());
  CHECK_EQ(Side<3>::north().getSides()[0], Side<3>::north());
  CHECK_EQ(Side<3>::bottom().getSides()[0], Side<3>::bottom());
  CHECK_EQ(Side<3>::top().getSides()[0], Side<3>::top());
}
TEST_CASE("Test isLowerOnAxis for Side<1>")
{
  CHECK_UNARY(Side<1>::west().isLowerOnAxis());
  CHECK_UNARY_FALSE(Side<1>::east().isLowerOnAxis());
}
TEST_CASE("Test isLowerOnAxis for Side<2>")
{
  CHECK_UNARY(Side<2>::west().isLowerOnAxis());
  CHECK_UNARY_FALSE(Side<2>::east().isLowerOnAxis());
  CHECK_UNARY(Side<2>::south().isLowerOnAxis());
  CHECK_UNARY_FALSE(Side<2>::north().isLowerOnAxis());
}
TEST_CASE("Test isLowerOnAxis for Side<3>")
{
  CHECK_UNARY(Side<3>::west().isLowerOnAxis());
  CHECK_UNARY_FALSE(Side<3>::east().isLowerOnAxis());
  CHECK_UNARY(Side<3>::south().isLowerOnAxis());
  CHECK_UNARY_FALSE(Side<3>::north().isLowerOnAxis());
  CHECK_UNARY(Side<3>::bottom().isLowerOnAxis());
  CHECK_UNARY_FALSE(Side<3>::top().isLowerOnAxis());
}
TEST_CASE("Test isHigherOnAxis for Side<1>")
{
  CHECK_UNARY_FALSE(Side<1>::west().isHigherOnAxis());
  CHECK_UNARY(Side<1>::east().isHigherOnAxis());
}
TEST_CASE("Test isHigherOnAxis for Side<2>")
{
  CHECK_UNARY_FALSE(Side<2>::west().isHigherOnAxis());
  CHECK_UNARY(Side<2>::east().isHigherOnAxis());
  CHECK_UNARY_FALSE(Side<2>::south().isHigherOnAxis());
  CHECK_UNARY(Side<2>::north().isHigherOnAxis());
}
TEST_CASE("Test isHigherOnAxis for Side<3>")
{
  CHECK_UNARY_FALSE(Side<3>::west().isHigherOnAxis());
  CHECK_UNARY(Side<3>::east().isHigherOnAxis());
  CHECK_UNARY_FALSE(Side<3>::south().isHigherOnAxis());
  CHECK_UNARY(Side<3>::north().isHigherOnAxis());
  CHECK_UNARY_FALSE(Side<3>::bottom().isHigherOnAxis());
  CHECK_UNARY(Side<3>::top().isHigherOnAxis());
}
TEST_CASE("Test iterator for Side<1>")
{
  auto iter = Side<1>::getValues().begin();
  CHECK_EQ(iter, Side<1>::getValues().begin());
  CHECK_NE(iter, Side<1>::getValues().end());
  CHECK_EQ(*iter, Side<1>::west());
  ++iter;
  CHECK_UNARY(iter->isHigherOnAxis());
  CHECK_EQ(*iter, Side<1>::east());
  ++iter;
  CHECK_EQ(*iter, Side<1>::null());
  CHECK_EQ(iter, Side<1>::getValues().end());
}
TEST_CASE("Test iterator for Side<2>")
{
  auto iter = Side<2>::getValues().begin();
  CHECK_EQ(iter, Side<2>::getValues().begin());
  CHECK_NE(iter, Side<2>::getValues().end());
  CHECK_EQ(*iter, Side<2>::west());
  ++iter;
  CHECK_UNARY(iter->isHigherOnAxis());
  CHECK_EQ(*iter, Side<2>::east());
  ++iter;
  CHECK_EQ(*iter, Side<2>::south());
  ++iter;
  CHECK_EQ(*iter, Side<2>::north());
  ++iter;
  CHECK_EQ(*iter, Side<2>::null());
  CHECK_EQ(iter, Side<2>::getValues().end());
}
TEST_CASE("Test iterator for Side<3>")
{
  auto iter = Side<3>::getValues().begin();
  CHECK_EQ(iter, Side<3>::getValues().begin());
  CHECK_NE(iter, Side<3>::getValues().end());
  CHECK_EQ(*iter, Side<3>::west());
  ++iter;
  CHECK_UNARY(iter->isHigherOnAxis());
  CHECK_EQ(*iter, Side<3>::east());
  ++iter;
  CHECK_EQ(*iter, Side<3>::south());
  ++iter;
  CHECK_EQ(*iter, Side<3>::north());
  ++iter;
  CHECK_EQ(*iter, Side<3>::bottom());
  ++iter;
  CHECK_EQ(*iter, Side<3>::top());
  ++iter;
  CHECK_EQ(*iter, Side<3>::null());
  CHECK_EQ(iter, Side<3>::getValues().end());
}
TEST_CASE("Test ostream for Side<1>")
{
  stringstream ss;
  ss << Side<1>::west();
  CHECK_EQ(ss.str(), "Side<1>::west()");
  ss.str("");
  ss << Side<1>::east();
  CHECK_EQ(ss.str(), "Side<1>::east()");
  ss.str("");
  ss << Side<1>::null();
  CHECK_EQ(ss.str(), "Side<1>::null()");
  ss.str("");
  ss << Side<1>(13);
  CHECK_EQ(ss.str(), "Side<1> undefined value: 13");
}
TEST_CASE("Test ostream for Side<2>")
{
  stringstream ss;
  ss << Side<2>::east();
  CHECK_EQ(ss.str(), "Side<2>::east()");
  ss.str("");
  ss << Side<2>::west();
  CHECK_EQ(ss.str(), "Side<2>::west()");
  ss.str("");
  ss << Side<2>::south();
  CHECK_EQ(ss.str(), "Side<2>::south()");
  ss.str("");
  ss << Side<2>::north();
  CHECK_EQ(ss.str(), "Side<2>::north()");
  ss.str("");
  ss << Side<2>::null();
  CHECK_EQ(ss.str(), "Side<2>::null()");
  ss.str("");
  ss << Side<2>(13);
  CHECK_EQ(ss.str(), "Side<2> undefined value: 13");
}
TEST_CASE("Test ostream for Side<3>")
{
  stringstream ss;
  ss << Side<3>::east();
  CHECK_EQ(ss.str(), "Side<3>::east()");
  ss.str("");
  ss << Side<3>::west();
  CHECK_EQ(ss.str(), "Side<3>::west()");
  ss.str("");
  ss << Side<3>::south();
  CHECK_EQ(ss.str(), "Side<3>::south()");
  ss.str("");
  ss << Side<3>::north();
  CHECK_EQ(ss.str(), "Side<3>::north()");
  ss.str("");
  ss << Side<3>::bottom();
  CHECK_EQ(ss.str(), "Side<3>::bottom()");
  ss.str("");
  ss << Side<3>::top();
  CHECK_EQ(ss.str(), "Side<3>::top()");
  ss.str("");
  ss << Side<3>::null();
  CHECK_EQ(ss.str(), "Side<3>::null()");
  ss.str("");
  ss << Side<3>(13);
  CHECK_EQ(ss.str(), "Side<3> undefined value: 13");
}
TEST_CASE("== operator works for Side<1>")
{
  CHECK_UNARY_FALSE(Side<1>::east() == Side<1>::west());
  CHECK_EQ(Side<1>::east(), Side<1>::east());
  CHECK_UNARY_FALSE(Side<1>::east() == Side<1>::null());

  CHECK_EQ(Side<1>::west(), Side<1>::west());
  CHECK_UNARY_FALSE(Side<1>::west() == Side<1>::east());
  CHECK_UNARY_FALSE(Side<1>::east() == Side<1>::null());

  CHECK_UNARY_FALSE(Side<1>::null() == Side<1>::west());
  CHECK_UNARY_FALSE(Side<1>::null() == Side<1>::east());
  CHECK_EQ(Side<1>::null(), Side<1>::null());
}
TEST_CASE("== operator works for Side<2>")
{
  CHECK_EQ(Side<2>::west(), Side<2>::west());
  CHECK_UNARY_FALSE(Side<2>::west() == Side<2>::east());
  CHECK_UNARY_FALSE(Side<2>::west() == Side<2>::south());
  CHECK_UNARY_FALSE(Side<2>::west() == Side<2>::north());
  CHECK_UNARY_FALSE(Side<2>::west() == Side<2>::null());

  CHECK_UNARY_FALSE(Side<2>::east() == Side<2>::west());
  CHECK_EQ(Side<2>::east(), Side<2>::east());
  CHECK_UNARY_FALSE(Side<2>::east() == Side<2>::south());
  CHECK_UNARY_FALSE(Side<2>::east() == Side<2>::north());
  CHECK_UNARY_FALSE(Side<2>::east() == Side<2>::null());

  CHECK_UNARY_FALSE(Side<2>::south() == Side<2>::west());
  CHECK_UNARY_FALSE(Side<2>::south() == Side<2>::east());
  CHECK_EQ(Side<2>::south(), Side<2>::south());
  CHECK_UNARY_FALSE(Side<2>::south() == Side<2>::north());
  CHECK_UNARY_FALSE(Side<2>::south() == Side<2>::null());

  CHECK_UNARY_FALSE(Side<2>::north() == Side<2>::west());
  CHECK_UNARY_FALSE(Side<2>::north() == Side<2>::east());
  CHECK_UNARY_FALSE(Side<2>::north() == Side<2>::south());
  CHECK_EQ(Side<2>::north(), Side<2>::north());
  CHECK_UNARY_FALSE(Side<2>::north() == Side<2>::null());

  CHECK_UNARY_FALSE(Side<2>::null() == Side<2>::west());
  CHECK_UNARY_FALSE(Side<2>::null() == Side<2>::east());
  CHECK_UNARY_FALSE(Side<2>::null() == Side<2>::south());
  CHECK_UNARY_FALSE(Side<2>::null() == Side<2>::north());
  CHECK_EQ(Side<2>::null(), Side<2>::null());
}
TEST_CASE("== operator works for Side<3>")
{
  CHECK_EQ(Side<3>::west(), Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::west() == Side<3>::east());
  CHECK_UNARY_FALSE(Side<3>::west() == Side<3>::south());
  CHECK_UNARY_FALSE(Side<3>::west() == Side<3>::north());
  CHECK_UNARY_FALSE(Side<3>::west() == Side<3>::bottom());
  CHECK_UNARY_FALSE(Side<3>::west() == Side<3>::top());
  CHECK_UNARY_FALSE(Side<3>::west() == Side<3>::null());

  CHECK_UNARY_FALSE(Side<3>::east() == Side<3>::west());
  CHECK_EQ(Side<3>::east(), Side<3>::east());
  CHECK_UNARY_FALSE(Side<3>::east() == Side<3>::south());
  CHECK_UNARY_FALSE(Side<3>::east() == Side<3>::north());
  CHECK_UNARY_FALSE(Side<3>::east() == Side<3>::bottom());
  CHECK_UNARY_FALSE(Side<3>::east() == Side<3>::top());
  CHECK_UNARY_FALSE(Side<3>::east() == Side<3>::null());

  CHECK_UNARY_FALSE(Side<3>::south() == Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::south() == Side<3>::east());
  CHECK_EQ(Side<3>::south(), Side<3>::south());
  CHECK_UNARY_FALSE(Side<3>::south() == Side<3>::north());
  CHECK_UNARY_FALSE(Side<3>::south() == Side<3>::bottom());
  CHECK_UNARY_FALSE(Side<3>::south() == Side<3>::top());
  CHECK_UNARY_FALSE(Side<3>::south() == Side<3>::null());

  CHECK_UNARY_FALSE(Side<3>::north() == Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::north() == Side<3>::east());
  CHECK_UNARY_FALSE(Side<3>::north() == Side<3>::south());
  CHECK_EQ(Side<3>::north(), Side<3>::north());
  CHECK_UNARY_FALSE(Side<3>::north() == Side<3>::bottom());
  CHECK_UNARY_FALSE(Side<3>::north() == Side<3>::top());
  CHECK_UNARY_FALSE(Side<3>::north() == Side<3>::null());

  CHECK_UNARY_FALSE(Side<3>::bottom() == Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::bottom() == Side<3>::east());
  CHECK_UNARY_FALSE(Side<3>::bottom() == Side<3>::south());
  CHECK_UNARY_FALSE(Side<3>::bottom() == Side<3>::north());
  CHECK_EQ(Side<3>::bottom(), Side<3>::bottom());
  CHECK_UNARY_FALSE(Side<3>::bottom() == Side<3>::top());
  CHECK_UNARY_FALSE(Side<3>::bottom() == Side<3>::null());

  CHECK_UNARY_FALSE(Side<3>::top() == Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::top() == Side<3>::east());
  CHECK_UNARY_FALSE(Side<3>::top() == Side<3>::south());
  CHECK_UNARY_FALSE(Side<3>::top() == Side<3>::north());
  CHECK_UNARY_FALSE(Side<3>::top() == Side<3>::bottom());
  CHECK_EQ(Side<3>::top(), Side<3>::top());
  CHECK_UNARY_FALSE(Side<3>::top() == Side<3>::null());

  CHECK_UNARY_FALSE(Side<3>::null() == Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::null() == Side<3>::east());
  CHECK_UNARY_FALSE(Side<3>::null() == Side<3>::south());
  CHECK_UNARY_FALSE(Side<3>::null() == Side<3>::north());
  CHECK_UNARY_FALSE(Side<3>::null() == Side<3>::bottom());
  CHECK_UNARY_FALSE(Side<3>::null() == Side<3>::top());
  CHECK_EQ(Side<3>::null(), Side<3>::null());
}
TEST_CASE("!= operator works for Side<1>")
{
  CHECK_UNARY_FALSE(Side<1>::west() != Side<1>::west());
  CHECK_NE(Side<1>::west(), Side<1>::east());
  CHECK_NE(Side<1>::west(), Side<1>::null());

  CHECK_NE(Side<1>::east(), Side<1>::west());
  CHECK_UNARY_FALSE(Side<1>::east() != Side<1>::east());
  CHECK_NE(Side<1>::east(), Side<1>::null());

  CHECK_NE(Side<1>::null(), Side<1>::west());
  CHECK_NE(Side<1>::null(), Side<1>::east());
  CHECK_UNARY_FALSE(Side<1>::null() != Side<1>::null());
}
TEST_CASE("!= operator works for Side<2>")
{
  CHECK_UNARY_FALSE(Side<2>::west() != Side<2>::west());
  CHECK_NE(Side<2>::west(), Side<2>::east());
  CHECK_NE(Side<2>::west(), Side<2>::south());
  CHECK_NE(Side<2>::west(), Side<2>::north());
  CHECK_NE(Side<2>::west(), Side<2>::null());

  CHECK_NE(Side<2>::east(), Side<2>::west());
  CHECK_UNARY_FALSE(Side<2>::east() != Side<2>::east());
  CHECK_NE(Side<2>::east(), Side<2>::south());
  CHECK_NE(Side<2>::east(), Side<2>::north());
  CHECK_NE(Side<2>::east(), Side<2>::null());

  CHECK_NE(Side<2>::south(), Side<2>::west());
  CHECK_NE(Side<2>::south(), Side<2>::east());
  CHECK_UNARY_FALSE(Side<2>::south() != Side<2>::south());
  CHECK_NE(Side<2>::south(), Side<2>::north());
  CHECK_NE(Side<2>::south(), Side<2>::null());

  CHECK_NE(Side<2>::north(), Side<2>::west());
  CHECK_NE(Side<2>::north(), Side<2>::east());
  CHECK_NE(Side<2>::north(), Side<2>::south());
  CHECK_UNARY_FALSE(Side<2>::north() != Side<2>::north());
  CHECK_NE(Side<2>::north(), Side<2>::null());

  CHECK_NE(Side<2>::null(), Side<2>::west());
  CHECK_NE(Side<2>::null(), Side<2>::east());
  CHECK_NE(Side<2>::null(), Side<2>::south());
  CHECK_NE(Side<2>::null(), Side<2>::north());
  CHECK_UNARY_FALSE(Side<2>::null() != Side<2>::null());
}
TEST_CASE("!= operator works for Side<3>")
{
  CHECK_UNARY_FALSE(Side<3>::west() != Side<3>::west());
  CHECK_NE(Side<3>::west(), Side<3>::east());
  CHECK_NE(Side<3>::west(), Side<3>::south());
  CHECK_NE(Side<3>::west(), Side<3>::north());
  CHECK_NE(Side<3>::west(), Side<3>::bottom());
  CHECK_NE(Side<3>::west(), Side<3>::top());
  CHECK_NE(Side<3>::west(), Side<3>::null());

  CHECK_NE(Side<3>::east(), Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::east() != Side<3>::east());
  CHECK_NE(Side<3>::east(), Side<3>::south());
  CHECK_NE(Side<3>::east(), Side<3>::north());
  CHECK_NE(Side<3>::east(), Side<3>::bottom());
  CHECK_NE(Side<3>::east(), Side<3>::top());
  CHECK_NE(Side<3>::east(), Side<3>::null());

  CHECK_NE(Side<3>::south(), Side<3>::west());
  CHECK_NE(Side<3>::south(), Side<3>::east());
  CHECK_UNARY_FALSE(Side<3>::south() != Side<3>::south());
  CHECK_NE(Side<3>::south(), Side<3>::north());
  CHECK_NE(Side<3>::south(), Side<3>::bottom());
  CHECK_NE(Side<3>::south(), Side<3>::top());
  CHECK_NE(Side<3>::south(), Side<3>::null());

  CHECK_NE(Side<3>::north(), Side<3>::west());
  CHECK_NE(Side<3>::north(), Side<3>::east());
  CHECK_NE(Side<3>::north(), Side<3>::south());
  CHECK_UNARY_FALSE(Side<3>::north() != Side<3>::north());
  CHECK_NE(Side<3>::north(), Side<3>::bottom());
  CHECK_NE(Side<3>::north(), Side<3>::top());
  CHECK_NE(Side<3>::north(), Side<3>::null());

  CHECK_NE(Side<3>::bottom(), Side<3>::west());
  CHECK_NE(Side<3>::bottom(), Side<3>::east());
  CHECK_NE(Side<3>::bottom(), Side<3>::south());
  CHECK_NE(Side<3>::bottom(), Side<3>::north());
  CHECK_UNARY_FALSE(Side<3>::bottom() != Side<3>::bottom());
  CHECK_NE(Side<3>::bottom(), Side<3>::top());
  CHECK_NE(Side<3>::bottom(), Side<3>::null());

  CHECK_NE(Side<3>::top(), Side<3>::west());
  CHECK_NE(Side<3>::top(), Side<3>::east());
  CHECK_NE(Side<3>::top(), Side<3>::south());
  CHECK_NE(Side<3>::top(), Side<3>::north());
  CHECK_NE(Side<3>::top(), Side<3>::bottom());
  CHECK_UNARY_FALSE(Side<3>::top() != Side<3>::top());
  CHECK_NE(Side<3>::top(), Side<3>::null());

  CHECK_NE(Side<3>::null(), Side<3>::west());
  CHECK_NE(Side<3>::null(), Side<3>::east());
  CHECK_NE(Side<3>::null(), Side<3>::south());
  CHECK_NE(Side<3>::null(), Side<3>::north());
  CHECK_NE(Side<3>::null(), Side<3>::bottom());
  CHECK_NE(Side<3>::null(), Side<3>::top());
  CHECK_UNARY_FALSE(Side<3>::null() != Side<3>::null());
}
TEST_CASE("< operator works for Side<1>")
{
  CHECK_UNARY_FALSE(Side<1>::west() < Side<1>::west());
  CHECK_LT(Side<1>::west(), Side<1>::east());
  CHECK_LT(Side<1>::west(), Side<1>::null());

  CHECK_UNARY_FALSE(Side<1>::east() < Side<1>::west());
  CHECK_UNARY_FALSE(Side<1>::east() < Side<1>::east());
  CHECK_LT(Side<1>::east(), Side<1>::null());

  CHECK_UNARY_FALSE(Side<1>::null() < Side<1>::west());
  CHECK_UNARY_FALSE(Side<1>::null() < Side<1>::east());
  CHECK_UNARY_FALSE(Side<1>::null() < Side<1>::null());
}
TEST_CASE("< operator works for Side<2>")
{
  CHECK_UNARY_FALSE(Side<2>::west() < Side<2>::west());
  CHECK_LT(Side<2>::west(), Side<2>::east());
  CHECK_LT(Side<2>::west(), Side<2>::south());
  CHECK_LT(Side<2>::west(), Side<2>::north());
  CHECK_LT(Side<2>::west(), Side<2>::null());

  CHECK_UNARY_FALSE(Side<2>::east() < Side<2>::west());
  CHECK_UNARY_FALSE(Side<2>::east() < Side<2>::east());
  CHECK_LT(Side<2>::east(), Side<2>::south());
  CHECK_LT(Side<2>::east(), Side<2>::north());
  CHECK_LT(Side<2>::east(), Side<2>::null());

  CHECK_UNARY_FALSE(Side<2>::south() < Side<2>::west());
  CHECK_UNARY_FALSE(Side<2>::south() < Side<2>::east());
  CHECK_UNARY_FALSE(Side<2>::south() < Side<2>::south());
  CHECK_LT(Side<2>::south(), Side<2>::north());
  CHECK_LT(Side<2>::south(), Side<2>::null());

  CHECK_UNARY_FALSE(Side<2>::north() < Side<2>::west());
  CHECK_UNARY_FALSE(Side<2>::north() < Side<2>::east());
  CHECK_UNARY_FALSE(Side<2>::north() < Side<2>::south());
  CHECK_UNARY_FALSE(Side<2>::north() < Side<2>::north());
  CHECK_LT(Side<2>::north(), Side<2>::null());

  CHECK_UNARY_FALSE(Side<2>::null() < Side<2>::west());
  CHECK_UNARY_FALSE(Side<2>::null() < Side<2>::east());
  CHECK_UNARY_FALSE(Side<2>::null() < Side<2>::south());
  CHECK_UNARY_FALSE(Side<2>::null() < Side<2>::north());
  CHECK_UNARY_FALSE(Side<2>::null() < Side<2>::null());
}
TEST_CASE("< operator works for Side<3>")
{
  CHECK_UNARY_FALSE(Side<3>::west() < Side<3>::west());
  CHECK_LT(Side<3>::west(), Side<3>::east());
  CHECK_LT(Side<3>::west(), Side<3>::south());
  CHECK_LT(Side<3>::west(), Side<3>::north());
  CHECK_LT(Side<3>::west(), Side<3>::bottom());
  CHECK_LT(Side<3>::west(), Side<3>::top());
  CHECK_LT(Side<3>::west(), Side<3>::null());

  CHECK_UNARY_FALSE(Side<3>::east() < Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::east() < Side<3>::east());
  CHECK_LT(Side<3>::east(), Side<3>::south());
  CHECK_LT(Side<3>::east(), Side<3>::north());
  CHECK_LT(Side<3>::east(), Side<3>::bottom());
  CHECK_LT(Side<3>::east(), Side<3>::top());
  CHECK_LT(Side<3>::east(), Side<3>::null());

  CHECK_UNARY_FALSE(Side<3>::south() < Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::south() < Side<3>::east());
  CHECK_UNARY_FALSE(Side<3>::south() < Side<3>::south());
  CHECK_LT(Side<3>::south(), Side<3>::north());
  CHECK_LT(Side<3>::south(), Side<3>::bottom());
  CHECK_LT(Side<3>::south(), Side<3>::top());
  CHECK_LT(Side<3>::south(), Side<3>::null());

  CHECK_UNARY_FALSE(Side<3>::north() < Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::north() < Side<3>::east());
  CHECK_UNARY_FALSE(Side<3>::north() < Side<3>::south());
  CHECK_UNARY_FALSE(Side<3>::north() < Side<3>::north());
  CHECK_LT(Side<3>::north(), Side<3>::bottom());
  CHECK_LT(Side<3>::north(), Side<3>::top());
  CHECK_LT(Side<3>::north(), Side<3>::null());

  CHECK_UNARY_FALSE(Side<3>::bottom() < Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::bottom() < Side<3>::east());
  CHECK_UNARY_FALSE(Side<3>::bottom() < Side<3>::south());
  CHECK_UNARY_FALSE(Side<3>::bottom() < Side<3>::north());
  CHECK_UNARY_FALSE(Side<3>::bottom() < Side<3>::bottom());
  CHECK_LT(Side<3>::bottom(), Side<3>::top());
  CHECK_LT(Side<3>::bottom(), Side<3>::null());

  CHECK_UNARY_FALSE(Side<3>::top() < Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::top() < Side<3>::east());
  CHECK_UNARY_FALSE(Side<3>::top() < Side<3>::south());
  CHECK_UNARY_FALSE(Side<3>::top() < Side<3>::north());
  CHECK_UNARY_FALSE(Side<3>::top() < Side<3>::bottom());
  CHECK_UNARY_FALSE(Side<3>::top() < Side<3>::top());
  CHECK_LT(Side<3>::top(), Side<3>::null());

  CHECK_UNARY_FALSE(Side<3>::null() < Side<3>::west());
  CHECK_UNARY_FALSE(Side<3>::null() < Side<3>::east());
  CHECK_UNARY_FALSE(Side<3>::null() < Side<3>::south());
  CHECK_UNARY_FALSE(Side<3>::null() < Side<3>::north());
  CHECK_UNARY_FALSE(Side<3>::null() < Side<3>::bottom());
  CHECK_UNARY_FALSE(Side<3>::null() < Side<3>::top());
  CHECK_UNARY_FALSE(Side<3>::null() < Side<3>::null());
}
TEST_CASE("LowerSideOnAxis works for Side<3>")
{
  CHECK_EQ(LowerSideOnAxis<3>(0), Side<3>::west());
  CHECK_EQ(LowerSideOnAxis<3>(1), Side<3>::south());
  CHECK_EQ(LowerSideOnAxis<3>(2), Side<3>::bottom());
}
TEST_CASE("LowerSideOnAxis works for Side<2>")
{
  CHECK_EQ(LowerSideOnAxis<2>(0), Side<2>::west());
  CHECK_EQ(LowerSideOnAxis<2>(1), Side<2>::south());
}
TEST_CASE("LowerSideOnAxis works for Side<1>")
{
  CHECK_EQ(LowerSideOnAxis<1>(0), Side<1>::west());
}
TEST_CASE("HigherSideOnAxis works for Side<3>")
{
  CHECK_EQ(HigherSideOnAxis<3>(0), Side<3>::east());
  CHECK_EQ(HigherSideOnAxis<3>(1), Side<3>::north());
  CHECK_EQ(HigherSideOnAxis<3>(2), Side<3>::top());
}
TEST_CASE("HigherSideOnAxis works for Side<2>")
{
  CHECK_EQ(HigherSideOnAxis<2>(0), Side<2>::east());
  CHECK_EQ(HigherSideOnAxis<2>(1), Side<2>::north());
}
TEST_CASE("HigherSideOnAxis works for Side<1>")
{
  CHECK_EQ(HigherSideOnAxis<1>(0), Side<1>::east());
}
