/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2021      Scott Aiton
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

TEST_CASE("Edge num_edges")
{
  CHECK(Edge::number_of == 12);
}
TEST_CASE("Edge dimensionality")
{
  CHECK(Edge::dimensionality == 1);
}
TEST_CASE("Edge unsigned char constructor works")
{
  Edge o(13);
  CHECK(o.getIndex() == 13);
}
TEST_CASE("Edge Default constructor works")
{
  Edge o;
  CHECK(o == Edge::null());
}
TEST_CASE("Edge named constructors give expected index values")
{
  CHECK(Edge::bs().getIndex() == 0);
  CHECK(Edge::bn().getIndex() == 1);
  CHECK(Edge::ts().getIndex() == 2);
  CHECK(Edge::tn().getIndex() == 3);
  CHECK(Edge::bw().getIndex() == 4);
  CHECK(Edge::be().getIndex() == 5);
  CHECK(Edge::tw().getIndex() == 6);
  CHECK(Edge::te().getIndex() == 7);
  CHECK(Edge::sw().getIndex() == 8);
  CHECK(Edge::se().getIndex() == 9);
  CHECK(Edge::nw().getIndex() == 10);
  CHECK(Edge::ne().getIndex() == 11);
  CHECK(Edge::null().getIndex() == 12);
}
TEST_CASE("Edge opposite")
{
  CHECK(Edge::bs().opposite() == Edge::tn());
  CHECK(Edge::tn().opposite() == Edge::bs());
  CHECK(Edge::bn().opposite() == Edge::ts());
  CHECK(Edge::ts().opposite() == Edge::bn());
  CHECK(Edge::bw().opposite() == Edge::te());
  CHECK(Edge::te().opposite() == Edge::bw());
  CHECK(Edge::be().opposite() == Edge::tw());
  CHECK(Edge::tw().opposite() == Edge::be());
  CHECK(Edge::sw().opposite() == Edge::ne());
  CHECK(Edge::ne().opposite() == Edge::sw());
  CHECK(Edge::se().opposite() == Edge::nw());
  CHECK(Edge::nw().opposite() == Edge::se());
}
TEST_CASE("Edge getSides")
{
  CHECK(Edge::bs().getSides()[0] == Side<3>::south());
  CHECK(Edge::bs().getSides()[1] == Side<3>::bottom());

  CHECK(Edge::tn().getSides()[0] == Side<3>::north());
  CHECK(Edge::tn().getSides()[1] == Side<3>::top());

  CHECK(Edge::bn().getSides()[0] == Side<3>::north());
  CHECK(Edge::bn().getSides()[1] == Side<3>::bottom());

  CHECK(Edge::ts().getSides()[0] == Side<3>::south());
  CHECK(Edge::ts().getSides()[1] == Side<3>::top());

  CHECK(Edge::bw().getSides()[0] == Side<3>::west());
  CHECK(Edge::bw().getSides()[1] == Side<3>::bottom());

  CHECK(Edge::te().getSides()[0] == Side<3>::east());
  CHECK(Edge::te().getSides()[1] == Side<3>::top());

  CHECK(Edge::be().getSides()[0] == Side<3>::east());
  CHECK(Edge::be().getSides()[1] == Side<3>::bottom());

  CHECK(Edge::tw().getSides()[0] == Side<3>::west());
  CHECK(Edge::tw().getSides()[1] == Side<3>::top());

  CHECK(Edge::sw().getSides()[0] == Side<3>::west());
  CHECK(Edge::sw().getSides()[1] == Side<3>::south());

  CHECK(Edge::ne().getSides()[0] == Side<3>::east());
  CHECK(Edge::ne().getSides()[1] == Side<3>::north());

  CHECK(Edge::se().getSides()[0] == Side<3>::east());
  CHECK(Edge::se().getSides()[1] == Side<3>::south());

  CHECK(Edge::nw().getSides()[0] == Side<3>::west());
  CHECK(Edge::nw().getSides()[1] == Side<3>::north());
}
TEST_CASE("Edge ==")
{
  CHECK(Edge::bs() == Edge::bs());
  CHECK_FALSE(Edge::bs() == Edge::tn());
  CHECK_FALSE(Edge::bs() == Edge::bn());
  CHECK_FALSE(Edge::bs() == Edge::ts());
  CHECK_FALSE(Edge::bs() == Edge::bw());
  CHECK_FALSE(Edge::bs() == Edge::te());
  CHECK_FALSE(Edge::bs() == Edge::be());
  CHECK_FALSE(Edge::bs() == Edge::tw());
  CHECK_FALSE(Edge::bs() == Edge::sw());
  CHECK_FALSE(Edge::bs() == Edge::ne());
  CHECK_FALSE(Edge::bs() == Edge::se());
  CHECK_FALSE(Edge::bs() == Edge::nw());
  CHECK_FALSE(Edge::bs() == Edge::null());

  CHECK_FALSE(Edge::tn() == Edge::bs());
  CHECK(Edge::tn() == Edge::tn());
  CHECK_FALSE(Edge::tn() == Edge::bn());
  CHECK_FALSE(Edge::tn() == Edge::ts());
  CHECK_FALSE(Edge::tn() == Edge::bw());
  CHECK_FALSE(Edge::tn() == Edge::te());
  CHECK_FALSE(Edge::tn() == Edge::be());
  CHECK_FALSE(Edge::tn() == Edge::tw());
  CHECK_FALSE(Edge::tn() == Edge::sw());
  CHECK_FALSE(Edge::tn() == Edge::ne());
  CHECK_FALSE(Edge::tn() == Edge::se());
  CHECK_FALSE(Edge::tn() == Edge::nw());
  CHECK_FALSE(Edge::tn() == Edge::null());

  CHECK_FALSE(Edge::bn() == Edge::bs());
  CHECK_FALSE(Edge::bn() == Edge::tn());
  CHECK(Edge::bn() == Edge::bn());
  CHECK_FALSE(Edge::bn() == Edge::ts());
  CHECK_FALSE(Edge::bn() == Edge::bw());
  CHECK_FALSE(Edge::bn() == Edge::te());
  CHECK_FALSE(Edge::bn() == Edge::be());
  CHECK_FALSE(Edge::bn() == Edge::tw());
  CHECK_FALSE(Edge::bn() == Edge::sw());
  CHECK_FALSE(Edge::bn() == Edge::ne());
  CHECK_FALSE(Edge::bn() == Edge::se());
  CHECK_FALSE(Edge::bn() == Edge::nw());
  CHECK_FALSE(Edge::bn() == Edge::null());

  CHECK_FALSE(Edge::ts() == Edge::bs());
  CHECK_FALSE(Edge::ts() == Edge::tn());
  CHECK_FALSE(Edge::ts() == Edge::bn());
  CHECK(Edge::ts() == Edge::ts());
  CHECK_FALSE(Edge::ts() == Edge::bw());
  CHECK_FALSE(Edge::ts() == Edge::te());
  CHECK_FALSE(Edge::ts() == Edge::be());
  CHECK_FALSE(Edge::ts() == Edge::tw());
  CHECK_FALSE(Edge::ts() == Edge::sw());
  CHECK_FALSE(Edge::ts() == Edge::ne());
  CHECK_FALSE(Edge::ts() == Edge::se());
  CHECK_FALSE(Edge::ts() == Edge::nw());
  CHECK_FALSE(Edge::ts() == Edge::null());

  CHECK_FALSE(Edge::bw() == Edge::bs());
  CHECK_FALSE(Edge::bw() == Edge::tn());
  CHECK_FALSE(Edge::bw() == Edge::bn());
  CHECK_FALSE(Edge::bw() == Edge::ts());
  CHECK(Edge::bw() == Edge::bw());
  CHECK_FALSE(Edge::bw() == Edge::te());
  CHECK_FALSE(Edge::bw() == Edge::be());
  CHECK_FALSE(Edge::bw() == Edge::tw());
  CHECK_FALSE(Edge::bw() == Edge::sw());
  CHECK_FALSE(Edge::bw() == Edge::ne());
  CHECK_FALSE(Edge::bw() == Edge::se());
  CHECK_FALSE(Edge::bw() == Edge::nw());
  CHECK_FALSE(Edge::bw() == Edge::null());

  CHECK_FALSE(Edge::te() == Edge::bs());
  CHECK_FALSE(Edge::te() == Edge::tn());
  CHECK_FALSE(Edge::te() == Edge::bn());
  CHECK_FALSE(Edge::te() == Edge::ts());
  CHECK_FALSE(Edge::te() == Edge::bw());
  CHECK(Edge::te() == Edge::te());
  CHECK_FALSE(Edge::te() == Edge::be());
  CHECK_FALSE(Edge::te() == Edge::tw());
  CHECK_FALSE(Edge::te() == Edge::sw());
  CHECK_FALSE(Edge::te() == Edge::ne());
  CHECK_FALSE(Edge::te() == Edge::se());
  CHECK_FALSE(Edge::te() == Edge::nw());
  CHECK_FALSE(Edge::te() == Edge::null());

  CHECK_FALSE(Edge::be() == Edge::bs());
  CHECK_FALSE(Edge::be() == Edge::tn());
  CHECK_FALSE(Edge::be() == Edge::bn());
  CHECK_FALSE(Edge::be() == Edge::ts());
  CHECK_FALSE(Edge::be() == Edge::bw());
  CHECK_FALSE(Edge::be() == Edge::te());
  CHECK(Edge::be() == Edge::be());
  CHECK_FALSE(Edge::be() == Edge::tw());
  CHECK_FALSE(Edge::be() == Edge::sw());
  CHECK_FALSE(Edge::be() == Edge::ne());
  CHECK_FALSE(Edge::be() == Edge::se());
  CHECK_FALSE(Edge::be() == Edge::nw());
  CHECK_FALSE(Edge::be() == Edge::null());

  CHECK_FALSE(Edge::tw() == Edge::bs());
  CHECK_FALSE(Edge::tw() == Edge::tn());
  CHECK_FALSE(Edge::tw() == Edge::bn());
  CHECK_FALSE(Edge::tw() == Edge::ts());
  CHECK_FALSE(Edge::tw() == Edge::bw());
  CHECK_FALSE(Edge::tw() == Edge::te());
  CHECK_FALSE(Edge::tw() == Edge::be());
  CHECK(Edge::tw() == Edge::tw());
  CHECK_FALSE(Edge::tw() == Edge::sw());
  CHECK_FALSE(Edge::tw() == Edge::ne());
  CHECK_FALSE(Edge::tw() == Edge::se());
  CHECK_FALSE(Edge::tw() == Edge::nw());
  CHECK_FALSE(Edge::tw() == Edge::null());

  CHECK_FALSE(Edge::sw() == Edge::bs());
  CHECK_FALSE(Edge::sw() == Edge::tn());
  CHECK_FALSE(Edge::sw() == Edge::bn());
  CHECK_FALSE(Edge::sw() == Edge::ts());
  CHECK_FALSE(Edge::sw() == Edge::bw());
  CHECK_FALSE(Edge::sw() == Edge::te());
  CHECK_FALSE(Edge::sw() == Edge::be());
  CHECK_FALSE(Edge::sw() == Edge::tw());
  CHECK(Edge::sw() == Edge::sw());
  CHECK_FALSE(Edge::sw() == Edge::ne());
  CHECK_FALSE(Edge::sw() == Edge::se());
  CHECK_FALSE(Edge::sw() == Edge::nw());
  CHECK_FALSE(Edge::sw() == Edge::null());

  CHECK_FALSE(Edge::ne() == Edge::bs());
  CHECK_FALSE(Edge::ne() == Edge::tn());
  CHECK_FALSE(Edge::ne() == Edge::bn());
  CHECK_FALSE(Edge::ne() == Edge::ts());
  CHECK_FALSE(Edge::ne() == Edge::bw());
  CHECK_FALSE(Edge::ne() == Edge::te());
  CHECK_FALSE(Edge::ne() == Edge::be());
  CHECK_FALSE(Edge::ne() == Edge::tw());
  CHECK_FALSE(Edge::ne() == Edge::sw());
  CHECK(Edge::ne() == Edge::ne());
  CHECK_FALSE(Edge::ne() == Edge::se());
  CHECK_FALSE(Edge::ne() == Edge::nw());
  CHECK_FALSE(Edge::ne() == Edge::null());

  CHECK_FALSE(Edge::se() == Edge::bs());
  CHECK_FALSE(Edge::se() == Edge::tn());
  CHECK_FALSE(Edge::se() == Edge::bn());
  CHECK_FALSE(Edge::se() == Edge::ts());
  CHECK_FALSE(Edge::se() == Edge::bw());
  CHECK_FALSE(Edge::se() == Edge::te());
  CHECK_FALSE(Edge::se() == Edge::be());
  CHECK_FALSE(Edge::se() == Edge::tw());
  CHECK_FALSE(Edge::se() == Edge::sw());
  CHECK_FALSE(Edge::se() == Edge::ne());
  CHECK(Edge::se() == Edge::se());
  CHECK_FALSE(Edge::se() == Edge::nw());
  CHECK_FALSE(Edge::se() == Edge::null());

  CHECK_FALSE(Edge::nw() == Edge::bs());
  CHECK_FALSE(Edge::nw() == Edge::tn());
  CHECK_FALSE(Edge::nw() == Edge::bn());
  CHECK_FALSE(Edge::nw() == Edge::ts());
  CHECK_FALSE(Edge::nw() == Edge::bw());
  CHECK_FALSE(Edge::nw() == Edge::te());
  CHECK_FALSE(Edge::nw() == Edge::be());
  CHECK_FALSE(Edge::nw() == Edge::tw());
  CHECK_FALSE(Edge::nw() == Edge::sw());
  CHECK_FALSE(Edge::nw() == Edge::ne());
  CHECK_FALSE(Edge::nw() == Edge::se());
  CHECK(Edge::nw() == Edge::nw());
  CHECK_FALSE(Edge::nw() == Edge::null());

  CHECK_FALSE(Edge::null() == Edge::bs());
  CHECK_FALSE(Edge::null() == Edge::tn());
  CHECK_FALSE(Edge::null() == Edge::bn());
  CHECK_FALSE(Edge::null() == Edge::ts());
  CHECK_FALSE(Edge::null() == Edge::bw());
  CHECK_FALSE(Edge::null() == Edge::te());
  CHECK_FALSE(Edge::null() == Edge::be());
  CHECK_FALSE(Edge::null() == Edge::tw());
  CHECK_FALSE(Edge::null() == Edge::sw());
  CHECK_FALSE(Edge::null() == Edge::ne());
  CHECK_FALSE(Edge::null() == Edge::se());
  CHECK_FALSE(Edge::null() == Edge::nw());
  CHECK(Edge::null() == Edge::null());
}
TEST_CASE("Edge !=")
{
  CHECK_FALSE(Edge::bs() != Edge::bs());
  CHECK(Edge::bs() != Edge::tn());
  CHECK(Edge::bs() != Edge::bn());
  CHECK(Edge::bs() != Edge::ts());
  CHECK(Edge::bs() != Edge::bw());
  CHECK(Edge::bs() != Edge::te());
  CHECK(Edge::bs() != Edge::be());
  CHECK(Edge::bs() != Edge::tw());
  CHECK(Edge::bs() != Edge::sw());
  CHECK(Edge::bs() != Edge::ne());
  CHECK(Edge::bs() != Edge::se());
  CHECK(Edge::bs() != Edge::nw());
  CHECK(Edge::bs() != Edge::null());

  CHECK(Edge::tn() != Edge::bs());
  CHECK_FALSE(Edge::tn() != Edge::tn());
  CHECK(Edge::tn() != Edge::bn());
  CHECK(Edge::tn() != Edge::ts());
  CHECK(Edge::tn() != Edge::bw());
  CHECK(Edge::tn() != Edge::te());
  CHECK(Edge::tn() != Edge::be());
  CHECK(Edge::tn() != Edge::tw());
  CHECK(Edge::tn() != Edge::sw());
  CHECK(Edge::tn() != Edge::ne());
  CHECK(Edge::tn() != Edge::se());
  CHECK(Edge::tn() != Edge::nw());
  CHECK(Edge::tn() != Edge::null());

  CHECK(Edge::bn() != Edge::bs());
  CHECK(Edge::bn() != Edge::tn());
  CHECK_FALSE(Edge::bn() != Edge::bn());
  CHECK(Edge::bn() != Edge::ts());
  CHECK(Edge::bn() != Edge::bw());
  CHECK(Edge::bn() != Edge::te());
  CHECK(Edge::bn() != Edge::be());
  CHECK(Edge::bn() != Edge::tw());
  CHECK(Edge::bn() != Edge::sw());
  CHECK(Edge::bn() != Edge::ne());
  CHECK(Edge::bn() != Edge::se());
  CHECK(Edge::bn() != Edge::nw());
  CHECK(Edge::bn() != Edge::null());

  CHECK(Edge::ts() != Edge::bs());
  CHECK(Edge::ts() != Edge::tn());
  CHECK(Edge::ts() != Edge::bn());
  CHECK_FALSE(Edge::ts() != Edge::ts());
  CHECK(Edge::ts() != Edge::bw());
  CHECK(Edge::ts() != Edge::te());
  CHECK(Edge::ts() != Edge::be());
  CHECK(Edge::ts() != Edge::tw());
  CHECK(Edge::ts() != Edge::sw());
  CHECK(Edge::ts() != Edge::ne());
  CHECK(Edge::ts() != Edge::se());
  CHECK(Edge::ts() != Edge::nw());
  CHECK(Edge::ts() != Edge::null());

  CHECK(Edge::bw() != Edge::bs());
  CHECK(Edge::bw() != Edge::tn());
  CHECK(Edge::bw() != Edge::bn());
  CHECK(Edge::bw() != Edge::ts());
  CHECK_FALSE(Edge::bw() != Edge::bw());
  CHECK(Edge::bw() != Edge::te());
  CHECK(Edge::bw() != Edge::be());
  CHECK(Edge::bw() != Edge::tw());
  CHECK(Edge::bw() != Edge::sw());
  CHECK(Edge::bw() != Edge::ne());
  CHECK(Edge::bw() != Edge::se());
  CHECK(Edge::bw() != Edge::nw());
  CHECK(Edge::bw() != Edge::null());

  CHECK(Edge::te() != Edge::bs());
  CHECK(Edge::te() != Edge::tn());
  CHECK(Edge::te() != Edge::bn());
  CHECK(Edge::te() != Edge::ts());
  CHECK(Edge::te() != Edge::bw());
  CHECK_FALSE(Edge::te() != Edge::te());
  CHECK(Edge::te() != Edge::be());
  CHECK(Edge::te() != Edge::tw());
  CHECK(Edge::te() != Edge::sw());
  CHECK(Edge::te() != Edge::ne());
  CHECK(Edge::te() != Edge::se());
  CHECK(Edge::te() != Edge::nw());
  CHECK(Edge::te() != Edge::null());

  CHECK(Edge::be() != Edge::bs());
  CHECK(Edge::be() != Edge::tn());
  CHECK(Edge::be() != Edge::bn());
  CHECK(Edge::be() != Edge::ts());
  CHECK(Edge::be() != Edge::bw());
  CHECK(Edge::be() != Edge::te());
  CHECK_FALSE(Edge::be() != Edge::be());
  CHECK(Edge::be() != Edge::tw());
  CHECK(Edge::be() != Edge::sw());
  CHECK(Edge::be() != Edge::ne());
  CHECK(Edge::be() != Edge::se());
  CHECK(Edge::be() != Edge::nw());
  CHECK(Edge::be() != Edge::null());

  CHECK(Edge::tw() != Edge::bs());
  CHECK(Edge::tw() != Edge::tn());
  CHECK(Edge::tw() != Edge::bn());
  CHECK(Edge::tw() != Edge::ts());
  CHECK(Edge::tw() != Edge::bw());
  CHECK(Edge::tw() != Edge::te());
  CHECK(Edge::tw() != Edge::be());
  CHECK_FALSE(Edge::tw() != Edge::tw());
  CHECK(Edge::tw() != Edge::sw());
  CHECK(Edge::tw() != Edge::ne());
  CHECK(Edge::tw() != Edge::se());
  CHECK(Edge::tw() != Edge::nw());
  CHECK(Edge::tw() != Edge::null());

  CHECK(Edge::sw() != Edge::bs());
  CHECK(Edge::sw() != Edge::tn());
  CHECK(Edge::sw() != Edge::bn());
  CHECK(Edge::sw() != Edge::ts());
  CHECK(Edge::sw() != Edge::bw());
  CHECK(Edge::sw() != Edge::te());
  CHECK(Edge::sw() != Edge::be());
  CHECK(Edge::sw() != Edge::tw());
  CHECK_FALSE(Edge::sw() != Edge::sw());
  CHECK(Edge::sw() != Edge::ne());
  CHECK(Edge::sw() != Edge::se());
  CHECK(Edge::sw() != Edge::nw());
  CHECK(Edge::sw() != Edge::null());

  CHECK(Edge::ne() != Edge::bs());
  CHECK(Edge::ne() != Edge::tn());
  CHECK(Edge::ne() != Edge::bn());
  CHECK(Edge::ne() != Edge::ts());
  CHECK(Edge::ne() != Edge::bw());
  CHECK(Edge::ne() != Edge::te());
  CHECK(Edge::ne() != Edge::be());
  CHECK(Edge::ne() != Edge::tw());
  CHECK(Edge::ne() != Edge::sw());
  CHECK_FALSE(Edge::ne() != Edge::ne());
  CHECK(Edge::ne() != Edge::se());
  CHECK(Edge::ne() != Edge::nw());
  CHECK(Edge::ne() != Edge::null());

  CHECK(Edge::se() != Edge::bs());
  CHECK(Edge::se() != Edge::tn());
  CHECK(Edge::se() != Edge::bn());
  CHECK(Edge::se() != Edge::ts());
  CHECK(Edge::se() != Edge::bw());
  CHECK(Edge::se() != Edge::te());
  CHECK(Edge::se() != Edge::be());
  CHECK(Edge::se() != Edge::tw());
  CHECK(Edge::se() != Edge::sw());
  CHECK(Edge::se() != Edge::ne());
  CHECK_FALSE(Edge::se() != Edge::se());
  CHECK(Edge::se() != Edge::nw());
  CHECK(Edge::se() != Edge::null());

  CHECK(Edge::nw() != Edge::bs());
  CHECK(Edge::nw() != Edge::tn());
  CHECK(Edge::nw() != Edge::bn());
  CHECK(Edge::nw() != Edge::ts());
  CHECK(Edge::nw() != Edge::bw());
  CHECK(Edge::nw() != Edge::te());
  CHECK(Edge::nw() != Edge::be());
  CHECK(Edge::nw() != Edge::tw());
  CHECK(Edge::nw() != Edge::sw());
  CHECK(Edge::nw() != Edge::ne());
  CHECK(Edge::nw() != Edge::se());
  CHECK_FALSE(Edge::nw() != Edge::nw());
  CHECK(Edge::nw() != Edge::null());

  CHECK(Edge::null() != Edge::bs());
  CHECK(Edge::null() != Edge::tn());
  CHECK(Edge::null() != Edge::bn());
  CHECK(Edge::null() != Edge::ts());
  CHECK(Edge::null() != Edge::bw());
  CHECK(Edge::null() != Edge::te());
  CHECK(Edge::null() != Edge::be());
  CHECK(Edge::null() != Edge::tw());
  CHECK(Edge::null() != Edge::sw());
  CHECK(Edge::null() != Edge::ne());
  CHECK(Edge::null() != Edge::se());
  CHECK(Edge::null() != Edge::nw());
  CHECK_FALSE(Edge::null() != Edge::null());
}
TEST_CASE("Edge <")
{
  for (unsigned char val = 0; val <= Edge::number_of; val++) {
    Edge edge(val);
    for (unsigned char other_val = 0; other_val <= Edge::number_of; other_val++) {
      Edge other_edge(other_val);
      CHECK((edge < other_edge) == (val < other_val));
    }
  }
}
TEST_CASE("Test ostream for Edge")
{
  stringstream ss;
  ss << Edge::bs();
  CHECK(ss.str() == "Edge::bs()");
  ss.str("");
  ss << Edge::tn();
  CHECK(ss.str() == "Edge::tn()");
  ss.str("");
  ss << Edge::bn();
  CHECK(ss.str() == "Edge::bn()");
  ss.str("");
  ss << Edge::ts();
  CHECK(ss.str() == "Edge::ts()");
  ss.str("");
  ss << Edge::bw();
  CHECK(ss.str() == "Edge::bw()");
  ss.str("");
  ss << Edge::te();
  CHECK(ss.str() == "Edge::te()");
  ss.str("");
  ss << Edge::be();
  CHECK(ss.str() == "Edge::be()");
  ss.str("");
  ss << Edge::tw();
  CHECK(ss.str() == "Edge::tw()");
  ss.str("");
  ss << Edge::sw();
  CHECK(ss.str() == "Edge::sw()");
  ss.str("");
  ss << Edge::ne();
  CHECK(ss.str() == "Edge::ne()");
  ss.str("");
  ss << Edge::se();
  CHECK(ss.str() == "Edge::se()");
  ss.str("");
  ss << Edge::nw();
  CHECK(ss.str() == "Edge::nw()");
  ss.str("");
  ss << Edge::null();
  CHECK(ss.str() == "Edge::null()");
  ss.str("");
  ss << Edge(64);
  CHECK(ss.str() == "Edge invalid value: 64");
}
TEST_CASE("Test iterator for Edge")
{
  auto iter = Edge::getValues().begin();
  CHECK(iter == Edge::getValues().begin());
  CHECK(iter != Edge::getValues().end());
  CHECK(*iter == Edge::bs());
  ++iter;
  CHECK(iter->getIndex() == 1);
  CHECK(*iter == Edge::bn());
  ++iter;
  CHECK(*iter == Edge::ts());
  ++iter;
  CHECK(*iter == Edge::tn());
  ++iter;
  CHECK(*iter == Edge::bw());
  ++iter;
  CHECK(*iter == Edge::be());
  ++iter;
  CHECK(*iter == Edge::tw());
  ++iter;
  CHECK(*iter == Edge::te());
  ++iter;
  CHECK(*iter == Edge::sw());
  ++iter;
  CHECK(*iter == Edge::se());
  ++iter;
  CHECK(*iter == Edge::nw());
  ++iter;
  CHECK(*iter == Edge::ne());
  ++iter;
  CHECK(*iter == Edge::null());
  CHECK(iter == Edge::getValues().end());
}
