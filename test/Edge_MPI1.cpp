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
  CHECK_EQ(Edge::number_of, 12);
}
TEST_CASE("Edge dimensionality")
{
  CHECK_EQ(Edge::dimensionality, 1);
}
TEST_CASE("Edge unsigned char constructor works")
{
  Edge o(13);
  CHECK_EQ(o.getIndex(), 13);
}
TEST_CASE("Edge Default constructor works")
{
  Edge o;
  CHECK_EQ(o, Edge::null());
}
TEST_CASE("Edge named constructors give expected index values")
{
  CHECK_EQ(Edge::bs().getIndex(), 0);
  CHECK_EQ(Edge::bn().getIndex(), 1);
  CHECK_EQ(Edge::ts().getIndex(), 2);
  CHECK_EQ(Edge::tn().getIndex(), 3);
  CHECK_EQ(Edge::bw().getIndex(), 4);
  CHECK_EQ(Edge::be().getIndex(), 5);
  CHECK_EQ(Edge::tw().getIndex(), 6);
  CHECK_EQ(Edge::te().getIndex(), 7);
  CHECK_EQ(Edge::sw().getIndex(), 8);
  CHECK_EQ(Edge::se().getIndex(), 9);
  CHECK_EQ(Edge::nw().getIndex(), 10);
  CHECK_EQ(Edge::ne().getIndex(), 11);
  CHECK_EQ(Edge::null().getIndex(), 12);
}
TEST_CASE("Edge opposite")
{
  CHECK_EQ(Edge::bs().opposite(), Edge::tn());
  CHECK_EQ(Edge::tn().opposite(), Edge::bs());
  CHECK_EQ(Edge::bn().opposite(), Edge::ts());
  CHECK_EQ(Edge::ts().opposite(), Edge::bn());
  CHECK_EQ(Edge::bw().opposite(), Edge::te());
  CHECK_EQ(Edge::te().opposite(), Edge::bw());
  CHECK_EQ(Edge::be().opposite(), Edge::tw());
  CHECK_EQ(Edge::tw().opposite(), Edge::be());
  CHECK_EQ(Edge::sw().opposite(), Edge::ne());
  CHECK_EQ(Edge::ne().opposite(), Edge::sw());
  CHECK_EQ(Edge::se().opposite(), Edge::nw());
  CHECK_EQ(Edge::nw().opposite(), Edge::se());
}
TEST_CASE("Edge getSides")
{
  CHECK_EQ(Edge::bs().getSides()[0], Side<3>::south());
  CHECK_EQ(Edge::bs().getSides()[1], Side<3>::bottom());

  CHECK_EQ(Edge::tn().getSides()[0], Side<3>::north());
  CHECK_EQ(Edge::tn().getSides()[1], Side<3>::top());

  CHECK_EQ(Edge::bn().getSides()[0], Side<3>::north());
  CHECK_EQ(Edge::bn().getSides()[1], Side<3>::bottom());

  CHECK_EQ(Edge::ts().getSides()[0], Side<3>::south());
  CHECK_EQ(Edge::ts().getSides()[1], Side<3>::top());

  CHECK_EQ(Edge::bw().getSides()[0], Side<3>::west());
  CHECK_EQ(Edge::bw().getSides()[1], Side<3>::bottom());

  CHECK_EQ(Edge::te().getSides()[0], Side<3>::east());
  CHECK_EQ(Edge::te().getSides()[1], Side<3>::top());

  CHECK_EQ(Edge::be().getSides()[0], Side<3>::east());
  CHECK_EQ(Edge::be().getSides()[1], Side<3>::bottom());

  CHECK_EQ(Edge::tw().getSides()[0], Side<3>::west());
  CHECK_EQ(Edge::tw().getSides()[1], Side<3>::top());

  CHECK_EQ(Edge::sw().getSides()[0], Side<3>::west());
  CHECK_EQ(Edge::sw().getSides()[1], Side<3>::south());

  CHECK_EQ(Edge::ne().getSides()[0], Side<3>::east());
  CHECK_EQ(Edge::ne().getSides()[1], Side<3>::north());

  CHECK_EQ(Edge::se().getSides()[0], Side<3>::east());
  CHECK_EQ(Edge::se().getSides()[1], Side<3>::south());

  CHECK_EQ(Edge::nw().getSides()[0], Side<3>::west());
  CHECK_EQ(Edge::nw().getSides()[1], Side<3>::north());
}
TEST_CASE("Edge ==")
{
  CHECK_EQ(Edge::bs(), Edge::bs());
  CHECK_UNARY_FALSE(Edge::bs() == Edge::tn());
  CHECK_UNARY_FALSE(Edge::bs() == Edge::bn());
  CHECK_UNARY_FALSE(Edge::bs() == Edge::ts());
  CHECK_UNARY_FALSE(Edge::bs() == Edge::bw());
  CHECK_UNARY_FALSE(Edge::bs() == Edge::te());
  CHECK_UNARY_FALSE(Edge::bs() == Edge::be());
  CHECK_UNARY_FALSE(Edge::bs() == Edge::tw());
  CHECK_UNARY_FALSE(Edge::bs() == Edge::sw());
  CHECK_UNARY_FALSE(Edge::bs() == Edge::ne());
  CHECK_UNARY_FALSE(Edge::bs() == Edge::se());
  CHECK_UNARY_FALSE(Edge::bs() == Edge::nw());
  CHECK_UNARY_FALSE(Edge::bs() == Edge::null());

  CHECK_UNARY_FALSE(Edge::tn() == Edge::bs());
  CHECK_EQ(Edge::tn(), Edge::tn());
  CHECK_UNARY_FALSE(Edge::tn() == Edge::bn());
  CHECK_UNARY_FALSE(Edge::tn() == Edge::ts());
  CHECK_UNARY_FALSE(Edge::tn() == Edge::bw());
  CHECK_UNARY_FALSE(Edge::tn() == Edge::te());
  CHECK_UNARY_FALSE(Edge::tn() == Edge::be());
  CHECK_UNARY_FALSE(Edge::tn() == Edge::tw());
  CHECK_UNARY_FALSE(Edge::tn() == Edge::sw());
  CHECK_UNARY_FALSE(Edge::tn() == Edge::ne());
  CHECK_UNARY_FALSE(Edge::tn() == Edge::se());
  CHECK_UNARY_FALSE(Edge::tn() == Edge::nw());
  CHECK_UNARY_FALSE(Edge::tn() == Edge::null());

  CHECK_UNARY_FALSE(Edge::bn() == Edge::bs());
  CHECK_UNARY_FALSE(Edge::bn() == Edge::tn());
  CHECK_EQ(Edge::bn(), Edge::bn());
  CHECK_UNARY_FALSE(Edge::bn() == Edge::ts());
  CHECK_UNARY_FALSE(Edge::bn() == Edge::bw());
  CHECK_UNARY_FALSE(Edge::bn() == Edge::te());
  CHECK_UNARY_FALSE(Edge::bn() == Edge::be());
  CHECK_UNARY_FALSE(Edge::bn() == Edge::tw());
  CHECK_UNARY_FALSE(Edge::bn() == Edge::sw());
  CHECK_UNARY_FALSE(Edge::bn() == Edge::ne());
  CHECK_UNARY_FALSE(Edge::bn() == Edge::se());
  CHECK_UNARY_FALSE(Edge::bn() == Edge::nw());
  CHECK_UNARY_FALSE(Edge::bn() == Edge::null());

  CHECK_UNARY_FALSE(Edge::ts() == Edge::bs());
  CHECK_UNARY_FALSE(Edge::ts() == Edge::tn());
  CHECK_UNARY_FALSE(Edge::ts() == Edge::bn());
  CHECK_EQ(Edge::ts(), Edge::ts());
  CHECK_UNARY_FALSE(Edge::ts() == Edge::bw());
  CHECK_UNARY_FALSE(Edge::ts() == Edge::te());
  CHECK_UNARY_FALSE(Edge::ts() == Edge::be());
  CHECK_UNARY_FALSE(Edge::ts() == Edge::tw());
  CHECK_UNARY_FALSE(Edge::ts() == Edge::sw());
  CHECK_UNARY_FALSE(Edge::ts() == Edge::ne());
  CHECK_UNARY_FALSE(Edge::ts() == Edge::se());
  CHECK_UNARY_FALSE(Edge::ts() == Edge::nw());
  CHECK_UNARY_FALSE(Edge::ts() == Edge::null());

  CHECK_UNARY_FALSE(Edge::bw() == Edge::bs());
  CHECK_UNARY_FALSE(Edge::bw() == Edge::tn());
  CHECK_UNARY_FALSE(Edge::bw() == Edge::bn());
  CHECK_UNARY_FALSE(Edge::bw() == Edge::ts());
  CHECK_EQ(Edge::bw(), Edge::bw());
  CHECK_UNARY_FALSE(Edge::bw() == Edge::te());
  CHECK_UNARY_FALSE(Edge::bw() == Edge::be());
  CHECK_UNARY_FALSE(Edge::bw() == Edge::tw());
  CHECK_UNARY_FALSE(Edge::bw() == Edge::sw());
  CHECK_UNARY_FALSE(Edge::bw() == Edge::ne());
  CHECK_UNARY_FALSE(Edge::bw() == Edge::se());
  CHECK_UNARY_FALSE(Edge::bw() == Edge::nw());
  CHECK_UNARY_FALSE(Edge::bw() == Edge::null());

  CHECK_UNARY_FALSE(Edge::te() == Edge::bs());
  CHECK_UNARY_FALSE(Edge::te() == Edge::tn());
  CHECK_UNARY_FALSE(Edge::te() == Edge::bn());
  CHECK_UNARY_FALSE(Edge::te() == Edge::ts());
  CHECK_UNARY_FALSE(Edge::te() == Edge::bw());
  CHECK_EQ(Edge::te(), Edge::te());
  CHECK_UNARY_FALSE(Edge::te() == Edge::be());
  CHECK_UNARY_FALSE(Edge::te() == Edge::tw());
  CHECK_UNARY_FALSE(Edge::te() == Edge::sw());
  CHECK_UNARY_FALSE(Edge::te() == Edge::ne());
  CHECK_UNARY_FALSE(Edge::te() == Edge::se());
  CHECK_UNARY_FALSE(Edge::te() == Edge::nw());
  CHECK_UNARY_FALSE(Edge::te() == Edge::null());

  CHECK_UNARY_FALSE(Edge::be() == Edge::bs());
  CHECK_UNARY_FALSE(Edge::be() == Edge::tn());
  CHECK_UNARY_FALSE(Edge::be() == Edge::bn());
  CHECK_UNARY_FALSE(Edge::be() == Edge::ts());
  CHECK_UNARY_FALSE(Edge::be() == Edge::bw());
  CHECK_UNARY_FALSE(Edge::be() == Edge::te());
  CHECK_EQ(Edge::be(), Edge::be());
  CHECK_UNARY_FALSE(Edge::be() == Edge::tw());
  CHECK_UNARY_FALSE(Edge::be() == Edge::sw());
  CHECK_UNARY_FALSE(Edge::be() == Edge::ne());
  CHECK_UNARY_FALSE(Edge::be() == Edge::se());
  CHECK_UNARY_FALSE(Edge::be() == Edge::nw());
  CHECK_UNARY_FALSE(Edge::be() == Edge::null());

  CHECK_UNARY_FALSE(Edge::tw() == Edge::bs());
  CHECK_UNARY_FALSE(Edge::tw() == Edge::tn());
  CHECK_UNARY_FALSE(Edge::tw() == Edge::bn());
  CHECK_UNARY_FALSE(Edge::tw() == Edge::ts());
  CHECK_UNARY_FALSE(Edge::tw() == Edge::bw());
  CHECK_UNARY_FALSE(Edge::tw() == Edge::te());
  CHECK_UNARY_FALSE(Edge::tw() == Edge::be());
  CHECK_EQ(Edge::tw(), Edge::tw());
  CHECK_UNARY_FALSE(Edge::tw() == Edge::sw());
  CHECK_UNARY_FALSE(Edge::tw() == Edge::ne());
  CHECK_UNARY_FALSE(Edge::tw() == Edge::se());
  CHECK_UNARY_FALSE(Edge::tw() == Edge::nw());
  CHECK_UNARY_FALSE(Edge::tw() == Edge::null());

  CHECK_UNARY_FALSE(Edge::sw() == Edge::bs());
  CHECK_UNARY_FALSE(Edge::sw() == Edge::tn());
  CHECK_UNARY_FALSE(Edge::sw() == Edge::bn());
  CHECK_UNARY_FALSE(Edge::sw() == Edge::ts());
  CHECK_UNARY_FALSE(Edge::sw() == Edge::bw());
  CHECK_UNARY_FALSE(Edge::sw() == Edge::te());
  CHECK_UNARY_FALSE(Edge::sw() == Edge::be());
  CHECK_UNARY_FALSE(Edge::sw() == Edge::tw());
  CHECK_EQ(Edge::sw(), Edge::sw());
  CHECK_UNARY_FALSE(Edge::sw() == Edge::ne());
  CHECK_UNARY_FALSE(Edge::sw() == Edge::se());
  CHECK_UNARY_FALSE(Edge::sw() == Edge::nw());
  CHECK_UNARY_FALSE(Edge::sw() == Edge::null());

  CHECK_UNARY_FALSE(Edge::ne() == Edge::bs());
  CHECK_UNARY_FALSE(Edge::ne() == Edge::tn());
  CHECK_UNARY_FALSE(Edge::ne() == Edge::bn());
  CHECK_UNARY_FALSE(Edge::ne() == Edge::ts());
  CHECK_UNARY_FALSE(Edge::ne() == Edge::bw());
  CHECK_UNARY_FALSE(Edge::ne() == Edge::te());
  CHECK_UNARY_FALSE(Edge::ne() == Edge::be());
  CHECK_UNARY_FALSE(Edge::ne() == Edge::tw());
  CHECK_UNARY_FALSE(Edge::ne() == Edge::sw());
  CHECK_EQ(Edge::ne(), Edge::ne());
  CHECK_UNARY_FALSE(Edge::ne() == Edge::se());
  CHECK_UNARY_FALSE(Edge::ne() == Edge::nw());
  CHECK_UNARY_FALSE(Edge::ne() == Edge::null());

  CHECK_UNARY_FALSE(Edge::se() == Edge::bs());
  CHECK_UNARY_FALSE(Edge::se() == Edge::tn());
  CHECK_UNARY_FALSE(Edge::se() == Edge::bn());
  CHECK_UNARY_FALSE(Edge::se() == Edge::ts());
  CHECK_UNARY_FALSE(Edge::se() == Edge::bw());
  CHECK_UNARY_FALSE(Edge::se() == Edge::te());
  CHECK_UNARY_FALSE(Edge::se() == Edge::be());
  CHECK_UNARY_FALSE(Edge::se() == Edge::tw());
  CHECK_UNARY_FALSE(Edge::se() == Edge::sw());
  CHECK_UNARY_FALSE(Edge::se() == Edge::ne());
  CHECK_EQ(Edge::se(), Edge::se());
  CHECK_UNARY_FALSE(Edge::se() == Edge::nw());
  CHECK_UNARY_FALSE(Edge::se() == Edge::null());

  CHECK_UNARY_FALSE(Edge::nw() == Edge::bs());
  CHECK_UNARY_FALSE(Edge::nw() == Edge::tn());
  CHECK_UNARY_FALSE(Edge::nw() == Edge::bn());
  CHECK_UNARY_FALSE(Edge::nw() == Edge::ts());
  CHECK_UNARY_FALSE(Edge::nw() == Edge::bw());
  CHECK_UNARY_FALSE(Edge::nw() == Edge::te());
  CHECK_UNARY_FALSE(Edge::nw() == Edge::be());
  CHECK_UNARY_FALSE(Edge::nw() == Edge::tw());
  CHECK_UNARY_FALSE(Edge::nw() == Edge::sw());
  CHECK_UNARY_FALSE(Edge::nw() == Edge::ne());
  CHECK_UNARY_FALSE(Edge::nw() == Edge::se());
  CHECK_EQ(Edge::nw(), Edge::nw());
  CHECK_UNARY_FALSE(Edge::nw() == Edge::null());

  CHECK_UNARY_FALSE(Edge::null() == Edge::bs());
  CHECK_UNARY_FALSE(Edge::null() == Edge::tn());
  CHECK_UNARY_FALSE(Edge::null() == Edge::bn());
  CHECK_UNARY_FALSE(Edge::null() == Edge::ts());
  CHECK_UNARY_FALSE(Edge::null() == Edge::bw());
  CHECK_UNARY_FALSE(Edge::null() == Edge::te());
  CHECK_UNARY_FALSE(Edge::null() == Edge::be());
  CHECK_UNARY_FALSE(Edge::null() == Edge::tw());
  CHECK_UNARY_FALSE(Edge::null() == Edge::sw());
  CHECK_UNARY_FALSE(Edge::null() == Edge::ne());
  CHECK_UNARY_FALSE(Edge::null() == Edge::se());
  CHECK_UNARY_FALSE(Edge::null() == Edge::nw());
  CHECK_EQ(Edge::null(), Edge::null());
}
TEST_CASE("Edge !=")
{
  CHECK_UNARY_FALSE(Edge::bs() != Edge::bs());
  CHECK_NE(Edge::bs(), Edge::tn());
  CHECK_NE(Edge::bs(), Edge::bn());
  CHECK_NE(Edge::bs(), Edge::ts());
  CHECK_NE(Edge::bs(), Edge::bw());
  CHECK_NE(Edge::bs(), Edge::te());
  CHECK_NE(Edge::bs(), Edge::be());
  CHECK_NE(Edge::bs(), Edge::tw());
  CHECK_NE(Edge::bs(), Edge::sw());
  CHECK_NE(Edge::bs(), Edge::ne());
  CHECK_NE(Edge::bs(), Edge::se());
  CHECK_NE(Edge::bs(), Edge::nw());
  CHECK_NE(Edge::bs(), Edge::null());

  CHECK_NE(Edge::tn(), Edge::bs());
  CHECK_UNARY_FALSE(Edge::tn() != Edge::tn());
  CHECK_NE(Edge::tn(), Edge::bn());
  CHECK_NE(Edge::tn(), Edge::ts());
  CHECK_NE(Edge::tn(), Edge::bw());
  CHECK_NE(Edge::tn(), Edge::te());
  CHECK_NE(Edge::tn(), Edge::be());
  CHECK_NE(Edge::tn(), Edge::tw());
  CHECK_NE(Edge::tn(), Edge::sw());
  CHECK_NE(Edge::tn(), Edge::ne());
  CHECK_NE(Edge::tn(), Edge::se());
  CHECK_NE(Edge::tn(), Edge::nw());
  CHECK_NE(Edge::tn(), Edge::null());

  CHECK_NE(Edge::bn(), Edge::bs());
  CHECK_NE(Edge::bn(), Edge::tn());
  CHECK_UNARY_FALSE(Edge::bn() != Edge::bn());
  CHECK_NE(Edge::bn(), Edge::ts());
  CHECK_NE(Edge::bn(), Edge::bw());
  CHECK_NE(Edge::bn(), Edge::te());
  CHECK_NE(Edge::bn(), Edge::be());
  CHECK_NE(Edge::bn(), Edge::tw());
  CHECK_NE(Edge::bn(), Edge::sw());
  CHECK_NE(Edge::bn(), Edge::ne());
  CHECK_NE(Edge::bn(), Edge::se());
  CHECK_NE(Edge::bn(), Edge::nw());
  CHECK_NE(Edge::bn(), Edge::null());

  CHECK_NE(Edge::ts(), Edge::bs());
  CHECK_NE(Edge::ts(), Edge::tn());
  CHECK_NE(Edge::ts(), Edge::bn());
  CHECK_UNARY_FALSE(Edge::ts() != Edge::ts());
  CHECK_NE(Edge::ts(), Edge::bw());
  CHECK_NE(Edge::ts(), Edge::te());
  CHECK_NE(Edge::ts(), Edge::be());
  CHECK_NE(Edge::ts(), Edge::tw());
  CHECK_NE(Edge::ts(), Edge::sw());
  CHECK_NE(Edge::ts(), Edge::ne());
  CHECK_NE(Edge::ts(), Edge::se());
  CHECK_NE(Edge::ts(), Edge::nw());
  CHECK_NE(Edge::ts(), Edge::null());

  CHECK_NE(Edge::bw(), Edge::bs());
  CHECK_NE(Edge::bw(), Edge::tn());
  CHECK_NE(Edge::bw(), Edge::bn());
  CHECK_NE(Edge::bw(), Edge::ts());
  CHECK_UNARY_FALSE(Edge::bw() != Edge::bw());
  CHECK_NE(Edge::bw(), Edge::te());
  CHECK_NE(Edge::bw(), Edge::be());
  CHECK_NE(Edge::bw(), Edge::tw());
  CHECK_NE(Edge::bw(), Edge::sw());
  CHECK_NE(Edge::bw(), Edge::ne());
  CHECK_NE(Edge::bw(), Edge::se());
  CHECK_NE(Edge::bw(), Edge::nw());
  CHECK_NE(Edge::bw(), Edge::null());

  CHECK_NE(Edge::te(), Edge::bs());
  CHECK_NE(Edge::te(), Edge::tn());
  CHECK_NE(Edge::te(), Edge::bn());
  CHECK_NE(Edge::te(), Edge::ts());
  CHECK_NE(Edge::te(), Edge::bw());
  CHECK_UNARY_FALSE(Edge::te() != Edge::te());
  CHECK_NE(Edge::te(), Edge::be());
  CHECK_NE(Edge::te(), Edge::tw());
  CHECK_NE(Edge::te(), Edge::sw());
  CHECK_NE(Edge::te(), Edge::ne());
  CHECK_NE(Edge::te(), Edge::se());
  CHECK_NE(Edge::te(), Edge::nw());
  CHECK_NE(Edge::te(), Edge::null());

  CHECK_NE(Edge::be(), Edge::bs());
  CHECK_NE(Edge::be(), Edge::tn());
  CHECK_NE(Edge::be(), Edge::bn());
  CHECK_NE(Edge::be(), Edge::ts());
  CHECK_NE(Edge::be(), Edge::bw());
  CHECK_NE(Edge::be(), Edge::te());
  CHECK_UNARY_FALSE(Edge::be() != Edge::be());
  CHECK_NE(Edge::be(), Edge::tw());
  CHECK_NE(Edge::be(), Edge::sw());
  CHECK_NE(Edge::be(), Edge::ne());
  CHECK_NE(Edge::be(), Edge::se());
  CHECK_NE(Edge::be(), Edge::nw());
  CHECK_NE(Edge::be(), Edge::null());

  CHECK_NE(Edge::tw(), Edge::bs());
  CHECK_NE(Edge::tw(), Edge::tn());
  CHECK_NE(Edge::tw(), Edge::bn());
  CHECK_NE(Edge::tw(), Edge::ts());
  CHECK_NE(Edge::tw(), Edge::bw());
  CHECK_NE(Edge::tw(), Edge::te());
  CHECK_NE(Edge::tw(), Edge::be());
  CHECK_UNARY_FALSE(Edge::tw() != Edge::tw());
  CHECK_NE(Edge::tw(), Edge::sw());
  CHECK_NE(Edge::tw(), Edge::ne());
  CHECK_NE(Edge::tw(), Edge::se());
  CHECK_NE(Edge::tw(), Edge::nw());
  CHECK_NE(Edge::tw(), Edge::null());

  CHECK_NE(Edge::sw(), Edge::bs());
  CHECK_NE(Edge::sw(), Edge::tn());
  CHECK_NE(Edge::sw(), Edge::bn());
  CHECK_NE(Edge::sw(), Edge::ts());
  CHECK_NE(Edge::sw(), Edge::bw());
  CHECK_NE(Edge::sw(), Edge::te());
  CHECK_NE(Edge::sw(), Edge::be());
  CHECK_NE(Edge::sw(), Edge::tw());
  CHECK_UNARY_FALSE(Edge::sw() != Edge::sw());
  CHECK_NE(Edge::sw(), Edge::ne());
  CHECK_NE(Edge::sw(), Edge::se());
  CHECK_NE(Edge::sw(), Edge::nw());
  CHECK_NE(Edge::sw(), Edge::null());

  CHECK_NE(Edge::ne(), Edge::bs());
  CHECK_NE(Edge::ne(), Edge::tn());
  CHECK_NE(Edge::ne(), Edge::bn());
  CHECK_NE(Edge::ne(), Edge::ts());
  CHECK_NE(Edge::ne(), Edge::bw());
  CHECK_NE(Edge::ne(), Edge::te());
  CHECK_NE(Edge::ne(), Edge::be());
  CHECK_NE(Edge::ne(), Edge::tw());
  CHECK_NE(Edge::ne(), Edge::sw());
  CHECK_UNARY_FALSE(Edge::ne() != Edge::ne());
  CHECK_NE(Edge::ne(), Edge::se());
  CHECK_NE(Edge::ne(), Edge::nw());
  CHECK_NE(Edge::ne(), Edge::null());

  CHECK_NE(Edge::se(), Edge::bs());
  CHECK_NE(Edge::se(), Edge::tn());
  CHECK_NE(Edge::se(), Edge::bn());
  CHECK_NE(Edge::se(), Edge::ts());
  CHECK_NE(Edge::se(), Edge::bw());
  CHECK_NE(Edge::se(), Edge::te());
  CHECK_NE(Edge::se(), Edge::be());
  CHECK_NE(Edge::se(), Edge::tw());
  CHECK_NE(Edge::se(), Edge::sw());
  CHECK_NE(Edge::se(), Edge::ne());
  CHECK_UNARY_FALSE(Edge::se() != Edge::se());
  CHECK_NE(Edge::se(), Edge::nw());
  CHECK_NE(Edge::se(), Edge::null());

  CHECK_NE(Edge::nw(), Edge::bs());
  CHECK_NE(Edge::nw(), Edge::tn());
  CHECK_NE(Edge::nw(), Edge::bn());
  CHECK_NE(Edge::nw(), Edge::ts());
  CHECK_NE(Edge::nw(), Edge::bw());
  CHECK_NE(Edge::nw(), Edge::te());
  CHECK_NE(Edge::nw(), Edge::be());
  CHECK_NE(Edge::nw(), Edge::tw());
  CHECK_NE(Edge::nw(), Edge::sw());
  CHECK_NE(Edge::nw(), Edge::ne());
  CHECK_NE(Edge::nw(), Edge::se());
  CHECK_UNARY_FALSE(Edge::nw() != Edge::nw());
  CHECK_NE(Edge::nw(), Edge::null());

  CHECK_NE(Edge::null(), Edge::bs());
  CHECK_NE(Edge::null(), Edge::tn());
  CHECK_NE(Edge::null(), Edge::bn());
  CHECK_NE(Edge::null(), Edge::ts());
  CHECK_NE(Edge::null(), Edge::bw());
  CHECK_NE(Edge::null(), Edge::te());
  CHECK_NE(Edge::null(), Edge::be());
  CHECK_NE(Edge::null(), Edge::tw());
  CHECK_NE(Edge::null(), Edge::sw());
  CHECK_NE(Edge::null(), Edge::ne());
  CHECK_NE(Edge::null(), Edge::se());
  CHECK_NE(Edge::null(), Edge::nw());
  CHECK_UNARY_FALSE(Edge::null() != Edge::null());
}
TEST_CASE("Edge <")
{
  for (unsigned char val = 0; val <= Edge::number_of; val++) {
    Edge edge(val);
    for (unsigned char other_val = 0; other_val <= Edge::number_of; other_val++) {
      Edge other_edge(other_val);
      CHECK_EQ((edge < other_edge), (val < other_val));
    }
  }
}
TEST_CASE("Test ostream for Edge")
{
  stringstream ss;
  ss << Edge::bs();
  CHECK_EQ(ss.str(), "Edge::bs()");
  ss.str("");
  ss << Edge::tn();
  CHECK_EQ(ss.str(), "Edge::tn()");
  ss.str("");
  ss << Edge::bn();
  CHECK_EQ(ss.str(), "Edge::bn()");
  ss.str("");
  ss << Edge::ts();
  CHECK_EQ(ss.str(), "Edge::ts()");
  ss.str("");
  ss << Edge::bw();
  CHECK_EQ(ss.str(), "Edge::bw()");
  ss.str("");
  ss << Edge::te();
  CHECK_EQ(ss.str(), "Edge::te()");
  ss.str("");
  ss << Edge::be();
  CHECK_EQ(ss.str(), "Edge::be()");
  ss.str("");
  ss << Edge::tw();
  CHECK_EQ(ss.str(), "Edge::tw()");
  ss.str("");
  ss << Edge::sw();
  CHECK_EQ(ss.str(), "Edge::sw()");
  ss.str("");
  ss << Edge::ne();
  CHECK_EQ(ss.str(), "Edge::ne()");
  ss.str("");
  ss << Edge::se();
  CHECK_EQ(ss.str(), "Edge::se()");
  ss.str("");
  ss << Edge::nw();
  CHECK_EQ(ss.str(), "Edge::nw()");
  ss.str("");
  ss << Edge::null();
  CHECK_EQ(ss.str(), "Edge::null()");
  ss.str("");
  ss << Edge(64);
  CHECK_EQ(ss.str(), "Edge invalid value: 64");
}
TEST_CASE("Test iterator for Edge")
{
  auto iter = Edge::getValues().begin();
  CHECK_EQ(iter, Edge::getValues().begin());
  CHECK_NE(iter, Edge::getValues().end());
  CHECK_EQ(*iter, Edge::bs());
  ++iter;
  CHECK_EQ(iter->getIndex(), 1);
  CHECK_EQ(*iter, Edge::bn());
  ++iter;
  CHECK_EQ(*iter, Edge::ts());
  ++iter;
  CHECK_EQ(*iter, Edge::tn());
  ++iter;
  CHECK_EQ(*iter, Edge::bw());
  ++iter;
  CHECK_EQ(*iter, Edge::be());
  ++iter;
  CHECK_EQ(*iter, Edge::tw());
  ++iter;
  CHECK_EQ(*iter, Edge::te());
  ++iter;
  CHECK_EQ(*iter, Edge::sw());
  ++iter;
  CHECK_EQ(*iter, Edge::se());
  ++iter;
  CHECK_EQ(*iter, Edge::nw());
  ++iter;
  CHECK_EQ(*iter, Edge::ne());
  ++iter;
  CHECK_EQ(*iter, Edge::null());
  CHECK_EQ(iter, Edge::getValues().end());
}
