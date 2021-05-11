#include <ThunderEgg/Edge.h>

#include <sstream>

#include <catch2/catch_test_macros.hpp>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Edge<1> unsigned char constructor works", "[Edge]")
{
	Edge<1> o(13);
	CHECK(o.getIndex() == 13);
}
TEST_CASE("Edge<2> unsigned char constructor works", "[Edge]")
{
	Edge<2> o(13);
	CHECK(o.getIndex() == 13);
}
TEST_CASE("Edge<3> unsigned char constructor works", "[Edge]")
{
	Edge<3> o(13);
	CHECK(o.getIndex() == 13);
}
TEST_CASE("Edge<1> Default constructor works", "[Edge]")
{
	Edge<1> o;
	CHECK(o == Edge<1>::null());
}
TEST_CASE("Edge<2> Default constructor works", "[Edge]")
{
	Edge<2> o;
	CHECK(o == Edge<2>::null());
}
TEST_CASE("Edge<3> Default constructor works", "[Edge]")
{
	Edge<3> o;
	CHECK(o == Edge<3>::null());
}
TEST_CASE("Edge<1> named constructors give expected index values", "[Edge]")
{
	CHECK(Edge<1>::null().getIndex() == 0);
}
TEST_CASE("Edge<2> named constructors give expected index values", "[Edge]")
{
	CHECK(Edge<2>::null().getIndex() == 0);
}
TEST_CASE("Edge<3> named constructors give expected index values", "[Edge]")
{
	CHECK(Edge<3>::bs().getIndex() == 0);
	CHECK(Edge<3>::tn().getIndex() == 1);
	CHECK(Edge<3>::bn().getIndex() == 2);
	CHECK(Edge<3>::ts().getIndex() == 3);
	CHECK(Edge<3>::bw().getIndex() == 4);
	CHECK(Edge<3>::te().getIndex() == 5);
	CHECK(Edge<3>::be().getIndex() == 6);
	CHECK(Edge<3>::tw().getIndex() == 7);
	CHECK(Edge<3>::sw().getIndex() == 8);
	CHECK(Edge<3>::ne().getIndex() == 9);
	CHECK(Edge<3>::se().getIndex() == 10);
	CHECK(Edge<3>::nw().getIndex() == 11);
	CHECK(Edge<3>::null().getIndex() == 12);
}
TEST_CASE("Edge<3> opposite", "[Edge]")
{
	CHECK(Edge<3>::bs().opposite() == Edge<3>::tn());
	CHECK(Edge<3>::tn().opposite() == Edge<3>::bs());
	CHECK(Edge<3>::bn().opposite() == Edge<3>::ts());
	CHECK(Edge<3>::ts().opposite() == Edge<3>::bn());
	CHECK(Edge<3>::bw().opposite() == Edge<3>::te());
	CHECK(Edge<3>::te().opposite() == Edge<3>::bw());
	CHECK(Edge<3>::be().opposite() == Edge<3>::tw());
	CHECK(Edge<3>::tw().opposite() == Edge<3>::be());
	CHECK(Edge<3>::sw().opposite() == Edge<3>::ne());
	CHECK(Edge<3>::ne().opposite() == Edge<3>::sw());
	CHECK(Edge<3>::se().opposite() == Edge<3>::nw());
	CHECK(Edge<3>::nw().opposite() == Edge<3>::se());
}
TEST_CASE("Edge<3> getAxisIndex", "[Edge]")
{
	CHECK(Edge<3>::bs().getAxisIndex() == 0);
	CHECK(Edge<3>::tn().getAxisIndex() == 0);
	CHECK(Edge<3>::bn().getAxisIndex() == 0);
	CHECK(Edge<3>::ts().getAxisIndex() == 0);
	CHECK(Edge<3>::bw().getAxisIndex() == 1);
	CHECK(Edge<3>::te().getAxisIndex() == 1);
	CHECK(Edge<3>::be().getAxisIndex() == 1);
	CHECK(Edge<3>::tw().getAxisIndex() == 1);
	CHECK(Edge<3>::sw().getAxisIndex() == 2);
	CHECK(Edge<3>::ne().getAxisIndex() == 2);
	CHECK(Edge<3>::se().getAxisIndex() == 2);
	CHECK(Edge<3>::nw().getAxisIndex() == 2);
}
TEST_CASE("Edge<3> getSides", "[Edge]")
{
	CHECK(Edge<3>::bs().getSides()[0] == Side<3>::south());
	CHECK(Edge<3>::bs().getSides()[1] == Side<3>::bottom());

	CHECK(Edge<3>::tn().getSides()[0] == Side<3>::north());
	CHECK(Edge<3>::tn().getSides()[1] == Side<3>::top());

	CHECK(Edge<3>::bn().getSides()[0] == Side<3>::north());
	CHECK(Edge<3>::bn().getSides()[1] == Side<3>::bottom());

	CHECK(Edge<3>::ts().getSides()[0] == Side<3>::south());
	CHECK(Edge<3>::ts().getSides()[1] == Side<3>::top());

	CHECK(Edge<3>::bw().getSides()[0] == Side<3>::west());
	CHECK(Edge<3>::bw().getSides()[1] == Side<3>::bottom());

	CHECK(Edge<3>::te().getSides()[0] == Side<3>::east());
	CHECK(Edge<3>::te().getSides()[1] == Side<3>::top());

	CHECK(Edge<3>::be().getSides()[0] == Side<3>::east());
	CHECK(Edge<3>::be().getSides()[1] == Side<3>::bottom());

	CHECK(Edge<3>::tw().getSides()[0] == Side<3>::west());
	CHECK(Edge<3>::tw().getSides()[1] == Side<3>::top());

	CHECK(Edge<3>::sw().getSides()[0] == Side<3>::west());
	CHECK(Edge<3>::sw().getSides()[1] == Side<3>::south());

	CHECK(Edge<3>::ne().getSides()[0] == Side<3>::east());
	CHECK(Edge<3>::ne().getSides()[1] == Side<3>::north());

	CHECK(Edge<3>::se().getSides()[0] == Side<3>::east());
	CHECK(Edge<3>::se().getSides()[1] == Side<3>::south());

	CHECK(Edge<3>::nw().getSides()[0] == Side<3>::west());
	CHECK(Edge<3>::nw().getSides()[1] == Side<3>::north());
}
TEST_CASE("Edge<1> ==", "[Edge]")
{
	CHECK(Edge<1>::null() == Edge<1>::null());
}
TEST_CASE("Edge<2> ==", "[Edge]")
{
	CHECK(Edge<2>::null() == Edge<2>::null());
}
TEST_CASE("Edge<3> ==", "[Edge]")
{
	CHECK(Edge<3>::bs() == Edge<3>::bs());
	CHECK_FALSE(Edge<3>::bs() == Edge<3>::tn());
	CHECK_FALSE(Edge<3>::bs() == Edge<3>::bn());
	CHECK_FALSE(Edge<3>::bs() == Edge<3>::ts());
	CHECK_FALSE(Edge<3>::bs() == Edge<3>::bw());
	CHECK_FALSE(Edge<3>::bs() == Edge<3>::te());
	CHECK_FALSE(Edge<3>::bs() == Edge<3>::be());
	CHECK_FALSE(Edge<3>::bs() == Edge<3>::tw());
	CHECK_FALSE(Edge<3>::bs() == Edge<3>::sw());
	CHECK_FALSE(Edge<3>::bs() == Edge<3>::ne());
	CHECK_FALSE(Edge<3>::bs() == Edge<3>::se());
	CHECK_FALSE(Edge<3>::bs() == Edge<3>::nw());
	CHECK_FALSE(Edge<3>::bs() == Edge<3>::null());

	CHECK_FALSE(Edge<3>::tn() == Edge<3>::bs());
	CHECK(Edge<3>::tn() == Edge<3>::tn());
	CHECK_FALSE(Edge<3>::tn() == Edge<3>::bn());
	CHECK_FALSE(Edge<3>::tn() == Edge<3>::ts());
	CHECK_FALSE(Edge<3>::tn() == Edge<3>::bw());
	CHECK_FALSE(Edge<3>::tn() == Edge<3>::te());
	CHECK_FALSE(Edge<3>::tn() == Edge<3>::be());
	CHECK_FALSE(Edge<3>::tn() == Edge<3>::tw());
	CHECK_FALSE(Edge<3>::tn() == Edge<3>::sw());
	CHECK_FALSE(Edge<3>::tn() == Edge<3>::ne());
	CHECK_FALSE(Edge<3>::tn() == Edge<3>::se());
	CHECK_FALSE(Edge<3>::tn() == Edge<3>::nw());
	CHECK_FALSE(Edge<3>::tn() == Edge<3>::null());

	CHECK_FALSE(Edge<3>::bn() == Edge<3>::bs());
	CHECK_FALSE(Edge<3>::bn() == Edge<3>::tn());
	CHECK(Edge<3>::bn() == Edge<3>::bn());
	CHECK_FALSE(Edge<3>::bn() == Edge<3>::ts());
	CHECK_FALSE(Edge<3>::bn() == Edge<3>::bw());
	CHECK_FALSE(Edge<3>::bn() == Edge<3>::te());
	CHECK_FALSE(Edge<3>::bn() == Edge<3>::be());
	CHECK_FALSE(Edge<3>::bn() == Edge<3>::tw());
	CHECK_FALSE(Edge<3>::bn() == Edge<3>::sw());
	CHECK_FALSE(Edge<3>::bn() == Edge<3>::ne());
	CHECK_FALSE(Edge<3>::bn() == Edge<3>::se());
	CHECK_FALSE(Edge<3>::bn() == Edge<3>::nw());
	CHECK_FALSE(Edge<3>::bn() == Edge<3>::null());

	CHECK_FALSE(Edge<3>::ts() == Edge<3>::bs());
	CHECK_FALSE(Edge<3>::ts() == Edge<3>::tn());
	CHECK_FALSE(Edge<3>::ts() == Edge<3>::bn());
	CHECK(Edge<3>::ts() == Edge<3>::ts());
	CHECK_FALSE(Edge<3>::ts() == Edge<3>::bw());
	CHECK_FALSE(Edge<3>::ts() == Edge<3>::te());
	CHECK_FALSE(Edge<3>::ts() == Edge<3>::be());
	CHECK_FALSE(Edge<3>::ts() == Edge<3>::tw());
	CHECK_FALSE(Edge<3>::ts() == Edge<3>::sw());
	CHECK_FALSE(Edge<3>::ts() == Edge<3>::ne());
	CHECK_FALSE(Edge<3>::ts() == Edge<3>::se());
	CHECK_FALSE(Edge<3>::ts() == Edge<3>::nw());
	CHECK_FALSE(Edge<3>::ts() == Edge<3>::null());

	CHECK_FALSE(Edge<3>::bw() == Edge<3>::bs());
	CHECK_FALSE(Edge<3>::bw() == Edge<3>::tn());
	CHECK_FALSE(Edge<3>::bw() == Edge<3>::bn());
	CHECK_FALSE(Edge<3>::bw() == Edge<3>::ts());
	CHECK(Edge<3>::bw() == Edge<3>::bw());
	CHECK_FALSE(Edge<3>::bw() == Edge<3>::te());
	CHECK_FALSE(Edge<3>::bw() == Edge<3>::be());
	CHECK_FALSE(Edge<3>::bw() == Edge<3>::tw());
	CHECK_FALSE(Edge<3>::bw() == Edge<3>::sw());
	CHECK_FALSE(Edge<3>::bw() == Edge<3>::ne());
	CHECK_FALSE(Edge<3>::bw() == Edge<3>::se());
	CHECK_FALSE(Edge<3>::bw() == Edge<3>::nw());
	CHECK_FALSE(Edge<3>::bw() == Edge<3>::null());

	CHECK_FALSE(Edge<3>::te() == Edge<3>::bs());
	CHECK_FALSE(Edge<3>::te() == Edge<3>::tn());
	CHECK_FALSE(Edge<3>::te() == Edge<3>::bn());
	CHECK_FALSE(Edge<3>::te() == Edge<3>::ts());
	CHECK_FALSE(Edge<3>::te() == Edge<3>::bw());
	CHECK(Edge<3>::te() == Edge<3>::te());
	CHECK_FALSE(Edge<3>::te() == Edge<3>::be());
	CHECK_FALSE(Edge<3>::te() == Edge<3>::tw());
	CHECK_FALSE(Edge<3>::te() == Edge<3>::sw());
	CHECK_FALSE(Edge<3>::te() == Edge<3>::ne());
	CHECK_FALSE(Edge<3>::te() == Edge<3>::se());
	CHECK_FALSE(Edge<3>::te() == Edge<3>::nw());
	CHECK_FALSE(Edge<3>::te() == Edge<3>::null());

	CHECK_FALSE(Edge<3>::be() == Edge<3>::bs());
	CHECK_FALSE(Edge<3>::be() == Edge<3>::tn());
	CHECK_FALSE(Edge<3>::be() == Edge<3>::bn());
	CHECK_FALSE(Edge<3>::be() == Edge<3>::ts());
	CHECK_FALSE(Edge<3>::be() == Edge<3>::bw());
	CHECK_FALSE(Edge<3>::be() == Edge<3>::te());
	CHECK(Edge<3>::be() == Edge<3>::be());
	CHECK_FALSE(Edge<3>::be() == Edge<3>::tw());
	CHECK_FALSE(Edge<3>::be() == Edge<3>::sw());
	CHECK_FALSE(Edge<3>::be() == Edge<3>::ne());
	CHECK_FALSE(Edge<3>::be() == Edge<3>::se());
	CHECK_FALSE(Edge<3>::be() == Edge<3>::nw());
	CHECK_FALSE(Edge<3>::be() == Edge<3>::null());

	CHECK_FALSE(Edge<3>::tw() == Edge<3>::bs());
	CHECK_FALSE(Edge<3>::tw() == Edge<3>::tn());
	CHECK_FALSE(Edge<3>::tw() == Edge<3>::bn());
	CHECK_FALSE(Edge<3>::tw() == Edge<3>::ts());
	CHECK_FALSE(Edge<3>::tw() == Edge<3>::bw());
	CHECK_FALSE(Edge<3>::tw() == Edge<3>::te());
	CHECK_FALSE(Edge<3>::tw() == Edge<3>::be());
	CHECK(Edge<3>::tw() == Edge<3>::tw());
	CHECK_FALSE(Edge<3>::tw() == Edge<3>::sw());
	CHECK_FALSE(Edge<3>::tw() == Edge<3>::ne());
	CHECK_FALSE(Edge<3>::tw() == Edge<3>::se());
	CHECK_FALSE(Edge<3>::tw() == Edge<3>::nw());
	CHECK_FALSE(Edge<3>::tw() == Edge<3>::null());

	CHECK_FALSE(Edge<3>::sw() == Edge<3>::bs());
	CHECK_FALSE(Edge<3>::sw() == Edge<3>::tn());
	CHECK_FALSE(Edge<3>::sw() == Edge<3>::bn());
	CHECK_FALSE(Edge<3>::sw() == Edge<3>::ts());
	CHECK_FALSE(Edge<3>::sw() == Edge<3>::bw());
	CHECK_FALSE(Edge<3>::sw() == Edge<3>::te());
	CHECK_FALSE(Edge<3>::sw() == Edge<3>::be());
	CHECK_FALSE(Edge<3>::sw() == Edge<3>::tw());
	CHECK(Edge<3>::sw() == Edge<3>::sw());
	CHECK_FALSE(Edge<3>::sw() == Edge<3>::ne());
	CHECK_FALSE(Edge<3>::sw() == Edge<3>::se());
	CHECK_FALSE(Edge<3>::sw() == Edge<3>::nw());
	CHECK_FALSE(Edge<3>::sw() == Edge<3>::null());

	CHECK_FALSE(Edge<3>::ne() == Edge<3>::bs());
	CHECK_FALSE(Edge<3>::ne() == Edge<3>::tn());
	CHECK_FALSE(Edge<3>::ne() == Edge<3>::bn());
	CHECK_FALSE(Edge<3>::ne() == Edge<3>::ts());
	CHECK_FALSE(Edge<3>::ne() == Edge<3>::bw());
	CHECK_FALSE(Edge<3>::ne() == Edge<3>::te());
	CHECK_FALSE(Edge<3>::ne() == Edge<3>::be());
	CHECK_FALSE(Edge<3>::ne() == Edge<3>::tw());
	CHECK_FALSE(Edge<3>::ne() == Edge<3>::sw());
	CHECK(Edge<3>::ne() == Edge<3>::ne());
	CHECK_FALSE(Edge<3>::ne() == Edge<3>::se());
	CHECK_FALSE(Edge<3>::ne() == Edge<3>::nw());
	CHECK_FALSE(Edge<3>::ne() == Edge<3>::null());

	CHECK_FALSE(Edge<3>::se() == Edge<3>::bs());
	CHECK_FALSE(Edge<3>::se() == Edge<3>::tn());
	CHECK_FALSE(Edge<3>::se() == Edge<3>::bn());
	CHECK_FALSE(Edge<3>::se() == Edge<3>::ts());
	CHECK_FALSE(Edge<3>::se() == Edge<3>::bw());
	CHECK_FALSE(Edge<3>::se() == Edge<3>::te());
	CHECK_FALSE(Edge<3>::se() == Edge<3>::be());
	CHECK_FALSE(Edge<3>::se() == Edge<3>::tw());
	CHECK_FALSE(Edge<3>::se() == Edge<3>::sw());
	CHECK_FALSE(Edge<3>::se() == Edge<3>::ne());
	CHECK(Edge<3>::se() == Edge<3>::se());
	CHECK_FALSE(Edge<3>::se() == Edge<3>::nw());
	CHECK_FALSE(Edge<3>::se() == Edge<3>::null());

	CHECK_FALSE(Edge<3>::nw() == Edge<3>::bs());
	CHECK_FALSE(Edge<3>::nw() == Edge<3>::tn());
	CHECK_FALSE(Edge<3>::nw() == Edge<3>::bn());
	CHECK_FALSE(Edge<3>::nw() == Edge<3>::ts());
	CHECK_FALSE(Edge<3>::nw() == Edge<3>::bw());
	CHECK_FALSE(Edge<3>::nw() == Edge<3>::te());
	CHECK_FALSE(Edge<3>::nw() == Edge<3>::be());
	CHECK_FALSE(Edge<3>::nw() == Edge<3>::tw());
	CHECK_FALSE(Edge<3>::nw() == Edge<3>::sw());
	CHECK_FALSE(Edge<3>::nw() == Edge<3>::ne());
	CHECK_FALSE(Edge<3>::nw() == Edge<3>::se());
	CHECK(Edge<3>::nw() == Edge<3>::nw());
	CHECK_FALSE(Edge<3>::nw() == Edge<3>::null());

	CHECK_FALSE(Edge<3>::null() == Edge<3>::bs());
	CHECK_FALSE(Edge<3>::null() == Edge<3>::tn());
	CHECK_FALSE(Edge<3>::null() == Edge<3>::bn());
	CHECK_FALSE(Edge<3>::null() == Edge<3>::ts());
	CHECK_FALSE(Edge<3>::null() == Edge<3>::bw());
	CHECK_FALSE(Edge<3>::null() == Edge<3>::te());
	CHECK_FALSE(Edge<3>::null() == Edge<3>::be());
	CHECK_FALSE(Edge<3>::null() == Edge<3>::tw());
	CHECK_FALSE(Edge<3>::null() == Edge<3>::sw());
	CHECK_FALSE(Edge<3>::null() == Edge<3>::ne());
	CHECK_FALSE(Edge<3>::null() == Edge<3>::se());
	CHECK_FALSE(Edge<3>::null() == Edge<3>::nw());
	CHECK(Edge<3>::null() == Edge<3>::null());
}
TEST_CASE("Edge<1> !=", "[Edge]")
{
	CHECK_FALSE(Edge<1>::null() != Edge<1>::null());
}
TEST_CASE("Edge<2> !=", "[Edge]")
{
	CHECK_FALSE(Edge<2>::null() != Edge<2>::null());
}
TEST_CASE("Edge<3> !=", "[Edge]")
{
	CHECK_FALSE(Edge<3>::bs() != Edge<3>::bs());
	CHECK(Edge<3>::bs() != Edge<3>::tn());
	CHECK(Edge<3>::bs() != Edge<3>::bn());
	CHECK(Edge<3>::bs() != Edge<3>::ts());
	CHECK(Edge<3>::bs() != Edge<3>::bw());
	CHECK(Edge<3>::bs() != Edge<3>::te());
	CHECK(Edge<3>::bs() != Edge<3>::be());
	CHECK(Edge<3>::bs() != Edge<3>::tw());
	CHECK(Edge<3>::bs() != Edge<3>::sw());
	CHECK(Edge<3>::bs() != Edge<3>::ne());
	CHECK(Edge<3>::bs() != Edge<3>::se());
	CHECK(Edge<3>::bs() != Edge<3>::nw());
	CHECK(Edge<3>::bs() != Edge<3>::null());

	CHECK(Edge<3>::tn() != Edge<3>::bs());
	CHECK_FALSE(Edge<3>::tn() != Edge<3>::tn());
	CHECK(Edge<3>::tn() != Edge<3>::bn());
	CHECK(Edge<3>::tn() != Edge<3>::ts());
	CHECK(Edge<3>::tn() != Edge<3>::bw());
	CHECK(Edge<3>::tn() != Edge<3>::te());
	CHECK(Edge<3>::tn() != Edge<3>::be());
	CHECK(Edge<3>::tn() != Edge<3>::tw());
	CHECK(Edge<3>::tn() != Edge<3>::sw());
	CHECK(Edge<3>::tn() != Edge<3>::ne());
	CHECK(Edge<3>::tn() != Edge<3>::se());
	CHECK(Edge<3>::tn() != Edge<3>::nw());
	CHECK(Edge<3>::tn() != Edge<3>::null());

	CHECK(Edge<3>::bn() != Edge<3>::bs());
	CHECK(Edge<3>::bn() != Edge<3>::tn());
	CHECK_FALSE(Edge<3>::bn() != Edge<3>::bn());
	CHECK(Edge<3>::bn() != Edge<3>::ts());
	CHECK(Edge<3>::bn() != Edge<3>::bw());
	CHECK(Edge<3>::bn() != Edge<3>::te());
	CHECK(Edge<3>::bn() != Edge<3>::be());
	CHECK(Edge<3>::bn() != Edge<3>::tw());
	CHECK(Edge<3>::bn() != Edge<3>::sw());
	CHECK(Edge<3>::bn() != Edge<3>::ne());
	CHECK(Edge<3>::bn() != Edge<3>::se());
	CHECK(Edge<3>::bn() != Edge<3>::nw());
	CHECK(Edge<3>::bn() != Edge<3>::null());

	CHECK(Edge<3>::ts() != Edge<3>::bs());
	CHECK(Edge<3>::ts() != Edge<3>::tn());
	CHECK(Edge<3>::ts() != Edge<3>::bn());
	CHECK_FALSE(Edge<3>::ts() != Edge<3>::ts());
	CHECK(Edge<3>::ts() != Edge<3>::bw());
	CHECK(Edge<3>::ts() != Edge<3>::te());
	CHECK(Edge<3>::ts() != Edge<3>::be());
	CHECK(Edge<3>::ts() != Edge<3>::tw());
	CHECK(Edge<3>::ts() != Edge<3>::sw());
	CHECK(Edge<3>::ts() != Edge<3>::ne());
	CHECK(Edge<3>::ts() != Edge<3>::se());
	CHECK(Edge<3>::ts() != Edge<3>::nw());
	CHECK(Edge<3>::ts() != Edge<3>::null());

	CHECK(Edge<3>::bw() != Edge<3>::bs());
	CHECK(Edge<3>::bw() != Edge<3>::tn());
	CHECK(Edge<3>::bw() != Edge<3>::bn());
	CHECK(Edge<3>::bw() != Edge<3>::ts());
	CHECK_FALSE(Edge<3>::bw() != Edge<3>::bw());
	CHECK(Edge<3>::bw() != Edge<3>::te());
	CHECK(Edge<3>::bw() != Edge<3>::be());
	CHECK(Edge<3>::bw() != Edge<3>::tw());
	CHECK(Edge<3>::bw() != Edge<3>::sw());
	CHECK(Edge<3>::bw() != Edge<3>::ne());
	CHECK(Edge<3>::bw() != Edge<3>::se());
	CHECK(Edge<3>::bw() != Edge<3>::nw());
	CHECK(Edge<3>::bw() != Edge<3>::null());

	CHECK(Edge<3>::te() != Edge<3>::bs());
	CHECK(Edge<3>::te() != Edge<3>::tn());
	CHECK(Edge<3>::te() != Edge<3>::bn());
	CHECK(Edge<3>::te() != Edge<3>::ts());
	CHECK(Edge<3>::te() != Edge<3>::bw());
	CHECK_FALSE(Edge<3>::te() != Edge<3>::te());
	CHECK(Edge<3>::te() != Edge<3>::be());
	CHECK(Edge<3>::te() != Edge<3>::tw());
	CHECK(Edge<3>::te() != Edge<3>::sw());
	CHECK(Edge<3>::te() != Edge<3>::ne());
	CHECK(Edge<3>::te() != Edge<3>::se());
	CHECK(Edge<3>::te() != Edge<3>::nw());
	CHECK(Edge<3>::te() != Edge<3>::null());

	CHECK(Edge<3>::be() != Edge<3>::bs());
	CHECK(Edge<3>::be() != Edge<3>::tn());
	CHECK(Edge<3>::be() != Edge<3>::bn());
	CHECK(Edge<3>::be() != Edge<3>::ts());
	CHECK(Edge<3>::be() != Edge<3>::bw());
	CHECK(Edge<3>::be() != Edge<3>::te());
	CHECK_FALSE(Edge<3>::be() != Edge<3>::be());
	CHECK(Edge<3>::be() != Edge<3>::tw());
	CHECK(Edge<3>::be() != Edge<3>::sw());
	CHECK(Edge<3>::be() != Edge<3>::ne());
	CHECK(Edge<3>::be() != Edge<3>::se());
	CHECK(Edge<3>::be() != Edge<3>::nw());
	CHECK(Edge<3>::be() != Edge<3>::null());

	CHECK(Edge<3>::tw() != Edge<3>::bs());
	CHECK(Edge<3>::tw() != Edge<3>::tn());
	CHECK(Edge<3>::tw() != Edge<3>::bn());
	CHECK(Edge<3>::tw() != Edge<3>::ts());
	CHECK(Edge<3>::tw() != Edge<3>::bw());
	CHECK(Edge<3>::tw() != Edge<3>::te());
	CHECK(Edge<3>::tw() != Edge<3>::be());
	CHECK_FALSE(Edge<3>::tw() != Edge<3>::tw());
	CHECK(Edge<3>::tw() != Edge<3>::sw());
	CHECK(Edge<3>::tw() != Edge<3>::ne());
	CHECK(Edge<3>::tw() != Edge<3>::se());
	CHECK(Edge<3>::tw() != Edge<3>::nw());
	CHECK(Edge<3>::tw() != Edge<3>::null());

	CHECK(Edge<3>::sw() != Edge<3>::bs());
	CHECK(Edge<3>::sw() != Edge<3>::tn());
	CHECK(Edge<3>::sw() != Edge<3>::bn());
	CHECK(Edge<3>::sw() != Edge<3>::ts());
	CHECK(Edge<3>::sw() != Edge<3>::bw());
	CHECK(Edge<3>::sw() != Edge<3>::te());
	CHECK(Edge<3>::sw() != Edge<3>::be());
	CHECK(Edge<3>::sw() != Edge<3>::tw());
	CHECK_FALSE(Edge<3>::sw() != Edge<3>::sw());
	CHECK(Edge<3>::sw() != Edge<3>::ne());
	CHECK(Edge<3>::sw() != Edge<3>::se());
	CHECK(Edge<3>::sw() != Edge<3>::nw());
	CHECK(Edge<3>::sw() != Edge<3>::null());

	CHECK(Edge<3>::ne() != Edge<3>::bs());
	CHECK(Edge<3>::ne() != Edge<3>::tn());
	CHECK(Edge<3>::ne() != Edge<3>::bn());
	CHECK(Edge<3>::ne() != Edge<3>::ts());
	CHECK(Edge<3>::ne() != Edge<3>::bw());
	CHECK(Edge<3>::ne() != Edge<3>::te());
	CHECK(Edge<3>::ne() != Edge<3>::be());
	CHECK(Edge<3>::ne() != Edge<3>::tw());
	CHECK(Edge<3>::ne() != Edge<3>::sw());
	CHECK_FALSE(Edge<3>::ne() != Edge<3>::ne());
	CHECK(Edge<3>::ne() != Edge<3>::se());
	CHECK(Edge<3>::ne() != Edge<3>::nw());
	CHECK(Edge<3>::ne() != Edge<3>::null());

	CHECK(Edge<3>::se() != Edge<3>::bs());
	CHECK(Edge<3>::se() != Edge<3>::tn());
	CHECK(Edge<3>::se() != Edge<3>::bn());
	CHECK(Edge<3>::se() != Edge<3>::ts());
	CHECK(Edge<3>::se() != Edge<3>::bw());
	CHECK(Edge<3>::se() != Edge<3>::te());
	CHECK(Edge<3>::se() != Edge<3>::be());
	CHECK(Edge<3>::se() != Edge<3>::tw());
	CHECK(Edge<3>::se() != Edge<3>::sw());
	CHECK(Edge<3>::se() != Edge<3>::ne());
	CHECK_FALSE(Edge<3>::se() != Edge<3>::se());
	CHECK(Edge<3>::se() != Edge<3>::nw());
	CHECK(Edge<3>::se() != Edge<3>::null());

	CHECK(Edge<3>::nw() != Edge<3>::bs());
	CHECK(Edge<3>::nw() != Edge<3>::tn());
	CHECK(Edge<3>::nw() != Edge<3>::bn());
	CHECK(Edge<3>::nw() != Edge<3>::ts());
	CHECK(Edge<3>::nw() != Edge<3>::bw());
	CHECK(Edge<3>::nw() != Edge<3>::te());
	CHECK(Edge<3>::nw() != Edge<3>::be());
	CHECK(Edge<3>::nw() != Edge<3>::tw());
	CHECK(Edge<3>::nw() != Edge<3>::sw());
	CHECK(Edge<3>::nw() != Edge<3>::ne());
	CHECK(Edge<3>::nw() != Edge<3>::se());
	CHECK_FALSE(Edge<3>::nw() != Edge<3>::nw());
	CHECK(Edge<3>::nw() != Edge<3>::null());

	CHECK(Edge<3>::null() != Edge<3>::bs());
	CHECK(Edge<3>::null() != Edge<3>::tn());
	CHECK(Edge<3>::null() != Edge<3>::bn());
	CHECK(Edge<3>::null() != Edge<3>::ts());
	CHECK(Edge<3>::null() != Edge<3>::bw());
	CHECK(Edge<3>::null() != Edge<3>::te());
	CHECK(Edge<3>::null() != Edge<3>::be());
	CHECK(Edge<3>::null() != Edge<3>::tw());
	CHECK(Edge<3>::null() != Edge<3>::sw());
	CHECK(Edge<3>::null() != Edge<3>::ne());
	CHECK(Edge<3>::null() != Edge<3>::se());
	CHECK(Edge<3>::null() != Edge<3>::nw());
	CHECK_FALSE(Edge<3>::null() != Edge<3>::null());
}
TEST_CASE("Edge<1> <", "[Edge]")
{
	CHECK_FALSE(Edge<1>::null() < Edge<1>::null());
}
TEST_CASE("Edge<2> <", "[Edge]")
{
	CHECK_FALSE(Edge<2>::null() < Edge<2>::null());
}
TEST_CASE("Edge<3> <", "[Edge]")
{
	CHECK_FALSE(Edge<3>::bs() < Edge<3>::bs());
	CHECK(Edge<3>::bs() < Edge<3>::tn());
	CHECK(Edge<3>::bs() < Edge<3>::bn());
	CHECK(Edge<3>::bs() < Edge<3>::ts());
	CHECK(Edge<3>::bs() < Edge<3>::bw());
	CHECK(Edge<3>::bs() < Edge<3>::te());
	CHECK(Edge<3>::bs() < Edge<3>::be());
	CHECK(Edge<3>::bs() < Edge<3>::tw());
	CHECK(Edge<3>::bs() < Edge<3>::sw());
	CHECK(Edge<3>::bs() < Edge<3>::ne());
	CHECK(Edge<3>::bs() < Edge<3>::se());
	CHECK(Edge<3>::bs() < Edge<3>::nw());
	CHECK(Edge<3>::bs() < Edge<3>::null());

	CHECK_FALSE(Edge<3>::tn() < Edge<3>::bs());
	CHECK_FALSE(Edge<3>::tn() < Edge<3>::tn());
	CHECK(Edge<3>::tn() < Edge<3>::bn());
	CHECK(Edge<3>::tn() < Edge<3>::ts());
	CHECK(Edge<3>::tn() < Edge<3>::bw());
	CHECK(Edge<3>::tn() < Edge<3>::te());
	CHECK(Edge<3>::tn() < Edge<3>::be());
	CHECK(Edge<3>::tn() < Edge<3>::tw());
	CHECK(Edge<3>::tn() < Edge<3>::sw());
	CHECK(Edge<3>::tn() < Edge<3>::ne());
	CHECK(Edge<3>::tn() < Edge<3>::se());
	CHECK(Edge<3>::tn() < Edge<3>::nw());
	CHECK(Edge<3>::tn() < Edge<3>::null());

	CHECK_FALSE(Edge<3>::bn() < Edge<3>::bs());
	CHECK_FALSE(Edge<3>::bn() < Edge<3>::tn());
	CHECK_FALSE(Edge<3>::bn() < Edge<3>::bn());
	CHECK(Edge<3>::bn() < Edge<3>::ts());
	CHECK(Edge<3>::bn() < Edge<3>::bw());
	CHECK(Edge<3>::bn() < Edge<3>::te());
	CHECK(Edge<3>::bn() < Edge<3>::be());
	CHECK(Edge<3>::bn() < Edge<3>::tw());
	CHECK(Edge<3>::bn() < Edge<3>::sw());
	CHECK(Edge<3>::bn() < Edge<3>::ne());
	CHECK(Edge<3>::bn() < Edge<3>::se());
	CHECK(Edge<3>::bn() < Edge<3>::nw());
	CHECK(Edge<3>::bn() < Edge<3>::null());

	CHECK_FALSE(Edge<3>::ts() < Edge<3>::bs());
	CHECK_FALSE(Edge<3>::ts() < Edge<3>::tn());
	CHECK_FALSE(Edge<3>::ts() < Edge<3>::bn());
	CHECK_FALSE(Edge<3>::ts() < Edge<3>::ts());
	CHECK(Edge<3>::ts() < Edge<3>::bw());
	CHECK(Edge<3>::ts() < Edge<3>::te());
	CHECK(Edge<3>::ts() < Edge<3>::be());
	CHECK(Edge<3>::ts() < Edge<3>::tw());
	CHECK(Edge<3>::ts() < Edge<3>::sw());
	CHECK(Edge<3>::ts() < Edge<3>::ne());
	CHECK(Edge<3>::ts() < Edge<3>::se());
	CHECK(Edge<3>::ts() < Edge<3>::nw());
	CHECK(Edge<3>::ts() < Edge<3>::null());

	CHECK_FALSE(Edge<3>::bw() < Edge<3>::bs());
	CHECK_FALSE(Edge<3>::bw() < Edge<3>::tn());
	CHECK_FALSE(Edge<3>::bw() < Edge<3>::bn());
	CHECK_FALSE(Edge<3>::bw() < Edge<3>::ts());
	CHECK_FALSE(Edge<3>::bw() < Edge<3>::bw());
	CHECK(Edge<3>::bw() < Edge<3>::te());
	CHECK(Edge<3>::bw() < Edge<3>::be());
	CHECK(Edge<3>::bw() < Edge<3>::tw());
	CHECK(Edge<3>::bw() < Edge<3>::sw());
	CHECK(Edge<3>::bw() < Edge<3>::ne());
	CHECK(Edge<3>::bw() < Edge<3>::se());
	CHECK(Edge<3>::bw() < Edge<3>::nw());
	CHECK(Edge<3>::bw() < Edge<3>::null());

	CHECK_FALSE(Edge<3>::te() < Edge<3>::bs());
	CHECK_FALSE(Edge<3>::te() < Edge<3>::tn());
	CHECK_FALSE(Edge<3>::te() < Edge<3>::bn());
	CHECK_FALSE(Edge<3>::te() < Edge<3>::ts());
	CHECK_FALSE(Edge<3>::te() < Edge<3>::bw());
	CHECK_FALSE(Edge<3>::te() < Edge<3>::te());
	CHECK(Edge<3>::te() < Edge<3>::be());
	CHECK(Edge<3>::te() < Edge<3>::tw());
	CHECK(Edge<3>::te() < Edge<3>::sw());
	CHECK(Edge<3>::te() < Edge<3>::ne());
	CHECK(Edge<3>::te() < Edge<3>::se());
	CHECK(Edge<3>::te() < Edge<3>::nw());
	CHECK(Edge<3>::te() < Edge<3>::null());

	CHECK_FALSE(Edge<3>::be() < Edge<3>::bs());
	CHECK_FALSE(Edge<3>::be() < Edge<3>::tn());
	CHECK_FALSE(Edge<3>::be() < Edge<3>::bn());
	CHECK_FALSE(Edge<3>::be() < Edge<3>::ts());
	CHECK_FALSE(Edge<3>::be() < Edge<3>::bw());
	CHECK_FALSE(Edge<3>::be() < Edge<3>::te());
	CHECK_FALSE(Edge<3>::be() < Edge<3>::be());
	CHECK(Edge<3>::be() < Edge<3>::tw());
	CHECK(Edge<3>::be() < Edge<3>::sw());
	CHECK(Edge<3>::be() < Edge<3>::ne());
	CHECK(Edge<3>::be() < Edge<3>::se());
	CHECK(Edge<3>::be() < Edge<3>::nw());
	CHECK(Edge<3>::be() < Edge<3>::null());

	CHECK_FALSE(Edge<3>::tw() < Edge<3>::bs());
	CHECK_FALSE(Edge<3>::tw() < Edge<3>::tn());
	CHECK_FALSE(Edge<3>::tw() < Edge<3>::bn());
	CHECK_FALSE(Edge<3>::tw() < Edge<3>::ts());
	CHECK_FALSE(Edge<3>::tw() < Edge<3>::bw());
	CHECK_FALSE(Edge<3>::tw() < Edge<3>::te());
	CHECK_FALSE(Edge<3>::tw() < Edge<3>::be());
	CHECK_FALSE(Edge<3>::tw() < Edge<3>::tw());
	CHECK(Edge<3>::tw() < Edge<3>::sw());
	CHECK(Edge<3>::tw() < Edge<3>::ne());
	CHECK(Edge<3>::tw() < Edge<3>::se());
	CHECK(Edge<3>::tw() < Edge<3>::nw());
	CHECK(Edge<3>::tw() < Edge<3>::null());

	CHECK_FALSE(Edge<3>::sw() < Edge<3>::bs());
	CHECK_FALSE(Edge<3>::sw() < Edge<3>::tn());
	CHECK_FALSE(Edge<3>::sw() < Edge<3>::bn());
	CHECK_FALSE(Edge<3>::sw() < Edge<3>::ts());
	CHECK_FALSE(Edge<3>::sw() < Edge<3>::bw());
	CHECK_FALSE(Edge<3>::sw() < Edge<3>::te());
	CHECK_FALSE(Edge<3>::sw() < Edge<3>::be());
	CHECK_FALSE(Edge<3>::sw() < Edge<3>::tw());
	CHECK_FALSE(Edge<3>::sw() < Edge<3>::sw());
	CHECK(Edge<3>::sw() < Edge<3>::ne());
	CHECK(Edge<3>::sw() < Edge<3>::se());
	CHECK(Edge<3>::sw() < Edge<3>::nw());
	CHECK(Edge<3>::sw() < Edge<3>::null());

	CHECK_FALSE(Edge<3>::ne() < Edge<3>::bs());
	CHECK_FALSE(Edge<3>::ne() < Edge<3>::tn());
	CHECK_FALSE(Edge<3>::ne() < Edge<3>::bn());
	CHECK_FALSE(Edge<3>::ne() < Edge<3>::ts());
	CHECK_FALSE(Edge<3>::ne() < Edge<3>::bw());
	CHECK_FALSE(Edge<3>::ne() < Edge<3>::te());
	CHECK_FALSE(Edge<3>::ne() < Edge<3>::be());
	CHECK_FALSE(Edge<3>::ne() < Edge<3>::tw());
	CHECK_FALSE(Edge<3>::ne() < Edge<3>::sw());
	CHECK_FALSE(Edge<3>::ne() < Edge<3>::ne());
	CHECK(Edge<3>::ne() < Edge<3>::se());
	CHECK(Edge<3>::ne() < Edge<3>::nw());
	CHECK(Edge<3>::ne() < Edge<3>::null());

	CHECK_FALSE(Edge<3>::se() < Edge<3>::bs());
	CHECK_FALSE(Edge<3>::se() < Edge<3>::tn());
	CHECK_FALSE(Edge<3>::se() < Edge<3>::bn());
	CHECK_FALSE(Edge<3>::se() < Edge<3>::ts());
	CHECK_FALSE(Edge<3>::se() < Edge<3>::bw());
	CHECK_FALSE(Edge<3>::se() < Edge<3>::te());
	CHECK_FALSE(Edge<3>::se() < Edge<3>::be());
	CHECK_FALSE(Edge<3>::se() < Edge<3>::tw());
	CHECK_FALSE(Edge<3>::se() < Edge<3>::sw());
	CHECK_FALSE(Edge<3>::se() < Edge<3>::ne());
	CHECK_FALSE(Edge<3>::se() < Edge<3>::se());
	CHECK(Edge<3>::se() < Edge<3>::nw());
	CHECK(Edge<3>::se() < Edge<3>::null());

	CHECK_FALSE(Edge<3>::nw() < Edge<3>::bs());
	CHECK_FALSE(Edge<3>::nw() < Edge<3>::tn());
	CHECK_FALSE(Edge<3>::nw() < Edge<3>::bn());
	CHECK_FALSE(Edge<3>::nw() < Edge<3>::ts());
	CHECK_FALSE(Edge<3>::nw() < Edge<3>::bw());
	CHECK_FALSE(Edge<3>::nw() < Edge<3>::te());
	CHECK_FALSE(Edge<3>::nw() < Edge<3>::be());
	CHECK_FALSE(Edge<3>::nw() < Edge<3>::tw());
	CHECK_FALSE(Edge<3>::nw() < Edge<3>::sw());
	CHECK_FALSE(Edge<3>::nw() < Edge<3>::ne());
	CHECK_FALSE(Edge<3>::nw() < Edge<3>::se());
	CHECK_FALSE(Edge<3>::nw() < Edge<3>::nw());
	CHECK(Edge<3>::nw() < Edge<3>::null());

	CHECK_FALSE(Edge<3>::null() < Edge<3>::bs());
	CHECK_FALSE(Edge<3>::null() < Edge<3>::tn());
	CHECK_FALSE(Edge<3>::null() < Edge<3>::bn());
	CHECK_FALSE(Edge<3>::null() < Edge<3>::ts());
	CHECK_FALSE(Edge<3>::null() < Edge<3>::bw());
	CHECK_FALSE(Edge<3>::null() < Edge<3>::te());
	CHECK_FALSE(Edge<3>::null() < Edge<3>::be());
	CHECK_FALSE(Edge<3>::null() < Edge<3>::tw());
	CHECK_FALSE(Edge<3>::null() < Edge<3>::sw());
	CHECK_FALSE(Edge<3>::null() < Edge<3>::ne());
	CHECK_FALSE(Edge<3>::null() < Edge<3>::se());
	CHECK_FALSE(Edge<3>::null() < Edge<3>::nw());
	CHECK_FALSE(Edge<3>::null() < Edge<3>::null());
}
TEST_CASE("Test ostream for Edge<1>", "[Edge]")
{
	stringstream ss;
	ss << Edge<1>::null();
	CHECK(ss.str() == "Edge<1>::null()");
	ss.str("");
	ss << Edge<1>(13);
	CHECK(ss.str() == "Edge<1> invalid value: 13");
}
TEST_CASE("Test ostream for Edge<2>", "[Edge]")
{
	stringstream ss;
	ss << Edge<2>::null();
	CHECK(ss.str() == "Edge<2>::null()");
	ss.str("");
	ss << Edge<2>(13);
	CHECK(ss.str() == "Edge<2> invalid value: 13");
}
TEST_CASE("Test ostream for Edge<3>", "[Edge]")
{
	stringstream ss;
	ss << Edge<3>::bs();
	CHECK(ss.str() == "Edge<3>::bs()");
	ss.str("");
	ss << Edge<3>::tn();
	CHECK(ss.str() == "Edge<3>::tn()");
	ss.str("");
	ss << Edge<3>::bn();
	CHECK(ss.str() == "Edge<3>::bn()");
	ss.str("");
	ss << Edge<3>::ts();
	CHECK(ss.str() == "Edge<3>::ts()");
	ss.str("");
	ss << Edge<3>::bw();
	CHECK(ss.str() == "Edge<3>::bw()");
	ss.str("");
	ss << Edge<3>::te();
	CHECK(ss.str() == "Edge<3>::te()");
	ss.str("");
	ss << Edge<3>::be();
	CHECK(ss.str() == "Edge<3>::be()");
	ss.str("");
	ss << Edge<3>::tw();
	CHECK(ss.str() == "Edge<3>::tw()");
	ss.str("");
	ss << Edge<3>::sw();
	CHECK(ss.str() == "Edge<3>::sw()");
	ss.str("");
	ss << Edge<3>::ne();
	CHECK(ss.str() == "Edge<3>::ne()");
	ss.str("");
	ss << Edge<3>::se();
	CHECK(ss.str() == "Edge<3>::se()");
	ss.str("");
	ss << Edge<3>::nw();
	CHECK(ss.str() == "Edge<3>::nw()");
	ss.str("");
	ss << Edge<3>::null();
	CHECK(ss.str() == "Edge<3>::null()");
	ss.str("");
	ss << Edge<3>(64);
	CHECK(ss.str() == "Edge<3> invalid value: 64");
}
TEST_CASE("Test iterator for Edge<1>", "[Edge]")
{
	auto iter = Edge<1>::getValues().begin();
	CHECK(*iter == Edge<1>::null());
	CHECK(iter == Edge<1>::getValues().end());
}
TEST_CASE("Test iterator for Edge<2>", "[Edge]")
{
	auto iter = Edge<2>::getValues().begin();
	CHECK(*iter == Edge<2>::null());
	CHECK(iter == Edge<2>::getValues().end());
}
TEST_CASE("Test iterator for Edge<3>", "[Edge]")
{
	auto iter = Edge<3>::getValues().begin();
	CHECK(iter == Edge<3>::getValues().begin());
	CHECK(iter != Edge<3>::getValues().end());
	CHECK(*iter == Edge<3>::bs());
	++iter;
	CHECK(iter->getIndex() == 1);
	CHECK(*iter == Edge<3>::tn());
	++iter;
	CHECK(*iter == Edge<3>::bn());
	++iter;
	CHECK(*iter == Edge<3>::ts());
	++iter;
	CHECK(*iter == Edge<3>::bw());
	++iter;
	CHECK(*iter == Edge<3>::te());
	++iter;
	CHECK(*iter == Edge<3>::be());
	++iter;
	CHECK(*iter == Edge<3>::tw());
	++iter;
	CHECK(*iter == Edge<3>::sw());
	++iter;
	CHECK(*iter == Edge<3>::ne());
	++iter;
	CHECK(*iter == Edge<3>::se());
	++iter;
	CHECK(*iter == Edge<3>::nw());
	++iter;
	CHECK(*iter == Edge<3>::null());
	CHECK(iter == Edge<3>::getValues().end());
}
TEST_CASE("Test from_json for Edge<1>", "[Edge]")
{
	nlohmann::json j;
	j["null"] = nullptr;
	CHECK(j["null"].get<Edge<1>>() == Edge<1>::null());
}
TEST_CASE("Test from_json for Edge<2>", "[Edge]")
{
	nlohmann::json j;
	j["null"] = nullptr;
	CHECK(j["null"].get<Edge<2>>() == Edge<2>::null());
}
TEST_CASE("Test from_json for Edge<3>", "[Edge]")
{
	nlohmann::json j;
	j["null"] = nullptr;
	j["bs"]   = "BS";
	j["tn"]   = "TN";
	j["bn"]   = "BN";
	j["ts"]   = "TS";
	j["bw"]   = "BW";
	j["te"]   = "TE";
	j["be"]   = "BE";
	j["tw"]   = "TW";
	j["sw"]   = "SW";
	j["ne"]   = "NE";
	j["se"]   = "SE";
	j["nw"]   = "NW";
	CHECK(j["null"].get<Edge<3>>() == Edge<3>::null());
	CHECK(j["bs"].get<Edge<3>>() == Edge<3>::bs());
	CHECK(j["tn"].get<Edge<3>>() == Edge<3>::tn());
	CHECK(j["bn"].get<Edge<3>>() == Edge<3>::bn());
	CHECK(j["ts"].get<Edge<3>>() == Edge<3>::ts());
	CHECK(j["bw"].get<Edge<3>>() == Edge<3>::bw());
	CHECK(j["te"].get<Edge<3>>() == Edge<3>::te());
	CHECK(j["be"].get<Edge<3>>() == Edge<3>::be());
	CHECK(j["tw"].get<Edge<3>>() == Edge<3>::tw());
	CHECK(j["sw"].get<Edge<3>>() == Edge<3>::sw());
	CHECK(j["ne"].get<Edge<3>>() == Edge<3>::ne());
	CHECK(j["nw"].get<Edge<3>>() == Edge<3>::nw());
}
TEST_CASE("Test to_json for Edge<1>", "[Edge]")
{
	nlohmann::json j;
	j["null"] = Edge<1>::null();
	CHECK(j["null"] == nullptr);
}
TEST_CASE("Test to_json for Edge<2>", "[Edge]")
{
	nlohmann::json j;
	j["null"] = Edge<2>::null();
	CHECK(j["null"] == nullptr);
}
TEST_CASE("Test to_json for Edge<3>", "[Edge]")
{
	nlohmann::json j;
	j["null"] = Edge<3>::null();
	j["bs"]   = Edge<3>::bs();
	j["tn"]   = Edge<3>::tn();
	j["bn"]   = Edge<3>::bn();
	j["ts"]   = Edge<3>::ts();
	j["bw"]   = Edge<3>::bw();
	j["te"]   = Edge<3>::te();
	j["be"]   = Edge<3>::be();
	j["tw"]   = Edge<3>::tw();
	j["sw"]   = Edge<3>::sw();
	j["ne"]   = Edge<3>::ne();
	j["se"]   = Edge<3>::se();
	j["nw"]   = Edge<3>::nw();
	CHECK(j["null"] == nullptr);
	CHECK(j["bs"] == "BS");
	CHECK(j["tn"] == "TN");
	CHECK(j["bn"] == "BN");
	CHECK(j["ts"] == "TS");
	CHECK(j["bw"] == "BW");
	CHECK(j["te"] == "TE");
	CHECK(j["be"] == "BE");
	CHECK(j["tw"] == "TW");
	CHECK(j["sw"] == "SW");
	CHECK(j["ne"] == "NE");
	CHECK(j["se"] == "SE");
	CHECK(j["nw"] == "NW");
}