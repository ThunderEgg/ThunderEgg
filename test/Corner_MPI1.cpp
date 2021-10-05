#include <ThunderEgg/Face.h>

#include <sstream>

#include <catch2/catch_test_macros.hpp>

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::tpl;

TEST_CASE("Corner<2> num_corners", "[Corner][Face]")
{
	CHECK(Corner<2>::number_of == 4);
}
TEST_CASE("Corner<3> num_corners", "[Corner][Face]")
{
	CHECK(Corner<3>::number_of == 8);
}

TEST_CASE("Corner<2> dimensionality", "[Corner][Face]")
{
	CHECK(Corner<2>::dimensionality == 0);
}
TEST_CASE("Corner<3> dimensionality", "[Corner][Face]")
{
	CHECK(Corner<3>::dimensionality == 0);
}

TEST_CASE("Corner<2> unsigned char constructor works", "[Corner][Face]")
{
	Corner<2> o(13);
	CHECK(o.getIndex() == 13);
}
TEST_CASE("Corner<3> unsigned char constructor works", "[Corner][Face]")
{
	Corner<3> o(13);
	CHECK(o.getIndex() == 13);
}
TEST_CASE("Corner<2> Default constructor works", "[Corner][Face]")
{
	Corner<2> o;
	CHECK(o == Corner<2>::null());
}
TEST_CASE("Corner<3> Default constructor works", "[Corner][Face]")
{
	Corner<3> o;
	CHECK(o == Corner<3>::null());
}
TEST_CASE("Corner<2> named constructors give expected index values", "[Corner][Face]")
{
	CHECK(Corner<2>::sw().getIndex() == 0);
	CHECK(Corner<2>::se().getIndex() == 1);
	CHECK(Corner<2>::nw().getIndex() == 2);
	CHECK(Corner<2>::ne().getIndex() == 3);
	CHECK(Corner<2>::null().getIndex() == 4);
}
TEST_CASE("Corner<3> named constructors give expected index values", "[Corner][Face]")
{
	CHECK(Corner<3>::bsw().getIndex() == 0);
	CHECK(Corner<3>::bse().getIndex() == 1);
	CHECK(Corner<3>::bnw().getIndex() == 2);
	CHECK(Corner<3>::bne().getIndex() == 3);
	CHECK(Corner<3>::tsw().getIndex() == 4);
	CHECK(Corner<3>::tse().getIndex() == 5);
	CHECK(Corner<3>::tnw().getIndex() == 6);
	CHECK(Corner<3>::tne().getIndex() == 7);
	CHECK(Corner<3>::null().getIndex() == 8);
}
TEST_CASE("Corner<2> opposite", "[Corner][Face]")
{
	CHECK(Corner<2>::sw().opposite() == Corner<2>::ne());
	CHECK(Corner<2>::se().opposite() == Corner<2>::nw());
	CHECK(Corner<2>::nw().opposite() == Corner<2>::se());
	CHECK(Corner<2>::ne().opposite() == Corner<2>::sw());
}
TEST_CASE("Corner<3> opposite", "[Corner][Face]")
{
	CHECK(Corner<3>::bsw().opposite() == Corner<3>::tne());
	CHECK(Corner<3>::bse().opposite() == Corner<3>::tnw());
	CHECK(Corner<3>::bnw().opposite() == Corner<3>::tse());
	CHECK(Corner<3>::bne().opposite() == Corner<3>::tsw());
	CHECK(Corner<3>::tsw().opposite() == Corner<3>::bne());
	CHECK(Corner<3>::tse().opposite() == Corner<3>::bnw());
	CHECK(Corner<3>::tnw().opposite() == Corner<3>::bse());
	CHECK(Corner<3>::tne().opposite() == Corner<3>::bsw());
}
TEST_CASE("Corner<2> getSides is as expected", "[Corner][Face]")
{
	{
		auto array = Corner<2>::sw().getSides();
		CHECK(array[0] == Side<2>::west());
		CHECK(array[1] == Side<2>::south());
	}
	{
		auto array = Corner<2>::se().getSides();
		CHECK(array[0] == Side<2>::east());
		CHECK(array[1] == Side<2>::south());
	}
	{
		auto array = Corner<2>::nw().getSides();
		CHECK(array[0] == Side<2>::west());
		CHECK(array[1] == Side<2>::north());
	}
	{
		auto array = Corner<2>::ne().getSides();
		CHECK(array[0] == Side<2>::east());
		CHECK(array[1] == Side<2>::north());
	}
}
TEST_CASE("Corner<3> getSides is as expected", "[Corner][Face]")
{
	{
		auto array = Corner<3>::bsw().getSides();
		CHECK(array[0] == Side<3>::west());
		CHECK(array[1] == Side<3>::south());
		CHECK(array[2] == Side<3>::bottom());
	}
	{
		auto array = Corner<3>::bse().getSides();
		CHECK(array[0] == Side<3>::east());
		CHECK(array[1] == Side<3>::south());
		CHECK(array[2] == Side<3>::bottom());
	}
	{
		auto array = Corner<3>::bnw().getSides();
		CHECK(array[0] == Side<3>::west());
		CHECK(array[1] == Side<3>::north());
		CHECK(array[2] == Side<3>::bottom());
	}
	{
		auto array = Corner<3>::bne().getSides();
		CHECK(array[0] == Side<3>::east());
		CHECK(array[1] == Side<3>::north());
		CHECK(array[2] == Side<3>::bottom());
	}
	{
		auto array = Corner<3>::tsw().getSides();
		CHECK(array[0] == Side<3>::west());
		CHECK(array[1] == Side<3>::south());
		CHECK(array[2] == Side<3>::top());
	}
	{
		auto array = Corner<3>::tse().getSides();
		CHECK(array[0] == Side<3>::east());
		CHECK(array[1] == Side<3>::south());
		CHECK(array[2] == Side<3>::top());
	}
	{
		auto array = Corner<3>::tnw().getSides();
		CHECK(array[0] == Side<3>::west());
		CHECK(array[1] == Side<3>::north());
		CHECK(array[2] == Side<3>::top());
	}
	{
		auto array = Corner<3>::tne().getSides();
		CHECK(array[0] == Side<3>::east());
		CHECK(array[1] == Side<3>::north());
		CHECK(array[2] == Side<3>::top());
	}
}
TEST_CASE("Corner<2> ==", "[Corner][Face]")
{
	CHECK(Corner<2>::sw() == Corner<2>::sw());
	CHECK_FALSE(Corner<2>::sw() == Corner<2>::se());
	CHECK_FALSE(Corner<2>::sw() == Corner<2>::nw());
	CHECK_FALSE(Corner<2>::sw() == Corner<2>::ne());
	CHECK_FALSE(Corner<2>::sw() == Corner<2>::null());

	CHECK_FALSE(Corner<2>::se() == Corner<2>::sw());
	CHECK(Corner<2>::se() == Corner<2>::se());
	CHECK_FALSE(Corner<2>::se() == Corner<2>::nw());
	CHECK_FALSE(Corner<2>::se() == Corner<2>::ne());
	CHECK_FALSE(Corner<2>::se() == Corner<2>::null());

	CHECK_FALSE(Corner<2>::nw() == Corner<2>::sw());
	CHECK_FALSE(Corner<2>::nw() == Corner<2>::se());
	CHECK(Corner<2>::nw() == Corner<2>::nw());
	CHECK_FALSE(Corner<2>::nw() == Corner<2>::ne());
	CHECK_FALSE(Corner<2>::nw() == Corner<2>::null());

	CHECK_FALSE(Corner<2>::ne() == Corner<2>::sw());
	CHECK_FALSE(Corner<2>::ne() == Corner<2>::se());
	CHECK_FALSE(Corner<2>::ne() == Corner<2>::nw());
	CHECK(Corner<2>::ne() == Corner<2>::ne());
	CHECK_FALSE(Corner<2>::ne() == Corner<2>::null());

	CHECK_FALSE(Corner<2>::null() == Corner<2>::sw());
	CHECK_FALSE(Corner<2>::null() == Corner<2>::se());
	CHECK_FALSE(Corner<2>::null() == Corner<2>::nw());
	CHECK_FALSE(Corner<2>::null() == Corner<2>::ne());
	CHECK(Corner<2>::null() == Corner<2>::null());
}
TEST_CASE("Corner<3> ==", "[Corner][Face]")
{
	CHECK(Corner<3>::bsw() == Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::bsw() == Corner<3>::bse());
	CHECK_FALSE(Corner<3>::bsw() == Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::bsw() == Corner<3>::bne());
	CHECK_FALSE(Corner<3>::bsw() == Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::bsw() == Corner<3>::tse());
	CHECK_FALSE(Corner<3>::bsw() == Corner<3>::tnw());
	CHECK_FALSE(Corner<3>::bsw() == Corner<3>::tne());
	CHECK_FALSE(Corner<3>::bsw() == Corner<3>::null());

	CHECK_FALSE(Corner<3>::bse() == Corner<3>::bsw());
	CHECK(Corner<3>::bse() == Corner<3>::bse());
	CHECK_FALSE(Corner<3>::bse() == Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::bse() == Corner<3>::bne());
	CHECK_FALSE(Corner<3>::bse() == Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::bse() == Corner<3>::tse());
	CHECK_FALSE(Corner<3>::bse() == Corner<3>::tnw());
	CHECK_FALSE(Corner<3>::bse() == Corner<3>::tne());
	CHECK_FALSE(Corner<3>::bse() == Corner<3>::null());

	CHECK_FALSE(Corner<3>::bnw() == Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::bnw() == Corner<3>::bse());
	CHECK(Corner<3>::bnw() == Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::bnw() == Corner<3>::bne());
	CHECK_FALSE(Corner<3>::bnw() == Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::bnw() == Corner<3>::tse());
	CHECK_FALSE(Corner<3>::bnw() == Corner<3>::tnw());
	CHECK_FALSE(Corner<3>::bnw() == Corner<3>::tne());
	CHECK_FALSE(Corner<3>::bnw() == Corner<3>::null());

	CHECK_FALSE(Corner<3>::bne() == Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::bne() == Corner<3>::bse());
	CHECK_FALSE(Corner<3>::bne() == Corner<3>::bnw());
	CHECK(Corner<3>::bne() == Corner<3>::bne());
	CHECK_FALSE(Corner<3>::bne() == Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::bne() == Corner<3>::tse());
	CHECK_FALSE(Corner<3>::bne() == Corner<3>::tnw());
	CHECK_FALSE(Corner<3>::bne() == Corner<3>::tne());
	CHECK_FALSE(Corner<3>::bne() == Corner<3>::null());

	CHECK_FALSE(Corner<3>::tsw() == Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::tsw() == Corner<3>::bse());
	CHECK_FALSE(Corner<3>::tsw() == Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::tsw() == Corner<3>::bne());
	CHECK(Corner<3>::tsw() == Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::tsw() == Corner<3>::tse());
	CHECK_FALSE(Corner<3>::tsw() == Corner<3>::tnw());
	CHECK_FALSE(Corner<3>::tsw() == Corner<3>::tne());
	CHECK_FALSE(Corner<3>::tsw() == Corner<3>::null());

	CHECK_FALSE(Corner<3>::tse() == Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::tse() == Corner<3>::bse());
	CHECK_FALSE(Corner<3>::tse() == Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::tse() == Corner<3>::bne());
	CHECK_FALSE(Corner<3>::tse() == Corner<3>::tsw());
	CHECK(Corner<3>::tse() == Corner<3>::tse());
	CHECK_FALSE(Corner<3>::tse() == Corner<3>::tnw());
	CHECK_FALSE(Corner<3>::tse() == Corner<3>::tne());
	CHECK_FALSE(Corner<3>::tse() == Corner<3>::null());

	CHECK_FALSE(Corner<3>::tnw() == Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::tnw() == Corner<3>::bse());
	CHECK_FALSE(Corner<3>::tnw() == Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::tnw() == Corner<3>::bne());
	CHECK_FALSE(Corner<3>::tnw() == Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::tnw() == Corner<3>::tse());
	CHECK(Corner<3>::tnw() == Corner<3>::tnw());
	CHECK_FALSE(Corner<3>::tnw() == Corner<3>::tne());
	CHECK_FALSE(Corner<3>::tnw() == Corner<3>::null());

	CHECK_FALSE(Corner<3>::tne() == Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::tne() == Corner<3>::bse());
	CHECK_FALSE(Corner<3>::tne() == Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::tne() == Corner<3>::bne());
	CHECK_FALSE(Corner<3>::tne() == Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::tne() == Corner<3>::tse());
	CHECK_FALSE(Corner<3>::tne() == Corner<3>::tnw());
	CHECK(Corner<3>::tne() == Corner<3>::tne());
	CHECK_FALSE(Corner<3>::tne() == Corner<3>::null());

	CHECK_FALSE(Corner<3>::null() == Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::null() == Corner<3>::bse());
	CHECK_FALSE(Corner<3>::null() == Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::null() == Corner<3>::bne());
	CHECK_FALSE(Corner<3>::null() == Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::null() == Corner<3>::tse());
	CHECK_FALSE(Corner<3>::null() == Corner<3>::tnw());
	CHECK_FALSE(Corner<3>::null() == Corner<3>::tne());
	CHECK(Corner<3>::null() == Corner<3>::null());
}
TEST_CASE("Corner<2> !=", "[Corner][Face]")
{
	CHECK_FALSE(Corner<2>::sw() != Corner<2>::sw());
	CHECK(Corner<2>::sw() != Corner<2>::se());
	CHECK(Corner<2>::sw() != Corner<2>::nw());
	CHECK(Corner<2>::sw() != Corner<2>::ne());
	CHECK(Corner<2>::sw() != Corner<2>::null());

	CHECK(Corner<2>::se() != Corner<2>::sw());
	CHECK_FALSE(Corner<2>::se() != Corner<2>::se());
	CHECK(Corner<2>::se() != Corner<2>::nw());
	CHECK(Corner<2>::se() != Corner<2>::ne());
	CHECK(Corner<2>::se() != Corner<2>::null());

	CHECK(Corner<2>::nw() != Corner<2>::sw());
	CHECK(Corner<2>::nw() != Corner<2>::se());
	CHECK_FALSE(Corner<2>::nw() != Corner<2>::nw());
	CHECK(Corner<2>::nw() != Corner<2>::ne());
	CHECK(Corner<2>::nw() != Corner<2>::null());

	CHECK(Corner<2>::ne() != Corner<2>::sw());
	CHECK(Corner<2>::ne() != Corner<2>::se());
	CHECK(Corner<2>::ne() != Corner<2>::nw());
	CHECK_FALSE(Corner<2>::ne() != Corner<2>::ne());
	CHECK(Corner<2>::ne() != Corner<2>::null());

	CHECK(Corner<2>::null() != Corner<2>::sw());
	CHECK(Corner<2>::null() != Corner<2>::se());
	CHECK(Corner<2>::null() != Corner<2>::nw());
	CHECK(Corner<2>::null() != Corner<2>::ne());
	CHECK_FALSE(Corner<2>::null() != Corner<2>::null());
}
TEST_CASE("Corner<3> !=", "[Corner][Face]")
{
	CHECK_FALSE(Corner<3>::bsw() != Corner<3>::bsw());
	CHECK(Corner<3>::bsw() != Corner<3>::bse());
	CHECK(Corner<3>::bsw() != Corner<3>::bnw());
	CHECK(Corner<3>::bsw() != Corner<3>::bne());
	CHECK(Corner<3>::bsw() != Corner<3>::tsw());
	CHECK(Corner<3>::bsw() != Corner<3>::tse());
	CHECK(Corner<3>::bsw() != Corner<3>::tnw());
	CHECK(Corner<3>::bsw() != Corner<3>::tne());
	CHECK(Corner<3>::bsw() != Corner<3>::null());

	CHECK(Corner<3>::bse() != Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::bse() != Corner<3>::bse());
	CHECK(Corner<3>::bse() != Corner<3>::bnw());
	CHECK(Corner<3>::bse() != Corner<3>::bne());
	CHECK(Corner<3>::bse() != Corner<3>::tsw());
	CHECK(Corner<3>::bse() != Corner<3>::tse());
	CHECK(Corner<3>::bse() != Corner<3>::tnw());
	CHECK(Corner<3>::bse() != Corner<3>::tne());
	CHECK(Corner<3>::bse() != Corner<3>::null());

	CHECK(Corner<3>::bnw() != Corner<3>::bsw());
	CHECK(Corner<3>::bnw() != Corner<3>::bse());
	CHECK_FALSE(Corner<3>::bnw() != Corner<3>::bnw());
	CHECK(Corner<3>::bnw() != Corner<3>::bne());
	CHECK(Corner<3>::bnw() != Corner<3>::tsw());
	CHECK(Corner<3>::bnw() != Corner<3>::tse());
	CHECK(Corner<3>::bnw() != Corner<3>::tnw());
	CHECK(Corner<3>::bnw() != Corner<3>::tne());
	CHECK(Corner<3>::bnw() != Corner<3>::null());

	CHECK(Corner<3>::bne() != Corner<3>::bsw());
	CHECK(Corner<3>::bne() != Corner<3>::bse());
	CHECK(Corner<3>::bne() != Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::bne() != Corner<3>::bne());
	CHECK(Corner<3>::bne() != Corner<3>::tsw());
	CHECK(Corner<3>::bne() != Corner<3>::tse());
	CHECK(Corner<3>::bne() != Corner<3>::tnw());
	CHECK(Corner<3>::bne() != Corner<3>::tne());
	CHECK(Corner<3>::bne() != Corner<3>::null());

	CHECK(Corner<3>::tsw() != Corner<3>::bsw());
	CHECK(Corner<3>::tsw() != Corner<3>::bse());
	CHECK(Corner<3>::tsw() != Corner<3>::bnw());
	CHECK(Corner<3>::tsw() != Corner<3>::bne());
	CHECK_FALSE(Corner<3>::tsw() != Corner<3>::tsw());
	CHECK(Corner<3>::tsw() != Corner<3>::tse());
	CHECK(Corner<3>::tsw() != Corner<3>::tnw());
	CHECK(Corner<3>::tsw() != Corner<3>::tne());
	CHECK(Corner<3>::tsw() != Corner<3>::null());

	CHECK(Corner<3>::tse() != Corner<3>::bsw());
	CHECK(Corner<3>::tse() != Corner<3>::bse());
	CHECK(Corner<3>::tse() != Corner<3>::bnw());
	CHECK(Corner<3>::tse() != Corner<3>::bne());
	CHECK(Corner<3>::tse() != Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::tse() != Corner<3>::tse());
	CHECK(Corner<3>::tse() != Corner<3>::tnw());
	CHECK(Corner<3>::tse() != Corner<3>::tne());
	CHECK(Corner<3>::tse() != Corner<3>::null());

	CHECK(Corner<3>::tnw() != Corner<3>::bsw());
	CHECK(Corner<3>::tnw() != Corner<3>::bse());
	CHECK(Corner<3>::tnw() != Corner<3>::bnw());
	CHECK(Corner<3>::tnw() != Corner<3>::bne());
	CHECK(Corner<3>::tnw() != Corner<3>::tsw());
	CHECK(Corner<3>::tnw() != Corner<3>::tse());
	CHECK_FALSE(Corner<3>::tnw() != Corner<3>::tnw());
	CHECK(Corner<3>::tnw() != Corner<3>::tne());
	CHECK(Corner<3>::tnw() != Corner<3>::null());

	CHECK(Corner<3>::tne() != Corner<3>::bsw());
	CHECK(Corner<3>::tne() != Corner<3>::bse());
	CHECK(Corner<3>::tne() != Corner<3>::bnw());
	CHECK(Corner<3>::tne() != Corner<3>::bne());
	CHECK(Corner<3>::tne() != Corner<3>::tsw());
	CHECK(Corner<3>::tne() != Corner<3>::tse());
	CHECK(Corner<3>::tne() != Corner<3>::tnw());
	CHECK_FALSE(Corner<3>::tne() != Corner<3>::tne());
	CHECK(Corner<3>::tne() != Corner<3>::null());

	CHECK(Corner<3>::null() != Corner<3>::bsw());
	CHECK(Corner<3>::null() != Corner<3>::bse());
	CHECK(Corner<3>::null() != Corner<3>::bnw());
	CHECK(Corner<3>::null() != Corner<3>::bne());
	CHECK(Corner<3>::null() != Corner<3>::tsw());
	CHECK(Corner<3>::null() != Corner<3>::tse());
	CHECK(Corner<3>::null() != Corner<3>::tnw());
	CHECK(Corner<3>::null() != Corner<3>::tne());
	CHECK_FALSE(Corner<3>::null() != Corner<3>::null());
}
TEST_CASE("Corner<2> <", "[Corner][Face]")
{
	CHECK_FALSE(Corner<2>::sw() < Corner<2>::sw());
	CHECK(Corner<2>::sw() < Corner<2>::se());
	CHECK(Corner<2>::sw() < Corner<2>::nw());
	CHECK(Corner<2>::sw() < Corner<2>::ne());
	CHECK(Corner<2>::sw() < Corner<2>::null());

	CHECK_FALSE(Corner<2>::se() < Corner<2>::sw());
	CHECK_FALSE(Corner<2>::se() < Corner<2>::se());
	CHECK(Corner<2>::se() < Corner<2>::nw());
	CHECK(Corner<2>::se() < Corner<2>::ne());
	CHECK(Corner<2>::se() < Corner<2>::null());

	CHECK_FALSE(Corner<2>::nw() < Corner<2>::sw());
	CHECK_FALSE(Corner<2>::nw() < Corner<2>::se());
	CHECK_FALSE(Corner<2>::nw() < Corner<2>::nw());
	CHECK(Corner<2>::nw() < Corner<2>::ne());
	CHECK(Corner<2>::nw() < Corner<2>::null());

	CHECK_FALSE(Corner<2>::ne() < Corner<2>::sw());
	CHECK_FALSE(Corner<2>::ne() < Corner<2>::se());
	CHECK_FALSE(Corner<2>::ne() < Corner<2>::nw());
	CHECK_FALSE(Corner<2>::ne() < Corner<2>::ne());
	CHECK(Corner<2>::ne() < Corner<2>::null());

	CHECK_FALSE(Corner<2>::null() < Corner<2>::sw());
	CHECK_FALSE(Corner<2>::null() < Corner<2>::se());
	CHECK_FALSE(Corner<2>::null() < Corner<2>::nw());
	CHECK_FALSE(Corner<2>::null() < Corner<2>::ne());
	CHECK_FALSE(Corner<2>::null() < Corner<2>::null());
}
TEST_CASE("Corner<3> <", "[Corner][Face]")
{
	CHECK_FALSE(Corner<3>::bsw() < Corner<3>::bsw());
	CHECK(Corner<3>::bsw() < Corner<3>::bse());
	CHECK(Corner<3>::bsw() < Corner<3>::bnw());
	CHECK(Corner<3>::bsw() < Corner<3>::bne());
	CHECK(Corner<3>::bsw() < Corner<3>::tsw());
	CHECK(Corner<3>::bsw() < Corner<3>::tse());
	CHECK(Corner<3>::bsw() < Corner<3>::tnw());
	CHECK(Corner<3>::bsw() < Corner<3>::tne());
	CHECK(Corner<3>::bsw() < Corner<3>::null());

	CHECK_FALSE(Corner<3>::bse() < Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::bse() < Corner<3>::bse());
	CHECK(Corner<3>::bse() < Corner<3>::bnw());
	CHECK(Corner<3>::bse() < Corner<3>::bne());
	CHECK(Corner<3>::bse() < Corner<3>::tsw());
	CHECK(Corner<3>::bse() < Corner<3>::tse());
	CHECK(Corner<3>::bse() < Corner<3>::tnw());
	CHECK(Corner<3>::bse() < Corner<3>::tne());
	CHECK(Corner<3>::bse() < Corner<3>::null());

	CHECK_FALSE(Corner<3>::bnw() < Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::bnw() < Corner<3>::bse());
	CHECK_FALSE(Corner<3>::bnw() < Corner<3>::bnw());
	CHECK(Corner<3>::bnw() < Corner<3>::bne());
	CHECK(Corner<3>::bnw() < Corner<3>::tsw());
	CHECK(Corner<3>::bnw() < Corner<3>::tse());
	CHECK(Corner<3>::bnw() < Corner<3>::tnw());
	CHECK(Corner<3>::bnw() < Corner<3>::tne());
	CHECK(Corner<3>::bnw() < Corner<3>::null());

	CHECK_FALSE(Corner<3>::bne() < Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::bne() < Corner<3>::bse());
	CHECK_FALSE(Corner<3>::bne() < Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::bne() < Corner<3>::bne());
	CHECK(Corner<3>::bne() < Corner<3>::tsw());
	CHECK(Corner<3>::bne() < Corner<3>::tse());
	CHECK(Corner<3>::bne() < Corner<3>::tnw());
	CHECK(Corner<3>::bne() < Corner<3>::tne());
	CHECK(Corner<3>::bne() < Corner<3>::null());

	CHECK_FALSE(Corner<3>::tsw() < Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::tsw() < Corner<3>::bse());
	CHECK_FALSE(Corner<3>::tsw() < Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::tsw() < Corner<3>::bne());
	CHECK_FALSE(Corner<3>::tsw() < Corner<3>::tsw());
	CHECK(Corner<3>::tsw() < Corner<3>::tse());
	CHECK(Corner<3>::tsw() < Corner<3>::tnw());
	CHECK(Corner<3>::tsw() < Corner<3>::tne());
	CHECK(Corner<3>::tsw() < Corner<3>::null());

	CHECK_FALSE(Corner<3>::tse() < Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::tse() < Corner<3>::bse());
	CHECK_FALSE(Corner<3>::tse() < Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::tse() < Corner<3>::bne());
	CHECK_FALSE(Corner<3>::tse() < Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::tse() < Corner<3>::tse());
	CHECK(Corner<3>::tse() < Corner<3>::tnw());
	CHECK(Corner<3>::tse() < Corner<3>::tne());
	CHECK(Corner<3>::tse() < Corner<3>::null());

	CHECK_FALSE(Corner<3>::tnw() < Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::tnw() < Corner<3>::bse());
	CHECK_FALSE(Corner<3>::tnw() < Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::tnw() < Corner<3>::bne());
	CHECK_FALSE(Corner<3>::tnw() < Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::tnw() < Corner<3>::tse());
	CHECK_FALSE(Corner<3>::tnw() < Corner<3>::tnw());
	CHECK(Corner<3>::tnw() < Corner<3>::tne());
	CHECK(Corner<3>::tnw() < Corner<3>::null());

	CHECK_FALSE(Corner<3>::tne() < Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::tne() < Corner<3>::bse());
	CHECK_FALSE(Corner<3>::tne() < Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::tne() < Corner<3>::bne());
	CHECK_FALSE(Corner<3>::tne() < Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::tne() < Corner<3>::tse());
	CHECK_FALSE(Corner<3>::tne() < Corner<3>::tnw());
	CHECK_FALSE(Corner<3>::tne() < Corner<3>::tne());
	CHECK(Corner<3>::tne() < Corner<3>::null());

	CHECK_FALSE(Corner<3>::null() < Corner<3>::bsw());
	CHECK_FALSE(Corner<3>::null() < Corner<3>::bse());
	CHECK_FALSE(Corner<3>::null() < Corner<3>::bnw());
	CHECK_FALSE(Corner<3>::null() < Corner<3>::bne());
	CHECK_FALSE(Corner<3>::null() < Corner<3>::tsw());
	CHECK_FALSE(Corner<3>::null() < Corner<3>::tse());
	CHECK_FALSE(Corner<3>::null() < Corner<3>::tnw());
	CHECK_FALSE(Corner<3>::null() < Corner<3>::tne());
	CHECK_FALSE(Corner<3>::null() < Corner<3>::null());
}
TEST_CASE("Test ostream for Corner<2>", "[Corner][Face]")
{
	stringstream ss;
	ss << Corner<2>::sw();
	CHECK(ss.str() == "Corner<2>::sw()");
	ss.str("");
	ss << Corner<2>::se();
	CHECK(ss.str() == "Corner<2>::se()");
	ss.str("");
	ss << Corner<2>::nw();
	CHECK(ss.str() == "Corner<2>::nw()");
	ss.str("");
	ss << Corner<2>::ne();
	CHECK(ss.str() == "Corner<2>::ne()");
	ss.str("");
	ss << Corner<2>::null();
	CHECK(ss.str() == "Corner<2>::null()");
	ss.str("");
	ss << Corner<2>(13);
	CHECK(ss.str() == "Corner<2> invalid value: 13");
}
TEST_CASE("Test ostream for Corner<3>", "[Corner][Face]")
{
	stringstream ss;
	ss << Corner<3>::bsw();
	CHECK(ss.str() == "Corner<3>::bsw()");
	ss.str("");
	ss << Corner<3>::bse();
	CHECK(ss.str() == "Corner<3>::bse()");
	ss.str("");
	ss << Corner<3>::bnw();
	CHECK(ss.str() == "Corner<3>::bnw()");
	ss.str("");
	ss << Corner<3>::bne();
	CHECK(ss.str() == "Corner<3>::bne()");
	ss.str("");
	ss << Corner<3>::tsw();
	CHECK(ss.str() == "Corner<3>::tsw()");
	ss.str("");
	ss << Corner<3>::tse();
	CHECK(ss.str() == "Corner<3>::tse()");
	ss.str("");
	ss << Corner<3>::tnw();
	CHECK(ss.str() == "Corner<3>::tnw()");
	ss.str("");
	ss << Corner<3>::tne();
	CHECK(ss.str() == "Corner<3>::tne()");
	ss.str("");
	ss << Corner<3>::null();
	CHECK(ss.str() == "Corner<3>::null()");
	ss.str("");
	ss << Corner<3>(13);
	CHECK(ss.str() == "Corner<3> invalid value: 13");
}
TEST_CASE("Test iterator for Corner<2>", "[Corner][Face]")
{
	auto iter = Corner<2>::getValues().begin();
	CHECK(iter == Corner<2>::getValues().begin());
	CHECK(iter != Corner<2>::getValues().end());
	CHECK(*iter == Corner<2>::sw());
	++iter;
	CHECK(iter->getIndex() == 1);
	CHECK(*iter == Corner<2>::se());
	++iter;
	CHECK(*iter == Corner<2>::nw());
	++iter;
	CHECK(*iter == Corner<2>::ne());
	++iter;
	CHECK(*iter == Corner<2>::null());
	CHECK(iter == Corner<2>::getValues().end());
}
TEST_CASE("Test iterator for Corner<3>", "[Corner][Face]")
{
	auto iter = Corner<3>::getValues().begin();
	CHECK(iter == Corner<3>::getValues().begin());
	CHECK(iter != Corner<3>::getValues().end());
	CHECK(*iter == Corner<3>::bsw());
	++iter;
	CHECK(iter->getIndex() == 1);
	CHECK(*iter == Corner<3>::bse());
	++iter;
	CHECK(*iter == Corner<3>::bnw());
	++iter;
	CHECK(*iter == Corner<3>::bne());
	++iter;
	CHECK(*iter == Corner<3>::tsw());
	++iter;
	CHECK(*iter == Corner<3>::tse());
	++iter;
	CHECK(*iter == Corner<3>::tnw());
	++iter;
	CHECK(*iter == Corner<3>::tne());
	++iter;
	CHECK(*iter == Corner<3>::null());
	CHECK(iter == Corner<3>::getValues().end());
}
TEST_CASE("Test from_json for Corner<2>", "[Corner][Face]")
{
	nlohmann::json j;
	j["null"] = nullptr;
	j["sw"]   = "SW";
	j["se"]   = "SE";
	j["nw"]   = "NW";
	j["ne"]   = "NE";
	CHECK(j["null"].get<Corner<2>>() == Corner<2>::null());
	CHECK(j["sw"].get<Corner<2>>() == Corner<2>::sw());
	CHECK(j["se"].get<Corner<2>>() == Corner<2>::se());
	CHECK(j["nw"].get<Corner<2>>() == Corner<2>::nw());
	CHECK(j["ne"].get<Corner<2>>() == Corner<2>::ne());
}
TEST_CASE("Test from_json for Corner<3>", "[Corner][Face]")
{
	nlohmann::json j;
	j["null"] = nullptr;
	j["bsw"]  = "BSW";
	j["bse"]  = "BSE";
	j["bnw"]  = "BNW";
	j["bne"]  = "BNE";
	j["tsw"]  = "TSW";
	j["tse"]  = "TSE";
	j["tnw"]  = "TNW";
	j["tne"]  = "TNE";
	CHECK(j["null"].get<Corner<3>>() == Corner<3>::null());
	CHECK(j["bsw"].get<Corner<3>>() == Corner<3>::bsw());
	CHECK(j["bse"].get<Corner<3>>() == Corner<3>::bse());
	CHECK(j["bnw"].get<Corner<3>>() == Corner<3>::bnw());
	CHECK(j["bne"].get<Corner<3>>() == Corner<3>::bne());
	CHECK(j["tsw"].get<Corner<3>>() == Corner<3>::tsw());
	CHECK(j["tse"].get<Corner<3>>() == Corner<3>::tse());
	CHECK(j["tnw"].get<Corner<3>>() == Corner<3>::tnw());
	CHECK(j["tne"].get<Corner<3>>() == Corner<3>::tne());
}
TEST_CASE("Test to_json for Corner<2>", "[Corner][Face]")
{
	nlohmann::json j;
	j["null"] = Corner<2>::null();
	j["sw"]   = Corner<2>::sw();
	j["se"]   = Corner<2>::se();
	j["nw"]   = Corner<2>::nw();
	j["ne"]   = Corner<2>::ne();
	CHECK(j["null"] == nullptr);
	CHECK(j["sw"] == "SW");
	CHECK(j["se"] == "SE");
	CHECK(j["nw"] == "NW");
	CHECK(j["ne"] == "NE");
}
TEST_CASE("Test to_json for Corner<3>", "[Corner][Face]")
{
	nlohmann::json j;
	j["null"] = Corner<3>::null();
	j["bsw"]  = Corner<3>::bsw();
	j["bse"]  = Corner<3>::bse();
	j["bnw"]  = Corner<3>::bnw();
	j["bne"]  = Corner<3>::bne();
	j["tsw"]  = Corner<3>::tsw();
	j["tse"]  = Corner<3>::tse();
	j["tnw"]  = Corner<3>::tnw();
	j["tne"]  = Corner<3>::tne();
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