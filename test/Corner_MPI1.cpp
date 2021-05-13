#include <ThunderEgg/Corner.h>

#include <sstream>

#include <catch2/catch_test_macros.hpp>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Corner<1> num_corners", "[Corner]")
{
	CHECK(Corner<1>::num_corners == 0);
}
TEST_CASE("Corner<2> num_corners", "[Corner]")
{
	CHECK(Corner<2>::num_corners == 4);
}
TEST_CASE("Corner<3> num_corners", "[Corner]")
{
	CHECK(Corner<3>::num_corners == 8);
}

TEST_CASE("Corner<1> dimensionality", "[Corner]")
{
	CHECK(Corner<1>::dimensionality == 0);
}
TEST_CASE("Corner<2> dimensionality", "[Corner]")
{
	CHECK(Corner<2>::dimensionality == 0);
}
TEST_CASE("Corner<3> dimensionality", "[Corner]")
{
	CHECK(Corner<3>::dimensionality == 0);
}

TEST_CASE("Corner<1> unsigned char constructor works", "[Corner]")
{
	Corner<1> o(13);
	CHECK(o.getIndex() == 13);
}
TEST_CASE("Corner<2> unsigned char constructor works", "[Corner]")
{
	Corner<2> o(13);
	CHECK(o.getIndex() == 13);
}
TEST_CASE("Corner<3> unsigned char constructor works", "[Corner]")
{
	Corner<3> o(13);
	CHECK(o.getIndex() == 13);
}
TEST_CASE("Corner<1> Default constructor works", "[Corner]")
{
	Corner<1> o;
	CHECK(o == Corner<1>::null());
}
TEST_CASE("Corner<2> Default constructor works", "[Corner]")
{
	Corner<2> o;
	CHECK(o == Corner<2>::null());
}
TEST_CASE("Corner<3> Default constructor works", "[Corner]")
{
	Corner<3> o;
	CHECK(o == Corner<3>::null());
}
TEST_CASE("Corner<1> named constructors give expected index values", "[Corner]")
{
	CHECK(Corner<1>::null().getIndex() == 0);
}
TEST_CASE("Corner<2> named constructors give expected index values", "[Corner]")
{
	CHECK(Corner<2>::sw().getIndex() == 0);
	CHECK(Corner<2>::se().getIndex() == 1);
	CHECK(Corner<2>::nw().getIndex() == 2);
	CHECK(Corner<2>::ne().getIndex() == 3);
	CHECK(Corner<2>::null().getIndex() == 4);
}
TEST_CASE("Corner<3> named constructors give expected index values", "[Corner]")
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
TEST_CASE("Corner<2> opposite", "[Corner]")
{
	CHECK(Corner<2>::sw().opposite() == Corner<2>::ne());
	CHECK(Corner<2>::se().opposite() == Corner<2>::nw());
	CHECK(Corner<2>::nw().opposite() == Corner<2>::se());
	CHECK(Corner<2>::ne().opposite() == Corner<2>::sw());
}
TEST_CASE("Corner<3> opposite", "[Corner]")
{
	CHECK(Corner<3>::bsw().opposite() == Corner<2>::tne());
	CHECK(Corner<3>::bse().opposite() == Corner<2>::tnw());
	CHECK(Corner<3>::bnw().opposite() == Corner<2>::tse());
	CHECK(Corner<3>::bne().opposite() == Corner<2>::tsw());
	CHECK(Corner<3>::tsw().opposite() == Corner<2>::bne());
	CHECK(Corner<3>::tse().opposite() == Corner<2>::bnw());
	CHECK(Corner<3>::tnw().opposite() == Corner<2>::bse());
	CHECK(Corner<3>::tne().opposite() == Corner<2>::bsw());
}
TEST_CASE("Corner<2> getNbrOnSide is as expected", "[Corner]")
{
	CHECK(Corner<2>::sw().getNbrOnSide(Side<2>::west()) == Corner<2>::se());
	CHECK(Corner<2>::sw().getNbrOnSide(Side<2>::east()) == Corner<2>::se());
	CHECK(Corner<2>::sw().getNbrOnSide(Side<2>::south()) == Corner<2>::nw());
	CHECK(Corner<2>::sw().getNbrOnSide(Side<2>::north()) == Corner<2>::nw());

	CHECK(Corner<2>::se().getNbrOnSide(Side<2>::west()) == Corner<2>::sw());
	CHECK(Corner<2>::se().getNbrOnSide(Side<2>::east()) == Corner<2>::sw());
	CHECK(Corner<2>::se().getNbrOnSide(Side<2>::south()) == Corner<2>::ne());
	CHECK(Corner<2>::se().getNbrOnSide(Side<2>::north()) == Corner<2>::ne());

	CHECK(Corner<2>::nw().getNbrOnSide(Side<2>::west()) == Corner<2>::ne());
	CHECK(Corner<2>::nw().getNbrOnSide(Side<2>::east()) == Corner<2>::ne());
	CHECK(Corner<2>::nw().getNbrOnSide(Side<2>::south()) == Corner<2>::sw());
	CHECK(Corner<2>::nw().getNbrOnSide(Side<2>::north()) == Corner<2>::sw());

	CHECK(Corner<2>::ne().getNbrOnSide(Side<2>::west()) == Corner<2>::nw());
	CHECK(Corner<2>::ne().getNbrOnSide(Side<2>::east()) == Corner<2>::nw());
	CHECK(Corner<2>::ne().getNbrOnSide(Side<2>::south()) == Corner<2>::se());
	CHECK(Corner<2>::ne().getNbrOnSide(Side<2>::north()) == Corner<2>::se());
}
TEST_CASE("Corner<3> getNbrOnSide is as expected", "[Corner]")
{
	CHECK(Corner<3>::bsw().getNbrOnSide(Side<3>::west()) == Corner<3>::bse());
	CHECK(Corner<3>::bsw().getNbrOnSide(Side<3>::east()) == Corner<3>::bse());
	CHECK(Corner<3>::bsw().getNbrOnSide(Side<3>::south()) == Corner<3>::bnw());
	CHECK(Corner<3>::bsw().getNbrOnSide(Side<3>::north()) == Corner<3>::bnw());
	CHECK(Corner<3>::bsw().getNbrOnSide(Side<3>::bottom()) == Corner<3>::tsw());
	CHECK(Corner<3>::bsw().getNbrOnSide(Side<3>::top()) == Corner<3>::tsw());

	CHECK(Corner<3>::bse().getNbrOnSide(Side<3>::west()) == Corner<3>::bsw());
	CHECK(Corner<3>::bse().getNbrOnSide(Side<3>::east()) == Corner<3>::bsw());
	CHECK(Corner<3>::bse().getNbrOnSide(Side<3>::south()) == Corner<3>::bne());
	CHECK(Corner<3>::bse().getNbrOnSide(Side<3>::north()) == Corner<3>::bne());
	CHECK(Corner<3>::bse().getNbrOnSide(Side<3>::bottom()) == Corner<3>::tse());
	CHECK(Corner<3>::bse().getNbrOnSide(Side<3>::top()) == Corner<3>::tse());

	CHECK(Corner<3>::bnw().getNbrOnSide(Side<3>::west()) == Corner<3>::bne());
	CHECK(Corner<3>::bnw().getNbrOnSide(Side<3>::east()) == Corner<3>::bne());
	CHECK(Corner<3>::bnw().getNbrOnSide(Side<3>::south()) == Corner<3>::bsw());
	CHECK(Corner<3>::bnw().getNbrOnSide(Side<3>::north()) == Corner<3>::bsw());
	CHECK(Corner<3>::bnw().getNbrOnSide(Side<3>::bottom()) == Corner<3>::tnw());
	CHECK(Corner<3>::bnw().getNbrOnSide(Side<3>::top()) == Corner<3>::tnw());

	CHECK(Corner<3>::bne().getNbrOnSide(Side<3>::west()) == Corner<3>::bnw());
	CHECK(Corner<3>::bne().getNbrOnSide(Side<3>::east()) == Corner<3>::bnw());
	CHECK(Corner<3>::bne().getNbrOnSide(Side<3>::south()) == Corner<3>::bse());
	CHECK(Corner<3>::bne().getNbrOnSide(Side<3>::north()) == Corner<3>::bse());
	CHECK(Corner<3>::bne().getNbrOnSide(Side<3>::bottom()) == Corner<3>::tne());
	CHECK(Corner<3>::bne().getNbrOnSide(Side<3>::top()) == Corner<3>::tne());

	CHECK(Corner<3>::tsw().getNbrOnSide(Side<3>::west()) == Corner<3>::tse());
	CHECK(Corner<3>::tsw().getNbrOnSide(Side<3>::east()) == Corner<3>::tse());
	CHECK(Corner<3>::tsw().getNbrOnSide(Side<3>::south()) == Corner<3>::tnw());
	CHECK(Corner<3>::tsw().getNbrOnSide(Side<3>::north()) == Corner<3>::tnw());
	CHECK(Corner<3>::tsw().getNbrOnSide(Side<3>::bottom()) == Corner<3>::bsw());
	CHECK(Corner<3>::tsw().getNbrOnSide(Side<3>::top()) == Corner<3>::bsw());

	CHECK(Corner<3>::tse().getNbrOnSide(Side<3>::east()) == Corner<3>::tsw());
	CHECK(Corner<3>::tse().getNbrOnSide(Side<3>::west()) == Corner<3>::tsw());
	CHECK(Corner<3>::tse().getNbrOnSide(Side<3>::south()) == Corner<3>::tne());
	CHECK(Corner<3>::tse().getNbrOnSide(Side<3>::north()) == Corner<3>::tne());
	CHECK(Corner<3>::tse().getNbrOnSide(Side<3>::bottom()) == Corner<3>::bse());
	CHECK(Corner<3>::tse().getNbrOnSide(Side<3>::top()) == Corner<3>::bse());

	CHECK(Corner<3>::tnw().getNbrOnSide(Side<3>::west()) == Corner<3>::tne());
	CHECK(Corner<3>::tnw().getNbrOnSide(Side<3>::east()) == Corner<3>::tne());
	CHECK(Corner<3>::tnw().getNbrOnSide(Side<3>::south()) == Corner<3>::tsw());
	CHECK(Corner<3>::tnw().getNbrOnSide(Side<3>::north()) == Corner<3>::tsw());
	CHECK(Corner<3>::tnw().getNbrOnSide(Side<3>::bottom()) == Corner<3>::bnw());
	CHECK(Corner<3>::tnw().getNbrOnSide(Side<3>::top()) == Corner<3>::bnw());

	CHECK(Corner<3>::tne().getNbrOnSide(Side<3>::east()) == Corner<3>::tnw());
	CHECK(Corner<3>::tne().getNbrOnSide(Side<3>::west()) == Corner<3>::tnw());
	CHECK(Corner<3>::tne().getNbrOnSide(Side<3>::south()) == Corner<3>::tse());
	CHECK(Corner<3>::tne().getNbrOnSide(Side<3>::north()) == Corner<3>::tse());
	CHECK(Corner<3>::tne().getNbrOnSide(Side<3>::bottom()) == Corner<3>::bne());
	CHECK(Corner<3>::tne().getNbrOnSide(Side<3>::top()) == Corner<3>::bne());
}
TEST_CASE("Corner<2> isLowerOnAxis ", "[Corner]")
{
	CHECK(Corner<2>::sw().isLowerOnAxis(0));
	CHECK(Corner<2>::sw().isLowerOnAxis(1));
	CHECK_FALSE(Corner<2>::se().isLowerOnAxis(0));
	CHECK(Corner<2>::se().isLowerOnAxis(1));
	CHECK(Corner<2>::nw().isLowerOnAxis(0));
	CHECK_FALSE(Corner<2>::nw().isLowerOnAxis(1));
	CHECK_FALSE(Corner<2>::ne().isLowerOnAxis(0));
	CHECK_FALSE(Corner<2>::ne().isLowerOnAxis(1));
}
TEST_CASE("Corner<3> isLowerOnAxis ", "[Corner]")
{
	CHECK(Corner<3>::bsw().isLowerOnAxis(0));
	CHECK(Corner<3>::bsw().isLowerOnAxis(1));
	CHECK(Corner<3>::bsw().isLowerOnAxis(2));

	CHECK_FALSE(Corner<3>::bse().isLowerOnAxis(0));
	CHECK(Corner<3>::bse().isLowerOnAxis(1));
	CHECK(Corner<3>::bse().isLowerOnAxis(2));

	CHECK(Corner<3>::bnw().isLowerOnAxis(0));
	CHECK_FALSE(Corner<3>::bnw().isLowerOnAxis(1));
	CHECK(Corner<3>::bnw().isLowerOnAxis(2));

	CHECK_FALSE(Corner<3>::bne().isLowerOnAxis(0));
	CHECK_FALSE(Corner<3>::bne().isLowerOnAxis(1));
	CHECK(Corner<3>::bne().isLowerOnAxis(2));

	CHECK(Corner<3>::tsw().isLowerOnAxis(0));
	CHECK(Corner<3>::tsw().isLowerOnAxis(1));
	CHECK_FALSE(Corner<3>::tsw().isLowerOnAxis(2));

	CHECK_FALSE(Corner<3>::tse().isLowerOnAxis(0));
	CHECK(Corner<3>::tse().isLowerOnAxis(1));
	CHECK_FALSE(Corner<3>::tse().isLowerOnAxis(2));

	CHECK(Corner<3>::tnw().isLowerOnAxis(0));
	CHECK_FALSE(Corner<3>::tnw().isLowerOnAxis(1));
	CHECK_FALSE(Corner<3>::tnw().isLowerOnAxis(2));

	CHECK_FALSE(Corner<3>::tne().isLowerOnAxis(0));
	CHECK_FALSE(Corner<3>::tne().isLowerOnAxis(1));
	CHECK_FALSE(Corner<3>::tne().isLowerOnAxis(2));
}
TEST_CASE("Corner<2> isHigherOnAxis ", "[Corner]")
{
	CHECK_FALSE(Corner<2>::sw().isHigherOnAxis(0));
	CHECK_FALSE(Corner<2>::sw().isHigherOnAxis(1));
	CHECK(Corner<2>::se().isHigherOnAxis(0));
	CHECK_FALSE(Corner<2>::se().isHigherOnAxis(1));
	CHECK_FALSE(Corner<2>::nw().isHigherOnAxis(0));
	CHECK(Corner<2>::nw().isHigherOnAxis(1));
	CHECK(Corner<2>::ne().isHigherOnAxis(0));
	CHECK(Corner<2>::ne().isHigherOnAxis(1));
}
TEST_CASE("Corner<3> isHigherOnAxis ", "[Corner]")
{
	CHECK_FALSE(Corner<3>::bsw().isHigherOnAxis(0));
	CHECK_FALSE(Corner<3>::bsw().isHigherOnAxis(1));
	CHECK_FALSE(Corner<3>::bsw().isHigherOnAxis(2));

	CHECK(Corner<3>::bse().isHigherOnAxis(0));
	CHECK_FALSE(Corner<3>::bse().isHigherOnAxis(1));
	CHECK_FALSE(Corner<3>::bse().isHigherOnAxis(2));

	CHECK_FALSE(Corner<3>::bnw().isHigherOnAxis(0));
	CHECK(Corner<3>::bnw().isHigherOnAxis(1));
	CHECK_FALSE(Corner<3>::bnw().isHigherOnAxis(2));

	CHECK(Corner<3>::bne().isHigherOnAxis(0));
	CHECK(Corner<3>::bne().isHigherOnAxis(1));
	CHECK_FALSE(Corner<3>::bne().isHigherOnAxis(2));

	CHECK_FALSE(Corner<3>::tsw().isHigherOnAxis(0));
	CHECK_FALSE(Corner<3>::tsw().isHigherOnAxis(1));
	CHECK(Corner<3>::tsw().isHigherOnAxis(2));

	CHECK(Corner<3>::tse().isHigherOnAxis(0));
	CHECK_FALSE(Corner<3>::tse().isHigherOnAxis(1));
	CHECK(Corner<3>::tse().isHigherOnAxis(2));

	CHECK_FALSE(Corner<3>::tnw().isHigherOnAxis(0));
	CHECK(Corner<3>::tnw().isHigherOnAxis(1));
	CHECK(Corner<3>::tnw().isHigherOnAxis(2));

	CHECK(Corner<3>::tne().isHigherOnAxis(0));
	CHECK(Corner<3>::tne().isHigherOnAxis(1));
	CHECK(Corner<3>::tne().isHigherOnAxis(2));
}
TEST_CASE("Corner<2> getInteriorSides is as expected", "[Corner]")
{
	{
		auto array = Corner<2>::sw().getInteriorSides();
		CHECK(array[0] == Side<2>::east());
		CHECK(array[1] == Side<2>::north());
	}
	{
		auto array = Corner<2>::se().getInteriorSides();
		CHECK(array[0] == Side<2>::west());
		CHECK(array[1] == Side<2>::north());
	}
	{
		auto array = Corner<2>::nw().getInteriorSides();
		CHECK(array[0] == Side<2>::east());
		CHECK(array[1] == Side<2>::south());
	}
	{
		auto array = Corner<2>::ne().getInteriorSides();
		CHECK(array[0] == Side<2>::west());
		CHECK(array[1] == Side<2>::south());
	}
}
TEST_CASE("Corner<3> getInteriorSides is as expected", "[Corner]")
{
	{
		auto array = Corner<3>::bsw().getInteriorSides();
		CHECK(array[0] == Side<3>::east());
		CHECK(array[1] == Side<3>::north());
		CHECK(array[2] == Side<3>::top());
	}
	{
		auto array = Corner<3>::bse().getInteriorSides();
		CHECK(array[0] == Side<3>::west());
		CHECK(array[1] == Side<3>::north());
		CHECK(array[2] == Side<3>::top());
	}
	{
		auto array = Corner<3>::bnw().getInteriorSides();
		CHECK(array[0] == Side<3>::east());
		CHECK(array[1] == Side<3>::south());
		CHECK(array[2] == Side<3>::top());
	}
	{
		auto array = Corner<3>::bne().getInteriorSides();
		CHECK(array[0] == Side<3>::west());
		CHECK(array[1] == Side<3>::south());
		CHECK(array[2] == Side<3>::top());
	}
	{
		auto array = Corner<3>::tsw().getInteriorSides();
		CHECK(array[0] == Side<3>::east());
		CHECK(array[1] == Side<3>::north());
		CHECK(array[2] == Side<3>::bottom());
	}
	{
		auto array = Corner<3>::tse().getInteriorSides();
		CHECK(array[0] == Side<3>::west());
		CHECK(array[1] == Side<3>::north());
		CHECK(array[2] == Side<3>::bottom());
	}
	{
		auto array = Corner<3>::tnw().getInteriorSides();
		CHECK(array[0] == Side<3>::east());
		CHECK(array[1] == Side<3>::south());
		CHECK(array[2] == Side<3>::bottom());
	}
	{
		auto array = Corner<3>::tne().getInteriorSides();
		CHECK(array[0] == Side<3>::west());
		CHECK(array[1] == Side<3>::south());
		CHECK(array[2] == Side<3>::bottom());
	}
}
TEST_CASE("Corner<2> getExteriorSides is as expected", "[Corner]")
{
	{
		auto array = Corner<2>::sw().getExteriorSides();
		CHECK(array[0] == Side<2>::west());
		CHECK(array[1] == Side<2>::south());
	}
	{
		auto array = Corner<2>::se().getExteriorSides();
		CHECK(array[0] == Side<2>::east());
		CHECK(array[1] == Side<2>::south());
	}
	{
		auto array = Corner<2>::nw().getExteriorSides();
		CHECK(array[0] == Side<2>::west());
		CHECK(array[1] == Side<2>::north());
	}
	{
		auto array = Corner<2>::ne().getExteriorSides();
		CHECK(array[0] == Side<2>::east());
		CHECK(array[1] == Side<2>::north());
	}
}
TEST_CASE("Corner<3> getExteriorSides is as expected", "[Corner]")
{
	{
		auto array = Corner<3>::bsw().getExteriorSides();
		CHECK(array[0] == Side<3>::west());
		CHECK(array[1] == Side<3>::south());
		CHECK(array[2] == Side<3>::bottom());
	}
	{
		auto array = Corner<3>::bse().getExteriorSides();
		CHECK(array[0] == Side<3>::east());
		CHECK(array[1] == Side<3>::south());
		CHECK(array[2] == Side<3>::bottom());
	}
	{
		auto array = Corner<3>::bnw().getExteriorSides();
		CHECK(array[0] == Side<3>::west());
		CHECK(array[1] == Side<3>::north());
		CHECK(array[2] == Side<3>::bottom());
	}
	{
		auto array = Corner<3>::bne().getExteriorSides();
		CHECK(array[0] == Side<3>::east());
		CHECK(array[1] == Side<3>::north());
		CHECK(array[2] == Side<3>::bottom());
	}
	{
		auto array = Corner<3>::tsw().getExteriorSides();
		CHECK(array[0] == Side<3>::west());
		CHECK(array[1] == Side<3>::south());
		CHECK(array[2] == Side<3>::top());
	}
	{
		auto array = Corner<3>::tse().getExteriorSides();
		CHECK(array[0] == Side<3>::east());
		CHECK(array[1] == Side<3>::south());
		CHECK(array[2] == Side<3>::top());
	}
	{
		auto array = Corner<3>::tnw().getExteriorSides();
		CHECK(array[0] == Side<3>::west());
		CHECK(array[1] == Side<3>::north());
		CHECK(array[2] == Side<3>::top());
	}
	{
		auto array = Corner<3>::tne().getExteriorSides();
		CHECK(array[0] == Side<3>::east());
		CHECK(array[1] == Side<3>::north());
		CHECK(array[2] == Side<3>::top());
	}
}
TEST_CASE("Corner<2> isOnSide is as expected", "[Corner]")
{
	CHECK(Corner<2>::sw().isOnSide(Side<2>::west()));
	CHECK_FALSE(Corner<2>::sw().isOnSide(Side<2>::east()));
	CHECK(Corner<2>::sw().isOnSide(Side<2>::south()));
	CHECK_FALSE(Corner<2>::sw().isOnSide(Side<2>::north()));

	CHECK_FALSE(Corner<2>::se().isOnSide(Side<2>::west()));
	CHECK(Corner<2>::se().isOnSide(Side<2>::east()));
	CHECK(Corner<2>::se().isOnSide(Side<2>::south()));
	CHECK_FALSE(Corner<2>::se().isOnSide(Side<2>::north()));

	CHECK(Corner<2>::nw().isOnSide(Side<2>::west()));
	CHECK_FALSE(Corner<2>::nw().isOnSide(Side<2>::east()));
	CHECK_FALSE(Corner<2>::nw().isOnSide(Side<2>::south()));
	CHECK(Corner<2>::nw().isOnSide(Side<2>::north()));

	CHECK_FALSE(Corner<2>::ne().isOnSide(Side<2>::west()));
	CHECK(Corner<2>::ne().isOnSide(Side<2>::east()));
	CHECK_FALSE(Corner<2>::ne().isOnSide(Side<2>::south()));
	CHECK(Corner<2>::ne().isOnSide(Side<2>::north()));
}
TEST_CASE("Corner<3> isOnSide is as expected", "[Corner]")
{
	CHECK(Corner<3>::bsw().isOnSide(Side<3>::west()));
	CHECK_FALSE(Corner<3>::bsw().isOnSide(Side<3>::east()));
	CHECK(Corner<3>::bsw().isOnSide(Side<3>::south()));
	CHECK_FALSE(Corner<3>::bsw().isOnSide(Side<3>::north()));
	CHECK(Corner<3>::bsw().isOnSide(Side<3>::bottom()));
	CHECK_FALSE(Corner<3>::bsw().isOnSide(Side<3>::top()));

	CHECK_FALSE(Corner<3>::bse().isOnSide(Side<3>::west()));
	CHECK(Corner<3>::bse().isOnSide(Side<3>::east()));
	CHECK(Corner<3>::bse().isOnSide(Side<3>::south()));
	CHECK_FALSE(Corner<3>::bse().isOnSide(Side<3>::north()));
	CHECK(Corner<3>::bse().isOnSide(Side<3>::bottom()));
	CHECK_FALSE(Corner<3>::bse().isOnSide(Side<3>::top()));

	CHECK(Corner<3>::bnw().isOnSide(Side<3>::west()));
	CHECK_FALSE(Corner<3>::bnw().isOnSide(Side<3>::east()));
	CHECK_FALSE(Corner<3>::bnw().isOnSide(Side<3>::south()));
	CHECK(Corner<3>::bnw().isOnSide(Side<3>::north()));
	CHECK(Corner<3>::bnw().isOnSide(Side<3>::bottom()));
	CHECK_FALSE(Corner<3>::bnw().isOnSide(Side<3>::top()));

	CHECK_FALSE(Corner<3>::bne().isOnSide(Side<3>::west()));
	CHECK(Corner<3>::bne().isOnSide(Side<3>::east()));
	CHECK_FALSE(Corner<3>::bne().isOnSide(Side<3>::south()));
	CHECK(Corner<3>::bne().isOnSide(Side<3>::north()));
	CHECK(Corner<3>::bne().isOnSide(Side<3>::bottom()));
	CHECK_FALSE(Corner<3>::bne().isOnSide(Side<3>::top()));

	CHECK(Corner<3>::tsw().isOnSide(Side<3>::west()));
	CHECK_FALSE(Corner<3>::tsw().isOnSide(Side<3>::east()));
	CHECK(Corner<3>::tsw().isOnSide(Side<3>::south()));
	CHECK_FALSE(Corner<3>::tsw().isOnSide(Side<3>::north()));
	CHECK_FALSE(Corner<3>::tsw().isOnSide(Side<3>::bottom()));
	CHECK(Corner<3>::tsw().isOnSide(Side<3>::top()));

	CHECK_FALSE(Corner<3>::tse().isOnSide(Side<3>::west()));
	CHECK(Corner<3>::tse().isOnSide(Side<3>::east()));
	CHECK(Corner<3>::tse().isOnSide(Side<3>::south()));
	CHECK_FALSE(Corner<3>::tse().isOnSide(Side<3>::north()));
	CHECK_FALSE(Corner<3>::tse().isOnSide(Side<3>::bottom()));
	CHECK(Corner<3>::tse().isOnSide(Side<3>::top()));

	CHECK(Corner<3>::tnw().isOnSide(Side<3>::west()));
	CHECK_FALSE(Corner<3>::tnw().isOnSide(Side<3>::east()));
	CHECK_FALSE(Corner<3>::tnw().isOnSide(Side<3>::south()));
	CHECK(Corner<3>::tnw().isOnSide(Side<3>::north()));
	CHECK_FALSE(Corner<3>::tnw().isOnSide(Side<3>::bottom()));
	CHECK(Corner<3>::tnw().isOnSide(Side<3>::top()));

	CHECK_FALSE(Corner<3>::tne().isOnSide(Side<3>::west()));
	CHECK(Corner<3>::tne().isOnSide(Side<3>::east()));
	CHECK_FALSE(Corner<3>::tne().isOnSide(Side<3>::south()));
	CHECK(Corner<3>::tne().isOnSide(Side<3>::north()));
	CHECK_FALSE(Corner<3>::tne().isOnSide(Side<3>::bottom()));
	CHECK(Corner<3>::tne().isOnSide(Side<3>::top()));
}
TEST_CASE("Corner<2> getValuesOnSide is as expected", "[Corner]")
{
	SECTION("Side<2>::west()")
	{
		std::array<Corner<2>, 2> values = Corner<2>::getValuesOnSide(Side<2>::west());
		CHECK(values[0] == Corner<2>::sw());
		CHECK(values[1] == Corner<2>::nw());
	}
	SECTION("Side<2>::east()")
	{
		std::array<Corner<2>, 2> values = Corner<2>::getValuesOnSide(Side<2>::east());
		CHECK(values[0] == Corner<2>::se());
		CHECK(values[1] == Corner<2>::ne());
	}
	SECTION("Side<2>::south()")
	{
		std::array<Corner<2>, 2> values = Corner<2>::getValuesOnSide(Side<2>::south());
		CHECK(values[0] == Corner<2>::sw());
		CHECK(values[1] == Corner<2>::se());
	}
	SECTION("Side<2>::north()")
	{
		std::array<Corner<2>, 2> values = Corner<2>::getValuesOnSide(Side<2>::north());
		CHECK(values[0] == Corner<2>::nw());
		CHECK(values[1] == Corner<2>::ne());
	}
}
TEST_CASE("Corner<3> getValuesOnSide is as expected", "[Corner]")
{
	SECTION("Side<3>::west()")
	{
		std::array<Corner<3>, 4> values = Corner<3>::getValuesOnSide(Side<3>::west());
		CHECK(values[0] == Corner<3>::bsw());
		CHECK(values[1] == Corner<3>::bnw());
		CHECK(values[2] == Corner<3>::tsw());
		CHECK(values[3] == Corner<3>::tnw());
	}
	SECTION("Side<3>::east()")
	{
		std::array<Corner<3>, 4> values = Corner<3>::getValuesOnSide(Side<3>::east());
		CHECK(values[0] == Corner<3>::bse());
		CHECK(values[1] == Corner<3>::bne());
		CHECK(values[2] == Corner<3>::tse());
		CHECK(values[3] == Corner<3>::tne());
	}
	SECTION("Side<3>::south()")
	{
		std::array<Corner<3>, 4> values = Corner<3>::getValuesOnSide(Side<3>::south());
		CHECK(values[0] == Corner<3>::bsw());
		CHECK(values[1] == Corner<3>::bse());
		CHECK(values[2] == Corner<3>::tsw());
		CHECK(values[3] == Corner<3>::tse());
	}
	SECTION("Side<3>::north()")
	{
		std::array<Corner<3>, 4> values = Corner<3>::getValuesOnSide(Side<3>::north());
		CHECK(values[0] == Corner<3>::bnw());
		CHECK(values[1] == Corner<3>::bne());
		CHECK(values[2] == Corner<3>::tnw());
		CHECK(values[3] == Corner<3>::tne());
	}
	SECTION("Side<3>::bottom()")
	{
		std::array<Corner<3>, 4> values = Corner<3>::getValuesOnSide(Side<3>::bottom());
		CHECK(values[0] == Corner<3>::bsw());
		CHECK(values[1] == Corner<3>::bse());
		CHECK(values[2] == Corner<3>::bnw());
		CHECK(values[3] == Corner<3>::bne());
	}
	SECTION("Side<3>::top()")
	{
		std::array<Corner<3>, 4> values = Corner<3>::getValuesOnSide(Side<3>::top());
		CHECK(values[0] == Corner<3>::tsw());
		CHECK(values[1] == Corner<3>::tse());
		CHECK(values[2] == Corner<3>::tnw());
		CHECK(values[3] == Corner<3>::tne());
	}
}
TEST_CASE("Corner<3> collapseOnAxis is as expected", "[Corner]")
{
	SECTION("Corner<3>::bsw()")
	{
		SECTION("x axis")
		{
			CHECK(Corner<3>::bsw().collapseOnAxis(0) == Corner<2>::sw());
		}
		SECTION("y axis")
		{
			CHECK(Corner<3>::bsw().collapseOnAxis(1) == Corner<2>::sw());
		}
		SECTION("z axis")
		{
			CHECK(Corner<3>::bsw().collapseOnAxis(2) == Corner<2>::sw());
		}
	}
	SECTION("Corner<3>::bse()")
	{
		SECTION("x axis")
		{
			CHECK(Corner<3>::bse().collapseOnAxis(0) == Corner<2>::sw());
		}
		SECTION("y axis")
		{
			CHECK(Corner<3>::bse().collapseOnAxis(1) == Corner<2>::se());
		}
		SECTION("z axis")
		{
			CHECK(Corner<3>::bse().collapseOnAxis(2) == Corner<2>::se());
		}
	}
	SECTION("Corner<3>::bnw()")
	{
		SECTION("x axis")
		{
			CHECK(Corner<3>::bnw().collapseOnAxis(0) == Corner<2>::se());
		}
		SECTION("y axis")
		{
			CHECK(Corner<3>::bnw().collapseOnAxis(1) == Corner<2>::sw());
		}
		SECTION("z axis")
		{
			CHECK(Corner<3>::bnw().collapseOnAxis(2) == Corner<2>::nw());
		}
	}
	SECTION("Corner<3>::bne()")
	{
		SECTION("x axis")
		{
			CHECK(Corner<3>::bne().collapseOnAxis(0) == Corner<2>::se());
		}
		SECTION("y axis")
		{
			CHECK(Corner<3>::bne().collapseOnAxis(1) == Corner<2>::se());
		}
		SECTION("z axis")
		{
			CHECK(Corner<3>::bne().collapseOnAxis(2) == Corner<2>::ne());
		}
	}
	SECTION("Corner<3>::tsw()")
	{
		SECTION("x axis")
		{
			CHECK(Corner<3>::tsw().collapseOnAxis(0) == Corner<2>::nw());
		}
		SECTION("y axis")
		{
			CHECK(Corner<3>::tsw().collapseOnAxis(1) == Corner<2>::nw());
		}
		SECTION("z axis")
		{
			CHECK(Corner<3>::tsw().collapseOnAxis(2) == Corner<2>::sw());
		}
	}
	SECTION("Corner<3>::tse()")
	{
		SECTION("x axis")
		{
			CHECK(Corner<3>::tse().collapseOnAxis(0) == Corner<2>::nw());
		}
		SECTION("y axis")
		{
			CHECK(Corner<3>::tse().collapseOnAxis(1) == Corner<2>::ne());
		}
		SECTION("z axis")
		{
			CHECK(Corner<3>::tse().collapseOnAxis(2) == Corner<2>::se());
		}
	}
	SECTION("Corner<3>::tnw()")
	{
		SECTION("x axis")
		{
			CHECK(Corner<3>::tnw().collapseOnAxis(0) == Corner<2>::ne());
		}
		SECTION("y axis")
		{
			CHECK(Corner<3>::tnw().collapseOnAxis(1) == Corner<2>::nw());
		}
		SECTION("z axis")
		{
			CHECK(Corner<3>::tnw().collapseOnAxis(2) == Corner<2>::nw());
		}
	}
	SECTION("Corner<3>::tne()")
	{
		SECTION("x axis")
		{
			CHECK(Corner<3>::tne().collapseOnAxis(0) == Corner<2>::ne());
		}
		SECTION("y axis")
		{
			CHECK(Corner<3>::tne().collapseOnAxis(1) == Corner<2>::ne());
		}
		SECTION("z axis")
		{
			CHECK(Corner<3>::tne().collapseOnAxis(2) == Corner<2>::ne());
		}
	}
}
TEST_CASE("Corner<1> ==", "[Corner]")
{
	CHECK(Corner<1>::null() == Corner<1>::null());
}
TEST_CASE("Corner<2> ==", "[Corner]")
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
TEST_CASE("Corner<3> ==", "[Corner]")
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
TEST_CASE("Corner<1> !=", "[Corner]")
{
	CHECK_FALSE(Corner<1>::null() != Corner<1>::null());
}
TEST_CASE("Corner<2> !=", "[Corner]")
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
TEST_CASE("Corner<3> !=", "[Corner]")
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
TEST_CASE("Corner<1> <", "[Corner]")
{
	CHECK_FALSE(Corner<1>::null() < Corner<1>::null());
}
TEST_CASE("Corner<2> <", "[Corner]")
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
TEST_CASE("Corner<3> <", "[Corner]")
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
TEST_CASE("Test ostream for Corner<1>", "[Corner]")
{
	stringstream ss;
	ss << Corner<1>::null();
	CHECK(ss.str() == "Corner<1>::null()");
	ss.str("");
	ss << Corner<1>(13);
	CHECK(ss.str() == "Corner<1> invalid value: 13");
}
TEST_CASE("Test ostream for Corner<2>", "[Corner]")
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
TEST_CASE("Test ostream for Corner<3>", "[Corner]")
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
TEST_CASE("Test iterator for Corner<1>", "[Corner]")
{
	auto iter = Corner<1>::getValues().begin();
	CHECK(*iter == Corner<1>::null());
	CHECK(iter == Corner<1>::getValues().end());
}
TEST_CASE("Test iterator for Corner<2>", "[Corner]")
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
TEST_CASE("Test iterator for Corner<3>", "[Corner]")
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
TEST_CASE("Test from_json for Corner<1>", "[Corner]")
{
	nlohmann::json j;
	j["null"] = nullptr;
	CHECK(j["null"].get<Corner<1>>() == Corner<1>::null());
}
TEST_CASE("Test from_json for Corner<2>", "[Corner]")
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
TEST_CASE("Test from_json for Corner<3>", "[Corner]")
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
TEST_CASE("Test to_json for Corner<1>", "[Corner]")
{
	nlohmann::json j;
	j["null"] = Corner<1>::null();
	CHECK(j["null"] == nullptr);
}
TEST_CASE("Test to_json for Corner<2>", "[Corner]")
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
TEST_CASE("Test to_json for Corner<3>", "[Corner]")
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