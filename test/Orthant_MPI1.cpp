#include <Thunderegg/Orthant.h>

#include "catch.hpp"
#include <sstream>

using namespace std;
using namespace Thunderegg;
TEST_CASE("Orthant<1> unsigned char constructor works", "[Octant]")
{
	Orthant<1> o(13);
	CHECK(o.getIndex() == 13);
}
TEST_CASE("Orthant<2> unsigned char constructor works", "[Octant]")
{
	Orthant<2> o(13);
	CHECK(o.getIndex() == 13);
}
TEST_CASE("Orthant<3> unsigned char constructor works", "[Octant]")
{
	Orthant<3> o(13);
	CHECK(o.getIndex() == 13);
}
TEST_CASE("Orthant<1> Default constructor works", "[Octant]")
{
	Orthant<1> o;
	CHECK(o == Orthant<1>::null());
}
TEST_CASE("Orthant<2> Default constructor works", "[Octant]")
{
	Orthant<2> o;
	CHECK(o == Orthant<2>::null());
}
TEST_CASE("Orthant<3> Default constructor works", "[Octant]")
{
	Orthant<3> o;
	CHECK(o == Orthant<3>::null());
}
TEST_CASE("Orthant<1> named constructors give expected index values", "[Octant]")
{
	CHECK(Orthant<1>::lower().getIndex() == 0);
	CHECK(Orthant<1>::upper().getIndex() == 1);
	CHECK(Orthant<1>::null().getIndex() == 2);
}
TEST_CASE("Orthant<2> named constructors give expected index values", "[Octant]")
{
	CHECK(Orthant<2>::sw().getIndex() == 0);
	CHECK(Orthant<2>::se().getIndex() == 1);
	CHECK(Orthant<2>::nw().getIndex() == 2);
	CHECK(Orthant<2>::ne().getIndex() == 3);
	CHECK(Orthant<2>::null().getIndex() == 4);
}
TEST_CASE("Orthant<3> named constructors give expected index values", "[Octant]")
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
TEST_CASE("Orthant<3> getInteriorNbrOnSide() is as expected", "[Octant]")
{
	CHECK(Orthant<3>::bsw().getInteriorNbrOnSide(Side<3>::east()) == Orthant<3>::bse());
	CHECK(Orthant<3>::bsw().getInteriorNbrOnSide(Side<3>::north()) == Orthant<3>::bnw());
	CHECK(Orthant<3>::bsw().getInteriorNbrOnSide(Side<3>::top()) == Orthant<3>::tsw());

	CHECK(Orthant<3>::bse().getInteriorNbrOnSide(Side<3>::west()) == Orthant<3>::bsw());
	CHECK(Orthant<3>::bse().getInteriorNbrOnSide(Side<3>::north()) == Orthant<3>::bne());
	CHECK(Orthant<3>::bse().getInteriorNbrOnSide(Side<3>::top()) == Orthant<3>::tse());

	CHECK(Orthant<3>::bnw().getInteriorNbrOnSide(Side<3>::east()) == Orthant<3>::bne());
	CHECK(Orthant<3>::bnw().getInteriorNbrOnSide(Side<3>::south()) == Orthant<3>::bsw());
	CHECK(Orthant<3>::bnw().getInteriorNbrOnSide(Side<3>::top()) == Orthant<3>::tnw());

	CHECK(Orthant<3>::bne().getInteriorNbrOnSide(Side<3>::west()) == Orthant<3>::bnw());
	CHECK(Orthant<3>::bne().getInteriorNbrOnSide(Side<3>::south()) == Orthant<3>::bse());
	CHECK(Orthant<3>::bne().getInteriorNbrOnSide(Side<3>::top()) == Orthant<3>::tne());

	CHECK(Orthant<3>::tsw().getInteriorNbrOnSide(Side<3>::east()) == Orthant<3>::tse());
	CHECK(Orthant<3>::tsw().getInteriorNbrOnSide(Side<3>::north()) == Orthant<3>::tnw());
	CHECK(Orthant<3>::tsw().getInteriorNbrOnSide(Side<3>::bottom()) == Orthant<3>::bsw());

	CHECK(Orthant<3>::tse().getInteriorNbrOnSide(Side<3>::west()) == Orthant<3>::tsw());
	CHECK(Orthant<3>::tse().getInteriorNbrOnSide(Side<3>::north()) == Orthant<3>::tne());
	CHECK(Orthant<3>::tse().getInteriorNbrOnSide(Side<3>::bottom()) == Orthant<3>::bse());

	CHECK(Orthant<3>::tnw().getInteriorNbrOnSide(Side<3>::east()) == Orthant<3>::tne());
	CHECK(Orthant<3>::tnw().getInteriorNbrOnSide(Side<3>::south()) == Orthant<3>::tsw());
	CHECK(Orthant<3>::tnw().getInteriorNbrOnSide(Side<3>::bottom()) == Orthant<3>::bnw());

	CHECK(Orthant<3>::tne().getInteriorNbrOnSide(Side<3>::west()) == Orthant<3>::tnw());
	CHECK(Orthant<3>::tne().getInteriorNbrOnSide(Side<3>::south()) == Orthant<3>::tse());
	CHECK(Orthant<3>::tne().getInteriorNbrOnSide(Side<3>::bottom()) == Orthant<3>::bne());
}
TEST_CASE("Orthant<3> getExteriorNbrOnSide() is as expected", "[Octant]")
{
	CHECK(Orthant<3>::bsw().getExteriorNbrOnSide(Side<3>::west()) == Orthant<3>::bse());
	CHECK(Orthant<3>::bsw().getExteriorNbrOnSide(Side<3>::south()) == Orthant<3>::bnw());
	CHECK(Orthant<3>::bsw().getExteriorNbrOnSide(Side<3>::bottom()) == Orthant<3>::tsw());

	CHECK(Orthant<3>::bse().getExteriorNbrOnSide(Side<3>::east()) == Orthant<3>::bsw());
	CHECK(Orthant<3>::bse().getExteriorNbrOnSide(Side<3>::south()) == Orthant<3>::bne());
	CHECK(Orthant<3>::bse().getExteriorNbrOnSide(Side<3>::bottom()) == Orthant<3>::tse());

	CHECK(Orthant<3>::bnw().getExteriorNbrOnSide(Side<3>::west()) == Orthant<3>::bne());
	CHECK(Orthant<3>::bnw().getExteriorNbrOnSide(Side<3>::north()) == Orthant<3>::bsw());
	CHECK(Orthant<3>::bnw().getExteriorNbrOnSide(Side<3>::bottom()) == Orthant<3>::tnw());

	CHECK(Orthant<3>::bne().getExteriorNbrOnSide(Side<3>::east()) == Orthant<3>::bnw());
	CHECK(Orthant<3>::bne().getExteriorNbrOnSide(Side<3>::north()) == Orthant<3>::bse());
	CHECK(Orthant<3>::bne().getExteriorNbrOnSide(Side<3>::bottom()) == Orthant<3>::tne());

	CHECK(Orthant<3>::tsw().getExteriorNbrOnSide(Side<3>::west()) == Orthant<3>::tse());
	CHECK(Orthant<3>::tsw().getExteriorNbrOnSide(Side<3>::south()) == Orthant<3>::tnw());
	CHECK(Orthant<3>::tsw().getExteriorNbrOnSide(Side<3>::top()) == Orthant<3>::bsw());

	CHECK(Orthant<3>::tse().getExteriorNbrOnSide(Side<3>::east()) == Orthant<3>::tsw());
	CHECK(Orthant<3>::tse().getExteriorNbrOnSide(Side<3>::south()) == Orthant<3>::tne());
	CHECK(Orthant<3>::tse().getExteriorNbrOnSide(Side<3>::top()) == Orthant<3>::bse());

	CHECK(Orthant<3>::tnw().getExteriorNbrOnSide(Side<3>::west()) == Orthant<3>::tne());
	CHECK(Orthant<3>::tnw().getExteriorNbrOnSide(Side<3>::north()) == Orthant<3>::tsw());
	CHECK(Orthant<3>::tnw().getExteriorNbrOnSide(Side<3>::top()) == Orthant<3>::bnw());

	CHECK(Orthant<3>::tne().getExteriorNbrOnSide(Side<3>::east()) == Orthant<3>::tnw());
	CHECK(Orthant<3>::tne().getExteriorNbrOnSide(Side<3>::north()) == Orthant<3>::tse());
	CHECK(Orthant<3>::tne().getExteriorNbrOnSide(Side<3>::top()) == Orthant<3>::bne());
}
TEST_CASE("Orthant<1> getInteriorSides() is as expected", "[Octant]")
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
TEST_CASE("Orthant<2> getInteriorSides() is as expected", "[Octant]")
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
TEST_CASE("Orthant<3> getInteriorSides() is as expected", "[Octant]")
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
TEST_CASE("Orthant<1> getExteriorSides() is as expected", "[Octant]")
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
TEST_CASE("Orthant<2> getExteriorSides() is as expected", "[Octant]")
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
TEST_CASE("Orthant<3> getExteriorSides() is as expected", "[Octant]")
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
TEST_CASE("Orthant<1> isOnSide() is as expected", "[Octant]")
{
	CHECK(Orthant<1>::lower().isOnSide(Side<1>::west()));
	CHECK_FALSE(Orthant<1>::lower().isOnSide(Side<1>::east()));

	CHECK_FALSE(Orthant<1>::upper().isOnSide(Side<1>::west()));
	CHECK(Orthant<1>::upper().isOnSide(Side<1>::east()));
}
TEST_CASE("Orthant<2> isOnSide() is as expected", "[Octant]")
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
TEST_CASE("Orthant<3> isOnSide() is as expected", "[Octant]")
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
TEST_CASE("Orthant<3> getValues() is as expected", "[Octant]")
{
	std::array<Orthant<3>, 8> values = Orthant<3>::getValues();
	CHECK(values[0] == Orthant<3>::bsw());
	CHECK(values[1] == Orthant<3>::bse());
	CHECK(values[2] == Orthant<3>::bnw());
	CHECK(values[3] == Orthant<3>::bne());
	CHECK(values[4] == Orthant<3>::tsw());
	CHECK(values[5] == Orthant<3>::tse());
	CHECK(values[6] == Orthant<3>::tnw());
	CHECK(values[7] == Orthant<3>::tne());
}
TEST_CASE("Orthant<1> getValuesOnSide() is as expected", "[Octant]")
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
TEST_CASE("Orthant<2> getValuesOnSide() is as expected", "[Octant]")
{
	SECTION("Side<2>::west()")
	{
		std::array<Orthant<2>, 2> values = Orthant<2>::getValuesOnSide(Side<2>::west());
		CHECK(values[0] == Orthant<2>::sw());
		CHECK(values[1] == Orthant<2>::nw());
	}
	SECTION("Side<2>::east()")
	{
		std::array<Orthant<2>, 2> values = Orthant<2>::getValuesOnSide(Side<2>::east());
		CHECK(values[0] == Orthant<2>::se());
		CHECK(values[1] == Orthant<2>::ne());
	}
	SECTION("Side<2>::south()")
	{
		std::array<Orthant<2>, 2> values = Orthant<2>::getValuesOnSide(Side<2>::south());
		CHECK(values[0] == Orthant<2>::sw());
		CHECK(values[1] == Orthant<2>::se());
	}
	SECTION("Side<2>::north()")
	{
		std::array<Orthant<2>, 2> values = Orthant<2>::getValuesOnSide(Side<2>::north());
		CHECK(values[0] == Orthant<2>::nw());
		CHECK(values[1] == Orthant<2>::ne());
	}
}
TEST_CASE("Orthant<3> getValuesOnSide() is as expected", "[Octant]")
{
	SECTION("Side<3>::west()")
	{
		std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::west());
		CHECK(values[0] == Orthant<3>::bsw());
		CHECK(values[1] == Orthant<3>::bnw());
		CHECK(values[2] == Orthant<3>::tsw());
		CHECK(values[3] == Orthant<3>::tnw());
	}
	SECTION("Side<3>::east()")
	{
		std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::east());
		CHECK(values[0] == Orthant<3>::bse());
		CHECK(values[1] == Orthant<3>::bne());
		CHECK(values[2] == Orthant<3>::tse());
		CHECK(values[3] == Orthant<3>::tne());
	}
	SECTION("Side<3>::south()")
	{
		std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::south());
		CHECK(values[0] == Orthant<3>::bsw());
		CHECK(values[1] == Orthant<3>::bse());
		CHECK(values[2] == Orthant<3>::tsw());
		CHECK(values[3] == Orthant<3>::tse());
	}
	SECTION("Side<3>::north()")
	{
		std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::north());
		CHECK(values[0] == Orthant<3>::bnw());
		CHECK(values[1] == Orthant<3>::bne());
		CHECK(values[2] == Orthant<3>::tnw());
		CHECK(values[3] == Orthant<3>::tne());
	}
	SECTION("Side<3>::bottom()")
	{
		std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::bottom());
		CHECK(values[0] == Orthant<3>::bsw());
		CHECK(values[1] == Orthant<3>::bse());
		CHECK(values[2] == Orthant<3>::bnw());
		CHECK(values[3] == Orthant<3>::bne());
	}
	SECTION("Side<3>::top()")
	{
		std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::top());
		CHECK(values[0] == Orthant<3>::tsw());
		CHECK(values[1] == Orthant<3>::tse());
		CHECK(values[2] == Orthant<3>::tnw());
		CHECK(values[3] == Orthant<3>::tne());
	}
}
TEST_CASE("Orthant<1> collapseOnAxis() is as expected", "[Octant]")
{
	SECTION("Orthant<1>::lower()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<1>::lower().collapseOnAxis(0) == Orthant<0>::null());
		}
	}
	SECTION("Orthant<1>::upper()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<1>::upper().collapseOnAxis(0) == Orthant<0>::null());
		}
	}
}
TEST_CASE("Orthant<2> collapseOnAxis() is as expected", "[Octant]")
{
	SECTION("Orthant<2>::sw()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<2>::sw().collapseOnAxis(0) == Orthant<1>::lower());
		}
		SECTION("y axis")
		{
			CHECK(Orthant<2>::sw().collapseOnAxis(1) == Orthant<1>::lower());
		}
	}
	SECTION("Orthant<2>::se()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<2>::se().collapseOnAxis(0) == Orthant<1>::lower());
		}
		SECTION("y axis")
		{
			CHECK(Orthant<2>::se().collapseOnAxis(1) == Orthant<1>::upper());
		}
	}
	SECTION("Orthant<2>::nw()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<2>::nw().collapseOnAxis(0) == Orthant<1>::upper());
		}
		SECTION("y axis")
		{
			CHECK(Orthant<2>::nw().collapseOnAxis(1) == Orthant<1>::lower());
		}
	}
	SECTION("Orthant<2>::ne()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<2>::ne().collapseOnAxis(0) == Orthant<1>::upper());
		}
		SECTION("y axis")
		{
			CHECK(Orthant<2>::ne().collapseOnAxis(1) == Orthant<1>::upper());
		}
	}
}
TEST_CASE("Orthant<3> collapseOnAxis() is as expected", "[Octant]")
{
	SECTION("Orthant<3>::bsw()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<3>::bsw().collapseOnAxis(0) == Orthant<2>::sw());
		}
		SECTION("y axis")
		{
			CHECK(Orthant<3>::bsw().collapseOnAxis(1) == Orthant<2>::sw());
		}
		SECTION("z axis")
		{
			CHECK(Orthant<3>::bsw().collapseOnAxis(2) == Orthant<2>::sw());
		}
	}
	SECTION("Orthant<3>::bse()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<3>::bse().collapseOnAxis(0) == Orthant<2>::sw());
		}
		SECTION("y axis")
		{
			CHECK(Orthant<3>::bse().collapseOnAxis(1) == Orthant<2>::se());
		}
		SECTION("z axis")
		{
			CHECK(Orthant<3>::bse().collapseOnAxis(2) == Orthant<2>::se());
		}
	}
	SECTION("Orthant<3>::bnw()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<3>::bnw().collapseOnAxis(0) == Orthant<2>::se());
		}
		SECTION("y axis")
		{
			CHECK(Orthant<3>::bnw().collapseOnAxis(1) == Orthant<2>::sw());
		}
		SECTION("z axis")
		{
			CHECK(Orthant<3>::bnw().collapseOnAxis(2) == Orthant<2>::nw());
		}
	}
	SECTION("Orthant<3>::bne()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<3>::bne().collapseOnAxis(0) == Orthant<2>::se());
		}
		SECTION("y axis")
		{
			CHECK(Orthant<3>::bne().collapseOnAxis(1) == Orthant<2>::se());
		}
		SECTION("z axis")
		{
			CHECK(Orthant<3>::bne().collapseOnAxis(2) == Orthant<2>::ne());
		}
	}
	SECTION("Orthant<3>::tsw()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<3>::tsw().collapseOnAxis(0) == Orthant<2>::nw());
		}
		SECTION("y axis")
		{
			CHECK(Orthant<3>::tsw().collapseOnAxis(1) == Orthant<2>::nw());
		}
		SECTION("z axis")
		{
			CHECK(Orthant<3>::tsw().collapseOnAxis(2) == Orthant<2>::sw());
		}
	}
	SECTION("Orthant<3>::tse()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<3>::tse().collapseOnAxis(0) == Orthant<2>::nw());
		}
		SECTION("y axis")
		{
			CHECK(Orthant<3>::tse().collapseOnAxis(1) == Orthant<2>::ne());
		}
		SECTION("z axis")
		{
			CHECK(Orthant<3>::tse().collapseOnAxis(2) == Orthant<2>::se());
		}
	}
	SECTION("Orthant<3>::tnw()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<3>::tnw().collapseOnAxis(0) == Orthant<2>::ne());
		}
		SECTION("y axis")
		{
			CHECK(Orthant<3>::tnw().collapseOnAxis(1) == Orthant<2>::nw());
		}
		SECTION("z axis")
		{
			CHECK(Orthant<3>::tnw().collapseOnAxis(2) == Orthant<2>::nw());
		}
	}
	SECTION("Orthant<3>::tne()")
	{
		SECTION("x axis")
		{
			CHECK(Orthant<3>::tne().collapseOnAxis(0) == Orthant<2>::ne());
		}
		SECTION("y axis")
		{
			CHECK(Orthant<3>::tne().collapseOnAxis(1) == Orthant<2>::ne());
		}
		SECTION("z axis")
		{
			CHECK(Orthant<3>::tne().collapseOnAxis(2) == Orthant<2>::ne());
		}
	}
}
TEST_CASE("Orthant<1> ==", "[Octant]")
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
TEST_CASE("Orthant<2> ==", "[Octant]")
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
TEST_CASE("Orthant<3> ==", "[Octant]")
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
TEST_CASE("Orthant<1> !=", "[Octant]")
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
TEST_CASE("Orthant<2> !=", "[Octant]")
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
TEST_CASE("Orthant<3> !=", "[Octant]")
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
TEST_CASE("Orthant<1> <", "[Octant]")
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
TEST_CASE("Orthant<2> <", "[Octant]")
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
TEST_CASE("Orthant<3> <", "[Octant]")
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