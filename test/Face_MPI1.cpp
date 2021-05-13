#include <ThunderEgg/Face.h>
#include <sstream>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Test number_of for Face<1,0>", "[Face]")
{
	size_t faces = Face<1, 0>::number_of;
	CHECK(faces == 2);
}

TEST_CASE("Test number_of for Face<2,0>", "[Face]")
{
	size_t faces = Face<2, 0>::number_of;
	CHECK(faces == 4);
}

TEST_CASE("Test number_of for Face<2,1>", "[Face]")
{
	size_t faces = Face<2, 1>::number_of;
	CHECK(faces == 4);
}

TEST_CASE("Test number_of for Face<3,0>", "[Face]")
{
	size_t faces = Face<3, 0>::number_of;
	CHECK(faces == 8);
}

TEST_CASE("Test number_of for Face<3,1>", "[Face]")
{
	size_t faces = Face<3, 1>::number_of;
	CHECK(faces == 12);
}

TEST_CASE("Test number_of for Face<3,2>", "[Face]")
{
	size_t faces = Face<3, 2>::number_of;
	CHECK(faces == 6);
}

TEST_CASE("getIndex returns value passed to constructor Face<3,2>", "[Face]")
{
	unsigned char val = GENERATE(1, 2, 3);
	Face<3, 2>    face(val);
	CHECK(face.getIndex() == val);
}