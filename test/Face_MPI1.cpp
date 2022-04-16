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

TEST_CASE("Test sum_of_faces for Face<1,0>", "[Face]")
{
  size_t sum = Face<1, 0>::sum_of_faces;
  CHECK(sum == 0);
}

TEST_CASE("Test sum_of_faces for Face<1,1>", "[Face]")
{
  size_t sum = Face<1, 1>::sum_of_faces;
  CHECK(sum == 2);
}

TEST_CASE("Test sum_of_faces for Face<2,0>", "[Face]")
{
  size_t sum = Face<2, 0>::sum_of_faces;
  CHECK(sum == 0);
}

TEST_CASE("Test sum_of_faces for Face<2,1>", "[Face]")
{
  size_t sum = Face<2, 1>::sum_of_faces;
  CHECK(sum == 4);
}

TEST_CASE("Test sum_of_faces for Face<2,2>", "[Face]")
{
  size_t sum = Face<2, 2>::sum_of_faces;
  CHECK(sum == 8);
}

TEST_CASE("Test sum_of_faces for Face<3,0>", "[Face]")
{
  size_t sum = Face<3, 0>::sum_of_faces;
  CHECK(sum == 0);
}

TEST_CASE("Test sum_of_faces for Face<3,1>", "[Face]")
{
  size_t sum = Face<3, 1>::sum_of_faces;
  CHECK(sum == 8);
}
TEST_CASE("Test sum_of_faces for Face<3,2>", "[Face]")
{
  size_t sum = Face<3, 2>::sum_of_faces;
  CHECK(sum == 20);
}
TEST_CASE("Test sum_of_faces for Face<3,3>", "[Face]")
{
  size_t sum = Face<3, 3>::sum_of_faces;
  CHECK(sum == 26);
}

TEST_CASE("getIndex returns value passed to constructor Face<3,2>", "[Face]")
{
  for (unsigned char val : { 1, 2, 3 }) {
    Face<3, 2> face(val);
    CHECK(face.getIndex() == val);
  }
}
