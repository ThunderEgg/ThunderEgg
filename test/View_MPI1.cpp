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
#include <ThunderEgg/View.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
using namespace std;
using namespace ThunderEgg;
TEST_CASE("View default constructor", "[View]")
{
	View<double, 2> v;

	for (int i = 0; i < 2; i++) {
		CHECK(v.getStrides()[i] == 0);
		CHECK(v.getStart()[i] == 0);
		CHECK(v.getEnd()[i] == -1);
		CHECK(v.getGhostStart()[i] == 0);
		CHECK(v.getGhostEnd()[i] == -1);
	}
	if constexpr (ENABLE_DEBUG) {
		CHECK_THROWS_AS((v[{0, 0}]), RuntimeError);
		CHECK_THROWS_AS(v(0, 0), RuntimeError);
	}
}
TEST_CASE("View constructor", "[View]")
{
	auto x_stride      = GENERATE(1, 2);
	auto y_stride      = GENERATE(1, 2);
	auto x_ghost_start = GENERATE(-1, 1);
	auto x_start       = GENERATE(2, 3);
	auto x_end         = GENERATE(3, 4);
	auto x_ghost_end   = GENERATE(5, 6);
	auto y_ghost_start = GENERATE(-1, 1);
	auto y_start       = GENERATE(2, 3);
	auto y_end         = GENERATE(3, 4);
	auto y_ghost_end   = GENERATE(5, 6);

	double data;

	View<double, 2> v(&data, {x_stride, y_stride}, {x_ghost_start, y_ghost_start}, {x_start, y_start}, {x_end, y_end}, {x_ghost_end, y_ghost_end});

	CHECK(v.getStrides()[0] == x_stride);
	CHECK(v.getStrides()[1] == y_stride);
	CHECK(v.getGhostStart()[0] == x_ghost_start);
	CHECK(v.getGhostStart()[1] == y_ghost_start);
	CHECK(v.getStart()[0] == x_start);
	CHECK(v.getStart()[1] == y_start);
	CHECK(v.getEnd()[0] == x_end);
	CHECK(v.getEnd()[1] == y_end);
	CHECK(v.getGhostEnd()[0] == x_ghost_end);
	CHECK(v.getGhostEnd()[1] == y_ghost_end);
	CHECK(&v[v.getGhostStart()] == &data);
}
TEST_CASE("View squarebracket operator", "[View]")
{
	auto x_stride      = GENERATE(1, 2);
	auto y_stride      = GENERATE(4, 5);
	auto x_ghost_start = GENERATE(-1, 1);
	auto x_start       = GENERATE(2, 3);
	auto x_end         = GENERATE(3, 4);
	auto x_ghost_end   = GENERATE(5, 6);
	auto y_ghost_start = GENERATE(-1, 1);
	auto y_start       = GENERATE(2, 3);
	auto y_end         = GENERATE(3, 4);
	auto y_ghost_end   = GENERATE(5, 6);

	double data[100];
	for (int i = 0; i < 100; i++) {
		data[i] = 0;
	}

	View<double, 2> v(data, {x_stride, y_stride}, {x_ghost_start, y_ghost_start}, {x_start, y_start}, {x_end, y_end}, {x_ghost_end, y_ghost_end});

	double *origin = data - (x_ghost_start * x_stride + y_ghost_start * y_stride);
	for (int yi = y_ghost_start - 1; yi <= y_ghost_end + 1; yi++) {
		for (int xi = x_ghost_start - 1; xi <= x_ghost_end + 1; xi++) {
			if (xi < x_ghost_start || xi > x_ghost_end || yi < y_ghost_start || yi > y_ghost_end) {
				//oob coord
				if constexpr (ENABLE_DEBUG) {
					CHECK_THROWS_AS((v[{xi, yi}]), RuntimeError);
				}
			} else {
				CHECK(&v[{xi, yi}] == origin + xi * x_stride + yi * y_stride);
			}
		}
	}
}
TEST_CASE("View parens operator", "[View]")
{
	auto x_stride      = GENERATE(1, 2);
	auto y_stride      = GENERATE(4, 5);
	auto x_ghost_start = GENERATE(-1, 1);
	auto x_start       = GENERATE(2, 3);
	auto x_end         = GENERATE(3, 4);
	auto x_ghost_end   = GENERATE(5, 6);
	auto y_ghost_start = GENERATE(-1, 1);
	auto y_start       = GENERATE(2, 3);
	auto y_end         = GENERATE(3, 4);
	auto y_ghost_end   = GENERATE(5, 6);

	double data[100];
	for (int i = 0; i < 100; i++) {
		data[i] = 0;
	}

	View<double, 2> v(data, {x_stride, y_stride}, {x_ghost_start, y_ghost_start}, {x_start, y_start}, {x_end, y_end}, {x_ghost_end, y_ghost_end});

	double *origin = data - (x_ghost_start * x_stride + y_ghost_start * y_stride);
	for (int yi = y_ghost_start - 1; yi <= y_ghost_end + 1; yi++) {
		for (int xi = x_ghost_start - 1; xi <= x_ghost_end + 1; xi++) {
			if (xi < x_ghost_start || xi > x_ghost_end || yi < y_ghost_start || yi > y_ghost_end) {
				//oob coord
				if constexpr (ENABLE_DEBUG) {
					CHECK_THROWS_AS((v(xi, yi)), RuntimeError);
				}
			} else {
				CHECK(&v(xi, yi) == origin + xi * x_stride + yi * y_stride);
			}
		}
	}
}
TEST_CASE("View set", "[View]")
{
	auto x_stride      = GENERATE(1, 2);
	auto y_stride      = GENERATE(4, 5);
	auto x_ghost_start = GENERATE(-1, 1);
	auto x_start       = GENERATE(2, 3);
	auto x_end         = GENERATE(3, 4);
	auto x_ghost_end   = GENERATE(5, 6);
	auto y_ghost_start = GENERATE(-1, 1);
	auto y_start       = GENERATE(2, 3);
	auto y_end         = GENERATE(3, 4);
	auto y_ghost_end   = GENERATE(5, 6);

	double data[1000];
	for (int i = 0; i < 1000; i++) {
		data[i] = 0;
	}

	View<double, 2> v(data, {x_stride, y_stride}, {x_ghost_start, y_ghost_start}, {x_start, y_start}, {x_end, y_end}, {x_ghost_end, y_ghost_end});

	double value = 0;
	for (int yi = y_ghost_start - 1; yi <= y_ghost_end + 1; yi++) {
		for (int xi = x_ghost_start - 1; xi <= x_ghost_end + 1; xi++) {
			if (xi < x_ghost_start || xi > x_ghost_end || yi < y_ghost_start || yi > y_ghost_end) {
				//oob coord
				if constexpr (ENABLE_DEBUG) {
					CHECK_THROWS_AS(v.set({xi, yi}, value), RuntimeError);
				}
			} else {
				v.set({xi, yi}, value);
				CHECK(v[{xi, yi}] == value);
			}
			value++;
		}
	}
}
TEST_CASE("View set const", "[View]")
{
	auto x_stride      = GENERATE(1, 2);
	auto y_stride      = GENERATE(4, 5);
	auto x_ghost_start = GENERATE(-1, 1);
	auto x_start       = GENERATE(2, 3);
	auto x_end         = GENERATE(3, 4);
	auto x_ghost_end   = GENERATE(5, 6);
	auto y_ghost_start = GENERATE(-1, 1);
	auto y_start       = GENERATE(2, 3);
	auto y_end         = GENERATE(3, 4);
	auto y_ghost_end   = GENERATE(5, 6);

	double data[1000];
	for (int i = 0; i < 1000; i++) {
		data[i] = 0;
	}

	View<const double, 2> v(data, {x_stride, y_stride}, {x_ghost_start, y_ghost_start}, {x_start, y_start}, {x_end, y_end}, {x_ghost_end, y_ghost_end});

	double value = 0;
	for (int yi = y_ghost_start - 1; yi <= y_ghost_end + 1; yi++) {
		for (int xi = x_ghost_start - 1; xi <= x_ghost_end + 1; xi++) {
			if (xi >= x_start && xi <= x_end && yi >= y_start && yi <= y_end) {
				//intertior coord
				if constexpr (ENABLE_DEBUG) {
					CHECK_THROWS_AS(v.set({xi, yi}, value), RuntimeError);
				}
			} else if (xi < x_ghost_start || xi > x_ghost_end || yi < y_ghost_start || yi > y_ghost_end) {
				//oob coord
				if constexpr (ENABLE_DEBUG) {
					CHECK_THROWS_AS(v.set({xi, yi}, value), RuntimeError);
				}
			} else {
				v.set({xi, yi}, value);
				CHECK(v[{xi, yi}] == value);
			}
			value++;
		}
	}
}
TEST_CASE("View implicit conversion to const type", "[View]")
{
	auto x_stride      = GENERATE(1, 2);
	auto y_stride      = GENERATE(4, 5);
	auto x_ghost_start = GENERATE(-1, 1);
	auto x_start       = GENERATE(2, 3);
	auto x_end         = GENERATE(3, 4);
	auto x_ghost_end   = GENERATE(5, 6);
	auto y_ghost_start = GENERATE(-1, 1);
	auto y_start       = GENERATE(2, 3);
	auto y_end         = GENERATE(3, 4);
	auto y_ghost_end   = GENERATE(5, 6);

	double data[1000];
	for (int i = 0; i < 1000; i++) {
		data[i] = 0;
	}

	View<double, 2>       v(data, {x_stride, y_stride}, {x_ghost_start, y_ghost_start}, {x_start, y_start}, {x_end, y_end}, {x_ghost_end, y_ghost_end});
	View<const double, 2> vc = v;

	CHECK(vc.getGhostStart() == v.getGhostStart());
	CHECK(vc.getStart() == v.getStart());
	CHECK(vc.getEnd() == v.getEnd());
	CHECK(vc.getGhostEnd() == v.getGhostEnd());
	CHECK(vc.getStrides() == v.getStrides());
	CHECK(&vc[vc.getGhostStart()] == &v[v.getGhostStart()]);
}