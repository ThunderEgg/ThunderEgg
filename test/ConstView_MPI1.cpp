#include <ThunderEgg/ConstView.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
using namespace std;
using namespace ThunderEgg;
TEST_CASE("ConstView default constructor", "[ConstView]")
{
	ConstView<2> v;

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
TEST_CASE("ConstView constructor", "[ConstView]")
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

	ConstView<2> v(&data, {x_stride, y_stride}, {x_ghost_start, y_ghost_start}, {x_start, y_start}, {x_end, y_end}, {x_ghost_end, y_ghost_end});

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
}
TEST_CASE("ConstView squarebracket operator const", "[ConstView]")
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

	const ConstView<2> v(data, {x_stride, y_stride}, {x_ghost_start, y_ghost_start}, {x_start, y_start}, {x_end, y_end}, {x_ghost_end, y_ghost_end});

	for (int yi = y_ghost_start - 1; yi <= y_ghost_end + 1; yi++) {
		for (int xi = x_ghost_start - 1; xi <= x_ghost_end + 1; xi++) {
			if (xi < x_ghost_start || xi > x_ghost_end || yi < y_ghost_start || yi > y_ghost_end) {
				//oob coord
				if constexpr (ENABLE_DEBUG) {
					CHECK_THROWS_AS((v[{xi, yi}]), RuntimeError);
				}
			} else {
				CHECK(&v[{xi, yi}] == data + xi * x_stride + yi * y_stride);
			}
		}
	}
}
TEST_CASE("ConstView parens operator const", "[ConstView]")
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

	const ConstView<2> v(data, {x_stride, y_stride}, {x_ghost_start, y_ghost_start}, {x_start, y_start}, {x_end, y_end}, {x_ghost_end, y_ghost_end});

	for (int yi = y_ghost_start - 1; yi <= y_ghost_end + 1; yi++) {
		for (int xi = x_ghost_start - 1; xi <= x_ghost_end + 1; xi++) {
			if (xi < x_ghost_start || xi > x_ghost_end || yi < y_ghost_start || yi > y_ghost_end) {
				//oob coord
				if constexpr (ENABLE_DEBUG) {
					CHECK_THROWS_AS((v(xi, yi)), RuntimeError);
				}
			} else {
				CHECK(&v(xi, yi) == data + xi * x_stride + yi * y_stride);
			}
		}
	}
}
TEST_CASE("ConstView set const", "[ConstView]")
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

	const ConstView<2> v(data, {x_stride, y_stride}, {x_ghost_start, y_ghost_start}, {x_start, y_start}, {x_end, y_end}, {x_ghost_end, y_ghost_end});

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