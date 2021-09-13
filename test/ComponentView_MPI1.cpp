#include <ThunderEgg/Vector.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
using namespace std;
using namespace ThunderEgg;
TEST_CASE("ComponentView constructor", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 2> lengths = {nx, ny};
	array<int, 2> strides = {1, nx + 2 * num_ghost};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	vector<double> vec(size);
	iota(vec.begin(), vec.end(), 0);

	ComponentView<double, 2> ld(vec.data(), strides, lengths, num_ghost);

	for (int i = 0; i < 2; i++) {
		CHECK(ld.getStrides()[i] == strides[i]);
		CHECK(ld.getStart()[i] == 0);
		CHECK(ld.getEnd()[i] == lengths[i] - 1);
		CHECK(ld.getGhostStart()[i] == -num_ghost);
		CHECK(ld.getGhostEnd()[i] == lengths[i] - 1 + num_ghost);
	}
}
TEST_CASE("ComponentView<double,2> getSliceOn<1>", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 2> lengths = {nx, ny};
	array<int, 2> strides = {1, nx + 2 * num_ghost};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	ComponentView<double, 2> ld(data, strides, lengths, num_ghost);

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		View<double, 1> slice = ld.getSliceOn(Side<2>::west(), {xi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			CHECK(&ld[{xi, yi}] == &slice[{yi}]);
		}
	}

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		View<double, 1> slice = ld.getSliceOn(Side<2>::east(), {xi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			CHECK(&ld[{nx - 1 - xi, yi}] == &slice[{yi}]);
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		View<double, 1> slice = ld.getSliceOn(Side<2>::south(), {yi});
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&ld[{xi, yi}] == &slice[{xi}]);
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		View<double, 1> slice = ld.getSliceOn(Side<2>::north(), {yi});
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&ld[{xi, ny - 1 - yi}] == &slice[{xi}]);
		}
	}
}
TEST_CASE("ComponentView<double,3> getSliceOn<2>", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, nz};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		View<double, 2> slice = ld.getSliceOn(Side<3>::west(), {xi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&ld[{xi, yi, zi}] == &slice[{yi, zi}]);
			}
		}
	}

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		View<double, 2> slice = ld.getSliceOn(Side<3>::east(), {xi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&ld[{nx - 1 - xi, yi, zi}] == &slice[{yi, zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		View<double, 2> slice = ld.getSliceOn(Side<3>::south(), {yi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, yi, zi}] == &slice[{xi, zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		View<double, 2> slice = ld.getSliceOn(Side<3>::north(), {yi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, ny - 1 - yi, zi}] == &slice[{xi, zi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		View<double, 2> slice = ld.getSliceOn(Side<3>::bottom(), {zi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, yi, zi}] == &slice[{xi, yi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		View<double, 2> slice = ld.getSliceOn(Side<3>::top(), {zi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, yi, nz - 1 - zi}] == &slice[{xi, yi}]);
			}
		}
	}
}
TEST_CASE("ComponentView<double,3> getSliceOn<1>", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, nz};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			View<double, 1> slice = ld.getSliceOn(Edge::bs(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, yi, zi}] == &slice[{xi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			View<double, 1> slice = ld.getSliceOn(Edge::tn(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, ny - 1 - yi, nz - 1 - zi}] == &slice[{xi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			View<double, 1> slice = ld.getSliceOn(Edge::bn(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, ny - 1 - yi, zi}] == &slice[{xi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			View<double, 1> slice = ld.getSliceOn(Edge::ts(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, yi, nz - 1 - zi}] == &slice[{xi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 1> slice = ld.getSliceOn(Edge::bw(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&ld[{xi, yi, zi}] == &slice[{yi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 1> slice = ld.getSliceOn(Edge::te(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&ld[{nx - 1 - xi, yi, nz - 1 - zi}] == &slice[{yi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 1> slice = ld.getSliceOn(Edge::be(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&ld[{nx - 1 - xi, yi, zi}] == &slice[{yi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 1> slice = ld.getSliceOn(Edge::tw(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&ld[{xi, yi, nz - 1 - zi}] == &slice[{yi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 1> slice = ld.getSliceOn(Edge::sw(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&ld[{xi, yi, zi}] == &slice[{zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 1> slice = ld.getSliceOn(Edge::ne(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&ld[{nx - 1 - xi, ny - 1 - yi, zi}] == &slice[{zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 1> slice = ld.getSliceOn(Edge::se(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&ld[{nx - 1 - xi, yi, zi}] == &slice[{zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 1> slice = ld.getSliceOn(Edge::nw(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&ld[{xi, ny - 1 - yi, zi}] == &slice[{zi}]);
			}
		}
	}
}
TEST_CASE("ComponentView<double,3> getSliceOn<0>", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, nz};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, yi, zi}] == &(ld.getSliceOn(Corner<3>::bsw(), {xi, yi, zi}))[{}]);
				CHECK(&ld[{nx - 1 - xi, yi, zi}] == &(ld.getSliceOn(Corner<3>::bse(), {xi, yi, zi})[{}]));
				CHECK(&ld[{xi, ny - 1 - yi, zi}] == &(ld.getSliceOn(Corner<3>::bnw(), {xi, yi, zi})[{}]));
				CHECK(&ld[{nx - 1 - xi, ny - 1 - yi, zi}] == &(ld.getSliceOn(Corner<3>::bne(), {xi, yi, zi})[{}]));
				CHECK(&ld[{xi, yi, nz - 1 - zi}] == &(ld.getSliceOn(Corner<3>::tsw(), {xi, yi, zi})[{}]));
				CHECK(&ld[{nx - 1 - xi, yi, nz - 1 - zi}] == &(ld.getSliceOn(Corner<3>::tse(), {xi, yi, zi})[{}]));
				CHECK(&ld[{xi, ny - 1 - yi, nz - 1 - zi}] == &(ld.getSliceOn(Corner<3>::tnw(), {xi, yi, zi})[{}]));
				CHECK(&ld[{nx - 1 - xi, ny - 1 - yi, nz - 1 - zi}] == &(ld.getSliceOn(Corner<3>::tne(), {xi, yi, zi})[{}]));
			}
		}
	}
}
TEST_CASE("ComponentView<double,3> getSliceOn<0> const", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, nz};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	const ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, yi, zi}] == &(ld.getSliceOn(Corner<3>::bsw(), {xi, yi, zi}))[{}]);
				CHECK(&ld[{nx - 1 - xi, yi, zi}] == &(ld.getSliceOn(Corner<3>::bse(), {xi, yi, zi})[{}]));
				CHECK(&ld[{xi, ny - 1 - yi, zi}] == &(ld.getSliceOn(Corner<3>::bnw(), {xi, yi, zi})[{}]));
				CHECK(&ld[{nx - 1 - xi, ny - 1 - yi, zi}] == &(ld.getSliceOn(Corner<3>::bne(), {xi, yi, zi})[{}]));
				CHECK(&ld[{xi, yi, nz - 1 - zi}] == &(ld.getSliceOn(Corner<3>::tsw(), {xi, yi, zi})[{}]));
				CHECK(&ld[{nx - 1 - xi, yi, nz - 1 - zi}] == &(ld.getSliceOn(Corner<3>::tse(), {xi, yi, zi})[{}]));
				CHECK(&ld[{xi, ny - 1 - yi, nz - 1 - zi}] == &(ld.getSliceOn(Corner<3>::tnw(), {xi, yi, zi})[{}]));
				CHECK(&ld[{nx - 1 - xi, ny - 1 - yi, nz - 1 - zi}] == &(ld.getSliceOn(Corner<3>::tne(), {xi, yi, zi})[{}]));
			}
		}
	}
}
TEST_CASE("ComponentView<double,2> getSliceOn<0>", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 2> lengths = {nx, ny};
	array<int, 2> strides = {1, nx + 2 * num_ghost};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	ComponentView<double, 2> ld(data, strides, lengths, num_ghost);

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&ld[{xi, yi}] == &(ld.getSliceOn(Corner<2>::sw(), {xi, yi})[{}]));
			CHECK(&ld[{nx - 1 - xi, yi}] == &(ld.getSliceOn(Corner<2>::se(), {xi, yi})[{}]));
			CHECK(&ld[{xi, ny - 1 - yi}] == &(ld.getSliceOn(Corner<2>::nw(), {xi, yi})[{}]));
			CHECK(&ld[{nx - 1 - xi, ny - 1 - yi}] == &(ld.getSliceOn(Corner<2>::ne(), {xi, yi})[{}]));
		}
	}
}
TEST_CASE("ComponentView<double,2> getSliceOn<0> const", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 2> lengths = {nx, ny};
	array<int, 2> strides = {1, nx + 2 * num_ghost};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	const ComponentView<double, 2> ld(data, strides, lengths, num_ghost);

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&ld[{xi, yi}] == &(ld.getSliceOn(Corner<2>::sw(), {xi, yi})[{}]));
			CHECK(&ld[{nx - 1 - xi, yi}] == &(ld.getSliceOn(Corner<2>::se(), {xi, yi})[{}]));
			CHECK(&ld[{xi, ny - 1 - yi}] == &(ld.getSliceOn(Corner<2>::nw(), {xi, yi})[{}]));
			CHECK(&ld[{nx - 1 - xi, ny - 1 - yi}] == &(ld.getSliceOn(Corner<2>::ne(), {xi, yi})[{}]));
		}
	}
}

TEST_CASE("ComponentView<double,2> getGhostSliceOn<1>", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 2> lengths = {nx, ny};
	array<int, 2> strides = {1, nx + 2 * num_ghost};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	const ComponentView<double, 2> ld(data, strides, lengths, num_ghost);

	for (unsigned char xi = 0; xi < num_ghost; xi++) {
		INFO("xi: " << xi);
		View<double, 1> slice_w = ld.getGhostSliceOn(Side<2>::west(), {xi});
		CHECK(slice_w.getGhostStart() == std::array<int, 1>({-num_ghost}));
		CHECK(slice_w.getStart() == std::array<int, 1>({0}));
		CHECK(slice_w.getEnd() == std::array<int, 1>({ny - 1}));
		CHECK(slice_w.getGhostEnd() == std::array<int, 1>({ny - 1 + num_ghost}));

		View<double, 1> slice_e = ld.getGhostSliceOn(Side<2>::east(), {xi});
		CHECK(slice_e.getGhostStart() == std::array<int, 1>({-num_ghost}));
		CHECK(slice_e.getStart() == std::array<int, 1>({0}));
		CHECK(slice_e.getEnd() == std::array<int, 1>({ny - 1}));
		CHECK(slice_e.getGhostEnd() == std::array<int, 1>({ny - 1 + num_ghost}));

		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			CHECK(&ld[{-1 - xi, yi}] == &slice_w[{yi}]);
			CHECK(&ld[{nx + xi, yi}] == &slice_e[{yi}]);
		}
	}

	for (unsigned char yi = 0; yi < num_ghost; yi++) {
		INFO("yi: " << yi);
		View<double, 1> slice_s = ld.getGhostSliceOn(Side<2>::south(), {yi});
		CHECK(slice_s.getGhostStart() == std::array<int, 1>({-num_ghost}));
		CHECK(slice_s.getStart() == std::array<int, 1>({0}));
		CHECK(slice_s.getEnd() == std::array<int, 1>({nx - 1}));
		CHECK(slice_s.getGhostEnd() == std::array<int, 1>({nx - 1 + num_ghost}));

		View<double, 1> slice_n = ld.getGhostSliceOn(Side<2>::north(), {yi});
		CHECK(slice_s.getGhostStart() == std::array<int, 1>({-num_ghost}));
		CHECK(slice_s.getStart() == std::array<int, 1>({0}));
		CHECK(slice_s.getEnd() == std::array<int, 1>({nx - 1}));
		CHECK(slice_s.getGhostEnd() == std::array<int, 1>({nx - 1 + num_ghost}));

		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&ld[{xi, -1 - yi}] == &slice_s[{xi}]);
			CHECK(&ld[{xi, ny + yi}] == &slice_n[{xi}]);
		}
	}
}
TEST_CASE("ComponentView<double,3> getGhostSliceOn<2>", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, nz};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	const ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

	for (unsigned char xi = 0; xi < num_ghost; xi++) {
		INFO("xi: " << xi);
		View<double, 2> slice_w = ld.getGhostSliceOn(Side<3>::west(), {xi});
		CHECK(slice_w.getGhostStart() == std::array<int, 2>({-num_ghost, -num_ghost}));
		CHECK(slice_w.getStart() == std::array<int, 2>({0, 0}));
		CHECK(slice_w.getEnd() == std::array<int, 2>({ny - 1, nz - 1}));
		CHECK(slice_w.getGhostEnd() == std::array<int, 2>({ny - 1 + num_ghost, nz - 1 + num_ghost}));

		View<double, 2> slice_e = ld.getGhostSliceOn(Side<3>::east(), {xi});
		CHECK(slice_e.getGhostStart() == std::array<int, 2>({-num_ghost, -num_ghost}));
		CHECK(slice_e.getStart() == std::array<int, 2>({0, 0}));
		CHECK(slice_e.getEnd() == std::array<int, 2>({ny - 1, nz - 1}));
		CHECK(slice_e.getGhostEnd() == std::array<int, 2>({ny - 1 + num_ghost, nz - 1 + num_ghost}));

		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&ld[{-1 - xi, yi, zi}] == &slice_w[{yi, zi}]);
				CHECK(&ld[{nx + xi, yi, zi}] == &slice_e[{yi, zi}]);
			}
		}
	}

	for (unsigned char yi = 0; yi < num_ghost; yi++) {
		INFO("yi: " << yi);
		View<double, 2> slice_s = ld.getGhostSliceOn(Side<3>::south(), {yi});
		CHECK(slice_s.getGhostStart() == std::array<int, 2>({-num_ghost, -num_ghost}));
		CHECK(slice_s.getStart() == std::array<int, 2>({0, 0}));
		CHECK(slice_s.getEnd() == std::array<int, 2>({nx - 1, nz - 1}));
		CHECK(slice_s.getGhostEnd() == std::array<int, 2>({nx - 1 + num_ghost, nz - 1 + num_ghost}));

		View<double, 2> slice_n = ld.getGhostSliceOn(Side<3>::north(), {yi});
		CHECK(slice_n.getGhostStart() == std::array<int, 2>({-num_ghost, -num_ghost}));
		CHECK(slice_n.getStart() == std::array<int, 2>({0, 0}));
		CHECK(slice_n.getEnd() == std::array<int, 2>({nx - 1, nz - 1}));
		CHECK(slice_n.getGhostEnd() == std::array<int, 2>({nx - 1 + num_ghost, nz - 1 + num_ghost}));

		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, -1 - yi, zi}] == &slice_s[{xi, zi}]);
				CHECK(&ld[{xi, ny + yi, zi}] == &slice_n[{xi, zi}]);
			}
		}
	}

	for (unsigned char zi = 0; zi < num_ghost; zi++) {
		INFO("zi: " << zi);
		View<double, 2> slice_b = ld.getGhostSliceOn(Side<3>::bottom(), {zi});
		CHECK(slice_b.getGhostStart() == std::array<int, 2>({-num_ghost, -num_ghost}));
		CHECK(slice_b.getStart() == std::array<int, 2>({0, 0}));
		CHECK(slice_b.getEnd() == std::array<int, 2>({nx - 1, ny - 1}));
		CHECK(slice_b.getGhostEnd() == std::array<int, 2>({nx - 1 + num_ghost, ny - 1 + num_ghost}));

		View<double, 2> slice_t = ld.getGhostSliceOn(Side<3>::top(), {zi});
		CHECK(slice_t.getGhostStart() == std::array<int, 2>({-num_ghost, -num_ghost}));
		CHECK(slice_t.getStart() == std::array<int, 2>({0, 0}));
		CHECK(slice_t.getEnd() == std::array<int, 2>({nx - 1, ny - 1}));
		CHECK(slice_t.getGhostEnd() == std::array<int, 2>({nx - 1 + num_ghost, ny - 1 + num_ghost}));

		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, yi, -1 - zi}] == &slice_b[{xi, yi}]);
				CHECK(&ld[{xi, yi, nz + zi}] == &slice_t[{xi, yi}]);
			}
		}
	}
}
TEST_CASE("ComponentView<double,3> getGhostSliceOn<1>", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, nz};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	const ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

	for (unsigned char zi = 0; zi < num_ghost; zi++) {
		INFO("zi: " << zi);
		for (unsigned char yi = 0; yi < num_ghost; yi++) {
			INFO("yi: " << yi);
			View<double, 1> slice_bs = ld.getGhostSliceOn(Edge::bs(), {yi, zi});
			CHECK(slice_bs.getGhostStart() == std::array<int, 1>({-num_ghost}));
			CHECK(slice_bs.getStart() == std::array<int, 1>({0}));
			CHECK(slice_bs.getEnd() == std::array<int, 1>({nx - 1}));
			CHECK(slice_bs.getGhostEnd() == std::array<int, 1>({nx - 1 + num_ghost}));

			View<double, 1> slice_tn = ld.getGhostSliceOn(Edge::tn(), {yi, zi});
			CHECK(slice_tn.getGhostStart() == std::array<int, 1>({-num_ghost}));
			CHECK(slice_tn.getStart() == std::array<int, 1>({0}));
			CHECK(slice_tn.getEnd() == std::array<int, 1>({nx - 1}));
			CHECK(slice_tn.getGhostEnd() == std::array<int, 1>({nx - 1 + num_ghost}));

			View<double, 1> slice_bn = ld.getGhostSliceOn(Edge::bn(), {yi, zi});
			CHECK(slice_bn.getGhostStart() == std::array<int, 1>({-num_ghost}));
			CHECK(slice_bn.getStart() == std::array<int, 1>({0}));
			CHECK(slice_bn.getEnd() == std::array<int, 1>({nx - 1}));
			CHECK(slice_bn.getGhostEnd() == std::array<int, 1>({nx - 1 + num_ghost}));

			View<double, 1> slice_ts = ld.getGhostSliceOn(Edge::ts(), {yi, zi});
			CHECK(slice_ts.getGhostStart() == std::array<int, 1>({-num_ghost}));
			CHECK(slice_ts.getStart() == std::array<int, 1>({0}));
			CHECK(slice_ts.getEnd() == std::array<int, 1>({nx - 1}));
			CHECK(slice_ts.getGhostEnd() == std::array<int, 1>({nx - 1 + num_ghost}));

			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{xi, -1 - yi, -1 - zi}] == &slice_bs[{xi}]);
				CHECK(&ld[{xi, ny + yi, nz + zi}] == &slice_tn[{xi}]);
				CHECK(&ld[{xi, ny + yi, -1 - zi}] == &slice_bn[{xi}]);
				CHECK(&ld[{xi, -1 - yi, nz + zi}] == &slice_ts[{xi}]);
			}
		}
	}

	for (unsigned char zi = 0; zi < num_ghost; zi++) {
		INFO("zi: " << zi);
		for (unsigned char xi = 0; xi < num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 1> slice_bw = ld.getGhostSliceOn(Edge::bw(), {xi, zi});
			CHECK(slice_bw.getGhostStart() == std::array<int, 1>({-num_ghost}));
			CHECK(slice_bw.getStart() == std::array<int, 1>({0}));
			CHECK(slice_bw.getEnd() == std::array<int, 1>({ny - 1}));
			CHECK(slice_bw.getGhostEnd() == std::array<int, 1>({ny - 1 + num_ghost}));

			View<double, 1> slice_te = ld.getGhostSliceOn(Edge::te(), {xi, zi});
			CHECK(slice_te.getGhostStart() == std::array<int, 1>({-num_ghost}));
			CHECK(slice_te.getStart() == std::array<int, 1>({0}));
			CHECK(slice_te.getEnd() == std::array<int, 1>({ny - 1}));
			CHECK(slice_te.getGhostEnd() == std::array<int, 1>({ny - 1 + num_ghost}));

			View<double, 1> slice_be = ld.getGhostSliceOn(Edge::be(), {xi, zi});
			CHECK(slice_be.getGhostStart() == std::array<int, 1>({-num_ghost}));
			CHECK(slice_be.getStart() == std::array<int, 1>({0}));
			CHECK(slice_be.getEnd() == std::array<int, 1>({ny - 1}));
			CHECK(slice_be.getGhostEnd() == std::array<int, 1>({ny - 1 + num_ghost}));

			View<double, 1> slice_tw = ld.getGhostSliceOn(Edge::tw(), {xi, zi});
			CHECK(slice_tw.getGhostStart() == std::array<int, 1>({-num_ghost}));
			CHECK(slice_tw.getStart() == std::array<int, 1>({0}));
			CHECK(slice_tw.getEnd() == std::array<int, 1>({ny - 1}));
			CHECK(slice_tw.getGhostEnd() == std::array<int, 1>({ny - 1 + num_ghost}));

			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&ld[{-1 - xi, yi, -1 - zi}] == &slice_bw[{yi}]);
				CHECK(&ld[{nx + xi, yi, nz + zi}] == &slice_te[{yi}]);
				CHECK(&ld[{nx + xi, yi, -1 - zi}] == &slice_be[{yi}]);
				CHECK(&ld[{-1 - xi, yi, nz + zi}] == &slice_tw[{yi}]);
			}
		}
	}

	for (unsigned char yi = 0; yi < num_ghost; yi++) {
		INFO("yi: " << yi);
		for (unsigned char xi = 0; xi < num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 1> slice_sw = ld.getGhostSliceOn(Edge::sw(), {xi, yi});
			CHECK(slice_sw.getGhostStart() == std::array<int, 1>({-num_ghost}));
			CHECK(slice_sw.getStart() == std::array<int, 1>({0}));
			CHECK(slice_sw.getEnd() == std::array<int, 1>({nz - 1}));
			CHECK(slice_sw.getGhostEnd() == std::array<int, 1>({nz - 1 + num_ghost}));

			View<double, 1> slice_ne = ld.getGhostSliceOn(Edge::ne(), {xi, yi});
			CHECK(slice_ne.getGhostStart() == std::array<int, 1>({-num_ghost}));
			CHECK(slice_ne.getStart() == std::array<int, 1>({0}));
			CHECK(slice_ne.getEnd() == std::array<int, 1>({nz - 1}));
			CHECK(slice_ne.getGhostEnd() == std::array<int, 1>({nz - 1 + num_ghost}));

			View<double, 1> slice_se = ld.getGhostSliceOn(Edge::se(), {xi, yi});
			CHECK(slice_se.getGhostStart() == std::array<int, 1>({-num_ghost}));
			CHECK(slice_se.getStart() == std::array<int, 1>({0}));
			CHECK(slice_se.getEnd() == std::array<int, 1>({nz - 1}));
			CHECK(slice_se.getGhostEnd() == std::array<int, 1>({nz - 1 + num_ghost}));

			View<double, 1> slice_nw = ld.getGhostSliceOn(Edge::nw(), {xi, yi});
			CHECK(slice_nw.getGhostStart() == std::array<int, 1>({-num_ghost}));
			CHECK(slice_nw.getStart() == std::array<int, 1>({0}));
			CHECK(slice_nw.getEnd() == std::array<int, 1>({nz - 1}));
			CHECK(slice_nw.getGhostEnd() == std::array<int, 1>({nz - 1 + num_ghost}));

			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&ld[{-1 - xi, -1 - yi, zi}] == &slice_sw[{zi}]);
				CHECK(&ld[{nx + xi, ny + yi, zi}] == &slice_ne[{zi}]);
				CHECK(&ld[{nx + xi, -1 - yi, zi}] == &slice_se[{zi}]);
				CHECK(&ld[{-1 - xi, ny + yi, zi}] == &slice_nw[{zi}]);
			}
		}
	}
}
TEST_CASE("ComponentView<double,3> getGhostSliceOn<0>", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, nz};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	const ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

	for (unsigned char zi = 0; zi < num_ghost; zi++) {
		INFO("zi: " << zi);
		for (unsigned char yi = 0; yi < num_ghost; yi++) {
			INFO("yi: " << yi);
			for (unsigned char xi = 0; xi < num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&ld[{-1 - xi, -1 - yi, -1 - zi}] == &(ld.getGhostSliceOn(Corner<3>::bsw(), {xi, yi, zi}))[{}]);
				CHECK(&ld[{nx + xi, -1 - yi, -1 - zi}] == &(ld.getGhostSliceOn(Corner<3>::bse(), {xi, yi, zi})[{}]));
				CHECK(&ld[{-1 - xi, ny + yi, -1 - zi}] == &(ld.getGhostSliceOn(Corner<3>::bnw(), {xi, yi, zi})[{}]));
				CHECK(&ld[{nx + xi, ny + yi, -1 - zi}] == &(ld.getGhostSliceOn(Corner<3>::bne(), {xi, yi, zi})[{}]));
				CHECK(&ld[{-1 - xi, -1 - yi, nz + zi}] == &(ld.getGhostSliceOn(Corner<3>::tsw(), {xi, yi, zi})[{}]));
				CHECK(&ld[{nx + xi, -1 - yi, nz + zi}] == &(ld.getGhostSliceOn(Corner<3>::tse(), {xi, yi, zi})[{}]));
				CHECK(&ld[{-1 - xi, ny + yi, nz + zi}] == &(ld.getGhostSliceOn(Corner<3>::tnw(), {xi, yi, zi})[{}]));
				CHECK(&ld[{nx + xi, ny + yi, nz + zi}] == &(ld.getGhostSliceOn(Corner<3>::tne(), {xi, yi, zi})[{}]));
			}
		}
	}
}
TEST_CASE("ComponentView<double,2> getGhostSliceOn<0>", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 2> lengths = {nx, ny};
	array<int, 2> strides = {1, nx + 2 * num_ghost};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	const ComponentView<double, 2> ld(data, strides, lengths, num_ghost);

	for (unsigned char yi = 0; yi < num_ghost; yi++) {
		INFO("yi: " << yi);
		for (unsigned char xi = 0; xi < num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&ld[{-xi - 1, -yi - 1}] == &(ld.getGhostSliceOn(Corner<2>::sw(), {xi, yi})[{}]));
			CHECK(&ld[{nx + xi, -yi - 1}] == &(ld.getGhostSliceOn(Corner<2>::se(), {xi, yi})[{}]));
			CHECK(&ld[{-xi - 1, ny + yi}] == &(ld.getGhostSliceOn(Corner<2>::nw(), {xi, yi})[{}]));
			CHECK(&ld[{nx + xi, ny + yi}] == &(ld.getGhostSliceOn(Corner<2>::ne(), {xi, yi})[{}]));
		}
	}
}
TEST_CASE("ComponentView squarebracket operator", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 2> lengths = {nx, ny};
	array<int, 2> strides = {1, nx + 2 * num_ghost};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	ComponentView<double, 2> v(data, strides, lengths, num_ghost);

	double *start = &v[{0, 0}];

	for (int yi = -num_ghost - 1; yi < ny + num_ghost + 1; yi++) {
		for (int xi = -num_ghost - 1; xi < nx + num_ghost + 1; xi++) {
			if (xi < -num_ghost || xi >= nx + num_ghost || yi < -num_ghost || yi >= ny + num_ghost) {
				//oob coord
				if constexpr (ENABLE_DEBUG) {
					CHECK_THROWS_AS((v[{xi, yi}]), RuntimeError);
				}
			} else {
				CHECK(&v[{xi, yi}] == start + xi + yi * (nx + 2 * num_ghost));
			}
		}
	}
}
TEST_CASE("ComponentView squarebracket operator const", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 2> lengths = {nx, ny};
	array<int, 2> strides = {1, nx + 2 * num_ghost};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	const ComponentView<double, 2> v(data, strides, lengths, num_ghost);

	const double *start = &v[{0, 0}];
	for (int yi = -num_ghost - 1; yi < ny + num_ghost + 1; yi++) {
		for (int xi = -num_ghost - 1; xi < nx + num_ghost + 1; xi++) {
			if (xi < -num_ghost || xi >= nx + num_ghost || yi < -num_ghost || yi >= ny + num_ghost) {
				//oob coord
				if constexpr (ENABLE_DEBUG) {
					CHECK_THROWS_AS((v[{xi, yi}]), RuntimeError);
				}
			} else {
				CHECK(&v[{xi, yi}] == start + xi + yi * (nx + 2 * num_ghost));
			}
		}
	}
}
TEST_CASE("ComponentView set", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 2> lengths = {nx, ny};
	array<int, 2> strides = {1, nx + 2 * num_ghost};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	ComponentView<double, 2> v(data, strides, lengths, num_ghost);

	double value = 0;
	for (int yi = -num_ghost - 1; yi < ny + num_ghost + 1; yi++) {
		for (int xi = -num_ghost - 1; xi < nx + num_ghost + 1; xi++) {
			if ((xi < -num_ghost) || (xi >= nx + num_ghost) || (yi < -num_ghost) || (yi >= ny + num_ghost)) {
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
TEST_CASE("ComponentView set const", "[ComponentView]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 2> lengths = {nx, ny};
	array<int, 2> strides = {1, nx + 2 * num_ghost};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	ComponentView<const double, 2> v(data, strides, lengths, num_ghost);

	double value = 0;
	for (int yi = -num_ghost - 1; yi < ny + num_ghost + 1; yi++) {
		for (int xi = -num_ghost - 1; xi < nx + num_ghost + 1; xi++) {
			if (xi >= 0 && xi < nx && yi >= 0 && yi < ny) {
				//intertior coord
				if constexpr (ENABLE_DEBUG) {
					CHECK_THROWS_AS(v.set({xi, yi}, value), RuntimeError);
				}
			} else if (xi < -num_ghost || xi >= nx + num_ghost || yi < -num_ghost || yi >= ny + num_ghost) {
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
TEST_CASE("ComponentView implicit conversion to const type", "[View]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	array<int, 2> lengths = {nx, ny};
	array<int, 2> strides = {1, nx + 2 * num_ghost};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	double data[size];

	ComponentView<double, 2>       v(data, strides, lengths, num_ghost);
	ComponentView<const double, 2> vc = v;

	CHECK(vc.getGhostStart() == v.getGhostStart());
	CHECK(vc.getStart() == v.getStart());
	CHECK(vc.getEnd() == v.getEnd());
	CHECK(vc.getGhostEnd() == v.getGhostEnd());
	CHECK(vc.getStrides() == v.getStrides());
	CHECK(&vc[vc.getGhostStart()] == &v[v.getGhostStart()]);
}