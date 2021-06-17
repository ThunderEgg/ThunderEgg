#include <ThunderEgg/PatchView.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
using namespace std;
using namespace ThunderEgg;
TEST_CASE("PatchView constructor", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, num_components};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (ny + 2 * num_ghost) * (nx + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	vector<double> vec(size);
	iota(vec.begin(), vec.end(), 0);

	PatchView<double, 2> pview(vec.data(), strides, lengths, num_ghost);

	for (int i = 0; i < 2; i++) {
		CHECK(pview.getStrides()[i] == strides[i]);
		CHECK(pview.getStart()[i] == 0);
		CHECK(pview.getEnd()[i] == lengths[i] - 1);
		CHECK(pview.getGhostStart()[i] == -num_ghost);
		CHECK(pview.getGhostEnd()[i] == lengths[i] - 1 + num_ghost);
	}
	CHECK(pview.getStrides()[2] == strides[2]);
	CHECK(pview.getStart()[2] == 0);
	CHECK(pview.getEnd()[2] == lengths[2] - 1);
	CHECK(pview.getGhostStart()[2] == 0);
	CHECK(pview.getGhostEnd()[2] == lengths[2] - 1);
}
TEST_CASE("PatchView<double,2> getSliceOn<1>", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, num_components};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (ny + 2 * num_ghost) * (nx + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 2> pview(data, strides, lengths, num_ghost);

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		View<double, 2> slice = pview.getSliceOn(Side<2>::west(), {xi});
		CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
		CHECK(slice.getStart() == array<int, 2>({0, 0}));
		CHECK(slice.getEnd() == array<int, 2>({ny - 1, num_components - 1}));
		CHECK(slice.getGhostEnd() == array<int, 2>({ny - 1 + num_ghost, num_components - 1}));
		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pview[{xi, yi, c}] == &slice[{yi, c}]);
			}
		}
	}

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		View<double, 2> slice = pview.getSliceOn(Side<2>::east(), {xi});
		CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
		CHECK(slice.getStart() == array<int, 2>({0, 0}));
		CHECK(slice.getEnd() == array<int, 2>({ny - 1, num_components - 1}));
		CHECK(slice.getGhostEnd() == array<int, 2>({ny - 1 + num_ghost, num_components - 1}));
		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pview[{nx - 1 - xi, yi, c}] == &slice[{yi, c}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		View<double, 2> slice = pview.getSliceOn(Side<2>::south(), {yi});
		CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
		CHECK(slice.getStart() == array<int, 2>({0, 0}));
		CHECK(slice.getEnd() == array<int, 2>({nx - 1, num_components - 1}));
		CHECK(slice.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, num_components - 1}));
		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pview[{xi, yi, c}] == &slice[{xi, c}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		View<double, 2> slice = pview.getSliceOn(Side<2>::north(), {yi});
		CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
		CHECK(slice.getStart() == array<int, 2>({0, 0}));
		CHECK(slice.getEnd() == array<int, 2>({nx - 1, num_components - 1}));
		CHECK(slice.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, num_components - 1}));
		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pview[{xi, ny - 1 - yi, c}] == &slice[{xi, c}]);
			}
		}
	}
}
TEST_CASE("PatchView<double,3> getSliceOn<2>", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto nz             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 4> lengths = {nx, ny, nz, num_components};
	array<int, 4> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost), (nx + 2 * num_ghost) * (ny + 2 * num_ghost) * (nz + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 3> pview(data, strides, lengths, num_ghost);

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		View<double, 3> slice = pview.getSliceOn(Side<3>::west(), {xi});
		CHECK(slice.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, 0}));
		CHECK(slice.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice.getEnd() == array<int, 3>({ny - 1, nz - 1, num_components - 1}));
		CHECK(slice.getGhostEnd() == array<int, 3>({ny - 1 + num_ghost, nz - 1 + num_ghost, num_components - 1}));

		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
					INFO("yi: " << yi);
					CHECK(&pview[{xi, yi, zi, c}] == &slice[{yi, zi, c}]);
				}
			}
		}
	}

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		View<double, 3> slice = pview.getSliceOn(Side<3>::east(), {xi});
		CHECK(slice.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, 0}));
		CHECK(slice.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice.getEnd() == array<int, 3>({ny - 1, nz - 1, num_components - 1}));
		CHECK(slice.getGhostEnd() == array<int, 3>({ny - 1 + num_ghost, nz - 1 + num_ghost, num_components - 1}));

		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
					INFO("yi: " << yi);
					CHECK(&pview[{nx - 1 - xi, yi, zi, c}] == &slice[{yi, zi, c}]);
				}
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		View<double, 3> slice = pview.getSliceOn(Side<3>::south(), {yi});
		CHECK(slice.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, 0}));
		CHECK(slice.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice.getEnd() == array<int, 3>({nx - 1, nz - 1, num_components - 1}));
		CHECK(slice.getGhostEnd() == array<int, 3>({nx - 1 + num_ghost, nz - 1 + num_ghost, num_components - 1}));

		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
					INFO("xi: " << xi);
					CHECK(&pview[{xi, yi, zi, c}] == &slice[{xi, zi, c}]);
				}
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		View<double, 3> slice = pview.getSliceOn(Side<3>::north(), {yi});
		CHECK(slice.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, 0}));
		CHECK(slice.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice.getEnd() == array<int, 3>({nx - 1, nz - 1, num_components - 1}));
		CHECK(slice.getGhostEnd() == array<int, 3>({nx - 1 + num_ghost, nz - 1 + num_ghost, num_components - 1}));

		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
					INFO("xi: " << xi);
					CHECK(&pview[{xi, ny - 1 - yi, zi, c}] == &slice[{xi, zi, c}]);
				}
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		View<double, 3> slice = pview.getSliceOn(Side<3>::bottom(), {zi});
		CHECK(slice.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, 0}));
		CHECK(slice.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice.getEnd() == array<int, 3>({nx - 1, ny - 1, num_components - 1}));
		CHECK(slice.getGhostEnd() == array<int, 3>({nx - 1 + num_ghost, ny - 1 + num_ghost, num_components - 1}));

		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
					INFO("xi: " << xi);
					CHECK(&pview[{xi, yi, zi, c}] == &slice[{xi, yi, c}]);
				}
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		View<double, 3> slice = pview.getSliceOn(Side<3>::top(), {zi});
		CHECK(slice.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, 0}));
		CHECK(slice.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice.getEnd() == array<int, 3>({nx - 1, ny - 1, num_components - 1}));
		CHECK(slice.getGhostEnd() == array<int, 3>({nx - 1 + num_ghost, ny - 1 + num_ghost, num_components - 1}));

		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
					INFO("xi: " << xi);
					CHECK(&pview[{xi, yi, nz - 1 - zi, c}] == &slice[{xi, yi, c}]);
				}
			}
		}
	}
}
TEST_CASE("PatchView<double,3> getSliceOn<1>", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto nz             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 4> lengths = {nx, ny, nz, num_components};
	array<int, 4> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost), (nx + 2 * num_ghost) * (ny + 2 * num_ghost) * (nz + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 3> pview(data, strides, lengths, num_ghost);

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			View<double, 2> slice = pview.getSliceOn(Edge::bs(), {yi, zi});
			CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice.getStart() == array<int, 2>({0, 0}));
			CHECK(slice.getEnd() == array<int, 2>({nx - 1, num_components - 1}));
			CHECK(slice.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
					INFO("xi: " << xi);
					CHECK(&pview[{xi, yi, zi, c}] == &slice[{xi, c}]);
				}
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			View<double, 2> slice = pview.getSliceOn(Edge::tn(), {yi, zi});
			CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice.getStart() == array<int, 2>({0, 0}));
			CHECK(slice.getEnd() == array<int, 2>({nx - 1, num_components - 1}));
			CHECK(slice.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
					INFO("xi: " << xi);
					CHECK(&pview[{xi, ny - 1 - yi, nz - 1 - zi, c}] == &slice[{xi, c}]);
				}
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			View<double, 2> slice = pview.getSliceOn(Edge::bn(), {yi, zi});
			CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice.getStart() == array<int, 2>({0, 0}));
			CHECK(slice.getEnd() == array<int, 2>({nx - 1, num_components - 1}));
			CHECK(slice.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
					INFO("xi: " << xi);
					CHECK(&pview[{xi, ny - 1 - yi, zi, c}] == &slice[{xi, c}]);
				}
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			View<double, 2> slice = pview.getSliceOn(Edge::ts(), {yi, zi});
			CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice.getStart() == array<int, 2>({0, 0}));
			CHECK(slice.getEnd() == array<int, 2>({nx - 1, num_components - 1}));
			CHECK(slice.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
					INFO("xi: " << xi);
					CHECK(&pview[{xi, yi, nz - 1 - zi, c}] == &slice[{xi, c}]);
				}
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 2> slice = pview.getSliceOn(Edge::bw(), {xi, zi});
			CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice.getStart() == array<int, 2>({0, 0}));
			CHECK(slice.getEnd() == array<int, 2>({ny - 1, num_components - 1}));
			CHECK(slice.getGhostEnd() == array<int, 2>({ny - 1 + num_ghost, num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
					INFO("yi: " << yi);
					CHECK(&pview[{xi, yi, zi, c}] == &slice[{yi, c}]);
				}
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 2> slice = pview.getSliceOn(Edge::te(), {xi, zi});
			CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice.getStart() == array<int, 2>({0, 0}));
			CHECK(slice.getEnd() == array<int, 2>({ny - 1, num_components - 1}));
			CHECK(slice.getGhostEnd() == array<int, 2>({ny - 1 + num_ghost, num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
					INFO("yi: " << yi);
					CHECK(&pview[{nx - 1 - xi, yi, nz - 1 - zi, c}] == &slice[{yi, c}]);
				}
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 2> slice = pview.getSliceOn(Edge::be(), {xi, zi});
			CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice.getStart() == array<int, 2>({0, 0}));
			CHECK(slice.getEnd() == array<int, 2>({ny - 1, num_components - 1}));
			CHECK(slice.getGhostEnd() == array<int, 2>({ny - 1 + num_ghost, num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
					INFO("yi: " << yi);
					CHECK(&pview[{nx - 1 - xi, yi, zi, c}] == &slice[{yi, c}]);
				}
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 2> slice = pview.getSliceOn(Edge::tw(), {xi, zi});
			CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice.getStart() == array<int, 2>({0, 0}));
			CHECK(slice.getEnd() == array<int, 2>({ny - 1, num_components - 1}));
			CHECK(slice.getGhostEnd() == array<int, 2>({ny - 1 + num_ghost, num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
					INFO("yi: " << yi);
					CHECK(&pview[{xi, yi, nz - 1 - zi, c}] == &slice[{yi, c}]);
				}
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 2> slice = pview.getSliceOn(Edge::sw(), {xi, yi});
			CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice.getStart() == array<int, 2>({0, 0}));
			CHECK(slice.getEnd() == array<int, 2>({nz - 1, num_components - 1}));
			CHECK(slice.getGhostEnd() == array<int, 2>({nz - 1 + num_ghost, num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
					INFO("zi: " << zi);
					CHECK(&pview[{xi, yi, zi, c}] == &slice[{zi, c}]);
				}
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 2> slice = pview.getSliceOn(Edge::ne(), {xi, yi});
			CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice.getStart() == array<int, 2>({0, 0}));
			CHECK(slice.getEnd() == array<int, 2>({nz - 1, num_components - 1}));
			CHECK(slice.getGhostEnd() == array<int, 2>({nz - 1 + num_ghost, num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
					INFO("zi: " << zi);
					CHECK(&pview[{nx - 1 - xi, ny - 1 - yi, zi, c}] == &slice[{zi, c}]);
				}
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 2> slice = pview.getSliceOn(Edge::se(), {xi, yi});
			CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice.getStart() == array<int, 2>({0, 0}));
			CHECK(slice.getEnd() == array<int, 2>({nz - 1, num_components - 1}));
			CHECK(slice.getGhostEnd() == array<int, 2>({nz - 1 + num_ghost, num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
					INFO("zi: " << zi);
					CHECK(&pview[{nx - 1 - xi, yi, zi, c}] == &slice[{zi, c}]);
				}
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 2> slice = pview.getSliceOn(Edge::nw(), {xi, yi});
			CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice.getStart() == array<int, 2>({0, 0}));
			CHECK(slice.getEnd() == array<int, 2>({nz - 1, num_components - 1}));
			CHECK(slice.getGhostEnd() == array<int, 2>({nz - 1 + num_ghost, num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
					INFO("zi: " << zi);
					CHECK(&pview[{xi, ny - 1 - yi, zi, c}] == &slice[{zi, c}]);
				}
			}
		}
	}
}
TEST_CASE("PatchView<double,3> getSliceOn<0>", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto nz             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 4> lengths = {nx, ny, nz, num_components};
	array<int, 4> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost), (nx + 2 * num_ghost) * (ny + 2 * num_ghost) * (nz + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 3> pview(data, strides, lengths, num_ghost);

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				View<double, 1> bsw = pview.getSliceOn(Corner<3>::bsw(), {xi, yi, zi});
				CHECK(bsw.getGhostStart() == array<int, 1>({0}));
				CHECK(bsw.getStart() == array<int, 1>({0}));
				CHECK(bsw.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(bsw.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> bse = pview.getSliceOn(Corner<3>::bse(), {xi, yi, zi});
				CHECK(bse.getGhostStart() == array<int, 1>({0}));
				CHECK(bse.getStart() == array<int, 1>({0}));
				CHECK(bse.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(bse.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> bnw = pview.getSliceOn(Corner<3>::bnw(), {xi, yi, zi});
				CHECK(bnw.getGhostStart() == array<int, 1>({0}));
				CHECK(bnw.getStart() == array<int, 1>({0}));
				CHECK(bnw.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(bnw.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> bne = pview.getSliceOn(Corner<3>::bne(), {xi, yi, zi});
				CHECK(bne.getGhostStart() == array<int, 1>({0}));
				CHECK(bne.getStart() == array<int, 1>({0}));
				CHECK(bne.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(bne.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> tsw = pview.getSliceOn(Corner<3>::tsw(), {xi, yi, zi});
				CHECK(tsw.getGhostStart() == array<int, 1>({0}));
				CHECK(tsw.getStart() == array<int, 1>({0}));
				CHECK(tsw.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(tsw.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> tse = pview.getSliceOn(Corner<3>::tse(), {xi, yi, zi});
				CHECK(tse.getGhostStart() == array<int, 1>({0}));
				CHECK(tse.getStart() == array<int, 1>({0}));
				CHECK(tse.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(tse.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> tnw = pview.getSliceOn(Corner<3>::tnw(), {xi, yi, zi});
				CHECK(tnw.getGhostStart() == array<int, 1>({0}));
				CHECK(tnw.getStart() == array<int, 1>({0}));
				CHECK(tnw.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(tnw.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> tne = pview.getSliceOn(Corner<3>::tne(), {xi, yi, zi});
				CHECK(tne.getGhostStart() == array<int, 1>({0}));
				CHECK(tne.getStart() == array<int, 1>({0}));
				CHECK(tne.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(tne.getGhostEnd() == array<int, 1>({num_components - 1}));

				for (int c = 0; c < num_components; c++) {
					INFO("c: " << c);
					CHECK(&pview[{xi, yi, zi, c}] == &bsw[{c}]);
					CHECK(&pview[{nx - 1 - xi, yi, zi, c}] == &bse[{c}]);
					CHECK(&pview[{xi, ny - 1 - yi, zi, c}] == &bnw[{c}]);
					CHECK(&pview[{nx - 1 - xi, ny - 1 - yi, zi, c}] == &bne[{c}]);
					CHECK(&pview[{xi, yi, nz - 1 - zi, c}] == &tsw[{c}]);
					CHECK(&pview[{nx - 1 - xi, yi, nz - 1 - zi, c}] == &tse[{c}]);
					CHECK(&pview[{xi, ny - 1 - yi, nz - 1 - zi, c}] == &tnw[{c}]);
					CHECK(&pview[{nx - 1 - xi, ny - 1 - yi, nz - 1 - zi, c}] == &tne[{c}]);
				}
			}
		}
	}
}
TEST_CASE("PatchView<double,2> getSliceOn<0>", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, num_components};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (ny + 2 * num_ghost) * (nx + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 2> pview(data, strides, lengths, num_ghost);

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 1> sw = pview.getSliceOn(Corner<2>::sw(), {xi, yi});
			CHECK(sw.getGhostStart() == array<int, 1>({0}));
			CHECK(sw.getStart() == array<int, 1>({0}));
			CHECK(sw.getEnd() == array<int, 1>({num_components - 1}));
			CHECK(sw.getGhostEnd() == array<int, 1>({num_components - 1}));

			View<double, 1> se = pview.getSliceOn(Corner<2>::se(), {xi, yi});
			CHECK(se.getGhostStart() == array<int, 1>({0}));
			CHECK(se.getStart() == array<int, 1>({0}));
			CHECK(se.getEnd() == array<int, 1>({num_components - 1}));
			CHECK(se.getGhostEnd() == array<int, 1>({num_components - 1}));

			View<double, 1> nw = pview.getSliceOn(Corner<2>::nw(), {xi, yi});
			CHECK(nw.getGhostStart() == array<int, 1>({0}));
			CHECK(nw.getStart() == array<int, 1>({0}));
			CHECK(nw.getEnd() == array<int, 1>({num_components - 1}));
			CHECK(nw.getGhostEnd() == array<int, 1>({num_components - 1}));

			View<double, 1> ne = pview.getSliceOn(Corner<2>::ne(), {xi, yi});
			CHECK(ne.getGhostStart() == array<int, 1>({0}));
			CHECK(ne.getStart() == array<int, 1>({0}));
			CHECK(ne.getEnd() == array<int, 1>({num_components - 1}));
			CHECK(ne.getGhostEnd() == array<int, 1>({num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				CHECK(&pview[{xi, yi, c}] == &sw[{c}]);
				CHECK(&pview[{nx - 1 - xi, yi, c}] == &se[{c}]);
				CHECK(&pview[{xi, ny - 1 - yi, c}] == &nw[{c}]);
				CHECK(&pview[{nx - 1 - xi, ny - 1 - yi, c}] == &ne[{c}]);
			}
		}
	}
}

TEST_CASE("PatchView<double,2> getGhostSliceOn<1>", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, num_components};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (ny + 2 * num_ghost) * (nx + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 2> pview(data, strides, lengths, num_ghost);

	for (unsigned char xi = 0; xi < num_ghost; xi++) {
		INFO("xi: " << xi);
		View<double, 2> slice_w = pview.getGhostSliceOn(Side<2>::west(), {xi});
		CHECK(slice_w.getGhostStart() == array<int, 2>({-num_ghost, 0}));
		CHECK(slice_w.getStart() == array<int, 2>({0, 0}));
		CHECK(slice_w.getEnd() == array<int, 2>({ny - 1, num_components - 1}));
		CHECK(slice_w.getGhostEnd() == array<int, 2>({ny - 1 + num_ghost, num_components - 1}));

		View<double, 2> slice_e = pview.getGhostSliceOn(Side<2>::east(), {xi});
		CHECK(slice_e.getGhostStart() == array<int, 2>({-num_ghost, 0}));
		CHECK(slice_e.getStart() == array<int, 2>({0, 0}));
		CHECK(slice_e.getEnd() == array<int, 2>({ny - 1, num_components - 1}));
		CHECK(slice_e.getGhostEnd() == array<int, 2>({ny - 1 + num_ghost, num_components - 1}));

		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pview[{-1 - xi, yi, c}] == &slice_w[{yi, c}]);
				CHECK(&pview[{nx + xi, yi, c}] == &slice_e[{yi, c}]);
			}
		}
	}

	for (unsigned char yi = 0; yi < num_ghost; yi++) {
		INFO("yi: " << yi);
		View<double, 2> slice_s = pview.getGhostSliceOn(Side<2>::south(), {yi});
		CHECK(slice_s.getGhostStart() == array<int, 2>({-num_ghost, 0}));
		CHECK(slice_s.getStart() == array<int, 2>({0, 0}));
		CHECK(slice_s.getEnd() == array<int, 2>({nx - 1, num_components - 1}));
		CHECK(slice_s.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, num_components - 1}));

		View<double, 2> slice_n = pview.getGhostSliceOn(Side<2>::north(), {yi});
		CHECK(slice_s.getGhostStart() == array<int, 2>({-num_ghost, 0}));
		CHECK(slice_s.getStart() == array<int, 2>({0, 0}));
		CHECK(slice_s.getEnd() == array<int, 2>({nx - 1, num_components - 1}));
		CHECK(slice_s.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, num_components - 1}));

		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pview[{xi, -1 - yi, c}] == &slice_s[{xi, c}]);
				CHECK(&pview[{xi, ny + yi, c}] == &slice_n[{xi, c}]);
			}
		}
	}
}
TEST_CASE("PatchView<double,3> getGhostSliceOn<2>", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto nz             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 4> lengths = {nx, ny, nz, num_components};
	array<int, 4> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost), (nx + 2 * num_ghost) * (ny + 2 * num_ghost) * (nz + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 3> pview(data, strides, lengths, num_ghost);

	for (unsigned char xi = 0; xi < num_ghost; xi++) {
		INFO("xi: " << xi);
		View<double, 3> slice_w = pview.getGhostSliceOn(Side<3>::west(), {xi});
		CHECK(slice_w.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, 0}));
		CHECK(slice_w.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice_w.getEnd() == array<int, 3>({ny - 1, nz - 1, num_components - 1}));
		CHECK(slice_w.getGhostEnd() == array<int, 3>({ny - 1 + num_ghost, nz - 1 + num_ghost, num_components - 1}));

		View<double, 3> slice_e = pview.getGhostSliceOn(Side<3>::east(), {xi});
		CHECK(slice_e.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, 0}));
		CHECK(slice_e.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice_e.getEnd() == array<int, 3>({ny - 1, nz - 1, num_components - 1}));
		CHECK(slice_e.getGhostEnd() == array<int, 3>({ny - 1 + num_ghost, nz - 1 + num_ghost, num_components - 1}));

		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
					INFO("yi: " << yi);
					CHECK(&pview[{-1 - xi, yi, zi, c}] == &slice_w[{yi, zi, c}]);
					CHECK(&pview[{nx + xi, yi, zi, c}] == &slice_e[{yi, zi, c}]);
				}
			}
		}
	}

	for (unsigned char yi = 0; yi < num_ghost; yi++) {
		INFO("yi: " << yi);
		View<double, 3> slice_s = pview.getGhostSliceOn(Side<3>::south(), {yi});
		CHECK(slice_s.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, 0}));
		CHECK(slice_s.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice_s.getEnd() == array<int, 3>({nx - 1, nz - 1, num_components - 1}));
		CHECK(slice_s.getGhostEnd() == array<int, 3>({nx - 1 + num_ghost, nz - 1 + num_ghost, num_components - 1}));

		View<double, 3> slice_n = pview.getGhostSliceOn(Side<3>::north(), {yi});
		CHECK(slice_n.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, 0}));
		CHECK(slice_n.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice_n.getEnd() == array<int, 3>({nx - 1, nz - 1, num_components - 1}));
		CHECK(slice_n.getGhostEnd() == array<int, 3>({nx - 1 + num_ghost, nz - 1 + num_ghost, num_components - 1}));

		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
					INFO("xi: " << xi);
					CHECK(&pview[{xi, -1 - yi, zi, c}] == &slice_s[{xi, zi, c}]);
					CHECK(&pview[{xi, ny + yi, zi, c}] == &slice_n[{xi, zi, c}]);
				}
			}
		}
	}

	for (unsigned char zi = 0; zi < num_ghost; zi++) {
		INFO("zi: " << zi);
		View<double, 3> slice_b = pview.getGhostSliceOn(Side<3>::bottom(), {zi});
		CHECK(slice_b.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, 0}));
		CHECK(slice_b.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice_b.getEnd() == array<int, 3>({nx - 1, ny - 1, num_components - 1}));
		CHECK(slice_b.getGhostEnd() == array<int, 3>({nx - 1 + num_ghost, ny - 1 + num_ghost, num_components - 1}));

		View<double, 3> slice_t = pview.getGhostSliceOn(Side<3>::top(), {zi});
		CHECK(slice_t.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, 0}));
		CHECK(slice_t.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice_t.getEnd() == array<int, 3>({nx - 1, ny - 1, num_components - 1}));
		CHECK(slice_t.getGhostEnd() == array<int, 3>({nx - 1 + num_ghost, ny - 1 + num_ghost, num_components - 1}));

		for (int c = 0; c < num_components; c++) {
			INFO("c: " << c);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
					INFO("xi: " << xi);
					CHECK(&pview[{xi, yi, -1 - zi, c}] == &slice_b[{xi, yi, c}]);
					CHECK(&pview[{xi, yi, nz + zi, c}] == &slice_t[{xi, yi, c}]);
				}
			}
		}
	}
}
TEST_CASE("PatchView<double,3> getGhostSliceOn<1>", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto nz             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 4> lengths = {nx, ny, nz, num_components};
	array<int, 4> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost), (nx + 2 * num_ghost) * (ny + 2 * num_ghost) * (nz + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 3> pview(data, strides, lengths, num_ghost);

	for (unsigned char zi = 0; zi < num_ghost; zi++) {
		INFO("zi: " << zi);
		for (unsigned char yi = 0; yi < num_ghost; yi++) {
			INFO("yi: " << yi);
			View<double, 2> slice_bs = pview.getGhostSliceOn(Edge::bs(), {yi, zi});
			CHECK(slice_bs.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice_bs.getStart() == array<int, 2>({0, 0}));
			CHECK(slice_bs.getEnd() == array<int, 2>({nx - 1, num_components - 1}));
			CHECK(slice_bs.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, num_components - 1}));

			View<double, 2> slice_tn = pview.getGhostSliceOn(Edge::tn(), {yi, zi});
			CHECK(slice_tn.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice_tn.getStart() == array<int, 2>({0, 0}));
			CHECK(slice_tn.getEnd() == array<int, 2>({nx - 1, num_components - 1}));
			CHECK(slice_tn.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, num_components - 1}));

			View<double, 2> slice_bn = pview.getGhostSliceOn(Edge::bn(), {yi, zi});
			CHECK(slice_bn.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice_bn.getStart() == array<int, 2>({0, 0}));
			CHECK(slice_bn.getEnd() == array<int, 2>({nx - 1, num_components - 1}));
			CHECK(slice_bn.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, num_components - 1}));

			View<double, 2> slice_ts = pview.getGhostSliceOn(Edge::ts(), {yi, zi});
			CHECK(slice_ts.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice_ts.getStart() == array<int, 2>({0, 0}));
			CHECK(slice_ts.getEnd() == array<int, 2>({nx - 1, num_components - 1}));
			CHECK(slice_ts.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, num_components - 1}));

			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pview[{xi, -1 - yi, -1 - zi}] == &slice_bs[{xi}]);
				CHECK(&pview[{xi, ny + yi, nz + zi}] == &slice_tn[{xi}]);
				CHECK(&pview[{xi, ny + yi, -1 - zi}] == &slice_bn[{xi}]);
				CHECK(&pview[{xi, -1 - yi, nz + zi}] == &slice_ts[{xi}]);
			}
		}
	}

	for (unsigned char zi = 0; zi < num_ghost; zi++) {
		INFO("zi: " << zi);
		for (unsigned char xi = 0; xi < num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 2> slice_bw = pview.getGhostSliceOn(Edge::bw(), {xi, zi});
			CHECK(slice_bw.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice_bw.getStart() == array<int, 2>({0, 0}));
			CHECK(slice_bw.getEnd() == array<int, 2>({ny - 1, num_components - 1}));
			CHECK(slice_bw.getGhostEnd() == array<int, 2>({ny - 1 + num_ghost, num_components - 1}));

			View<double, 2> slice_te = pview.getGhostSliceOn(Edge::te(), {xi, zi});
			CHECK(slice_te.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice_te.getStart() == array<int, 2>({0, 0}));
			CHECK(slice_te.getEnd() == array<int, 2>({ny - 1, num_components - 1}));
			CHECK(slice_te.getGhostEnd() == array<int, 2>({ny - 1 + num_ghost, num_components - 1}));

			View<double, 2> slice_be = pview.getGhostSliceOn(Edge::be(), {xi, zi});
			CHECK(slice_be.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice_be.getStart() == array<int, 2>({0, 0}));
			CHECK(slice_be.getEnd() == array<int, 2>({ny - 1, num_components - 1}));
			CHECK(slice_be.getGhostEnd() == array<int, 2>({ny - 1 + num_ghost, num_components - 1}));

			View<double, 2> slice_tw = pview.getGhostSliceOn(Edge::tw(), {xi, zi});
			CHECK(slice_tw.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice_tw.getStart() == array<int, 2>({0, 0}));
			CHECK(slice_tw.getEnd() == array<int, 2>({ny - 1, num_components - 1}));
			CHECK(slice_tw.getGhostEnd() == array<int, 2>({ny - 1 + num_ghost, num_components - 1}));

			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pview[{-1 - xi, yi, -1 - zi}] == &slice_bw[{yi}]);
				CHECK(&pview[{nx + xi, yi, nz + zi}] == &slice_te[{yi}]);
				CHECK(&pview[{nx + xi, yi, -1 - zi}] == &slice_be[{yi}]);
				CHECK(&pview[{-1 - xi, yi, nz + zi}] == &slice_tw[{yi}]);
			}
		}
	}

	for (unsigned char yi = 0; yi < num_ghost; yi++) {
		INFO("yi: " << yi);
		for (unsigned char xi = 0; xi < num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 2> slice_sw = pview.getGhostSliceOn(Edge::sw(), {xi, yi});
			CHECK(slice_sw.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice_sw.getStart() == array<int, 2>({0, 0}));
			CHECK(slice_sw.getEnd() == array<int, 2>({nz - 1, num_components - 1}));
			CHECK(slice_sw.getGhostEnd() == array<int, 2>({nz - 1 + num_ghost, num_components - 1}));

			View<double, 2> slice_ne = pview.getGhostSliceOn(Edge::ne(), {xi, yi});
			CHECK(slice_ne.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice_ne.getStart() == array<int, 2>({0, 0}));
			CHECK(slice_ne.getEnd() == array<int, 2>({nz - 1, num_components - 1}));
			CHECK(slice_ne.getGhostEnd() == array<int, 2>({nz - 1 + num_ghost, num_components - 1}));

			View<double, 2> slice_se = pview.getGhostSliceOn(Edge::se(), {xi, yi});
			CHECK(slice_se.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice_se.getStart() == array<int, 2>({0, 0}));
			CHECK(slice_se.getEnd() == array<int, 2>({nz - 1, num_components - 1}));
			CHECK(slice_se.getGhostEnd() == array<int, 2>({nz - 1 + num_ghost, num_components - 1}));

			View<double, 2> slice_nw = pview.getGhostSliceOn(Edge::nw(), {xi, yi});
			CHECK(slice_nw.getGhostStart() == array<int, 2>({-num_ghost, 0}));
			CHECK(slice_nw.getStart() == array<int, 2>({0, 0}));
			CHECK(slice_nw.getEnd() == array<int, 2>({nz - 1, num_components - 1}));
			CHECK(slice_nw.getGhostEnd() == array<int, 2>({nz - 1 + num_ghost, num_components - 1}));

			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&pview[{-1 - xi, -1 - yi, zi}] == &slice_sw[{zi}]);
				CHECK(&pview[{nx + xi, ny + yi, zi}] == &slice_ne[{zi}]);
				CHECK(&pview[{nx + xi, -1 - yi, zi}] == &slice_se[{zi}]);
				CHECK(&pview[{-1 - xi, ny + yi, zi}] == &slice_nw[{zi}]);
			}
		}
	}
}
TEST_CASE("PatchView<double,3> getGhostSliceOn<0>", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto nz             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 4> lengths = {nx, ny, nz, num_components};
	array<int, 4> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost), (nx + 2 * num_ghost) * (ny + 2 * num_ghost) * (nz + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 3> pview(data, strides, lengths, num_ghost);

	for (unsigned char zi = 0; zi < num_ghost; zi++) {
		INFO("zi: " << zi);
		for (unsigned char yi = 0; yi < num_ghost; yi++) {
			INFO("yi: " << yi);
			for (unsigned char xi = 0; xi < num_ghost; xi++) {
				INFO("xi: " << xi);
				View<double, 1> bsw = pview.getGhostSliceOn(Corner<3>::bsw(), {xi, yi, zi});
				CHECK(bsw.getGhostStart() == array<int, 1>({0}));
				CHECK(bsw.getStart() == array<int, 1>({0}));
				CHECK(bsw.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(bsw.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> bse = pview.getGhostSliceOn(Corner<3>::bse(), {xi, yi, zi});
				CHECK(bse.getGhostStart() == array<int, 1>({0}));
				CHECK(bse.getStart() == array<int, 1>({0}));
				CHECK(bse.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(bse.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> bnw = pview.getGhostSliceOn(Corner<3>::bnw(), {xi, yi, zi});
				CHECK(bnw.getGhostStart() == array<int, 1>({0}));
				CHECK(bnw.getStart() == array<int, 1>({0}));
				CHECK(bnw.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(bnw.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> bne = pview.getGhostSliceOn(Corner<3>::bne(), {xi, yi, zi});
				CHECK(bne.getGhostStart() == array<int, 1>({0}));
				CHECK(bne.getStart() == array<int, 1>({0}));
				CHECK(bne.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(bne.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> tsw = pview.getGhostSliceOn(Corner<3>::tsw(), {xi, yi, zi});
				CHECK(tsw.getGhostStart() == array<int, 1>({0}));
				CHECK(tsw.getStart() == array<int, 1>({0}));
				CHECK(tsw.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(tsw.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> tse = pview.getGhostSliceOn(Corner<3>::tse(), {xi, yi, zi});
				CHECK(tse.getGhostStart() == array<int, 1>({0}));
				CHECK(tse.getStart() == array<int, 1>({0}));
				CHECK(tse.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(tse.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> tnw = pview.getGhostSliceOn(Corner<3>::tnw(), {xi, yi, zi});
				CHECK(tnw.getGhostStart() == array<int, 1>({0}));
				CHECK(tnw.getStart() == array<int, 1>({0}));
				CHECK(tnw.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(tnw.getGhostEnd() == array<int, 1>({num_components - 1}));

				View<double, 1> tne = pview.getGhostSliceOn(Corner<3>::tne(), {xi, yi, zi});
				CHECK(tne.getGhostStart() == array<int, 1>({0}));
				CHECK(tne.getStart() == array<int, 1>({0}));
				CHECK(tne.getEnd() == array<int, 1>({num_components - 1}));
				CHECK(tne.getGhostEnd() == array<int, 1>({num_components - 1}));

				for (int c = 0; c < num_components; c++) {
					INFO("c: " << c);
					CHECK(&pview[{-1 - xi, -1 - yi, -1 - zi, c}] == &bsw[{c}]);
					CHECK(&pview[{nx + xi, -1 - yi, -1 - zi, c}] == &bse[{c}]);
					CHECK(&pview[{-1 - xi, ny + yi, -1 - zi, c}] == &bnw[{c}]);
					CHECK(&pview[{nx + xi, ny + yi, -1 - zi, c}] == &bne[{c}]);
					CHECK(&pview[{-1 - xi, -1 - yi, nz + zi, c}] == &tsw[{c}]);
					CHECK(&pview[{nx + xi, -1 - yi, nz + zi, c}] == &tse[{c}]);
					CHECK(&pview[{-1 - xi, ny + yi, nz + zi, c}] == &tnw[{c}]);
					CHECK(&pview[{nx + xi, ny + yi, nz + zi, c}] == &tne[{c}]);
				}
			}
		}
	}
}
TEST_CASE("PatchView<double,2> getGhostSliceOn<0>", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, num_components};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (ny + 2 * num_ghost) * (nx + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 2> pview(data, strides, lengths, num_ghost);

	for (unsigned char yi = 0; yi < num_ghost; yi++) {
		INFO("yi: " << yi);
		for (unsigned char xi = 0; xi < num_ghost; xi++) {
			INFO("xi: " << xi);
			View<double, 1> sw = pview.getGhostSliceOn(Corner<2>::sw(), {xi, yi});
			CHECK(sw.getGhostStart() == array<int, 1>({0}));
			CHECK(sw.getStart() == array<int, 1>({0}));
			CHECK(sw.getEnd() == array<int, 1>({num_components - 1}));
			CHECK(sw.getGhostEnd() == array<int, 1>({num_components - 1}));

			View<double, 1> se = pview.getGhostSliceOn(Corner<2>::se(), {xi, yi});
			CHECK(se.getGhostStart() == array<int, 1>({0}));
			CHECK(se.getStart() == array<int, 1>({0}));
			CHECK(se.getEnd() == array<int, 1>({num_components - 1}));
			CHECK(se.getGhostEnd() == array<int, 1>({num_components - 1}));

			View<double, 1> nw = pview.getGhostSliceOn(Corner<2>::nw(), {xi, yi});
			CHECK(nw.getGhostStart() == array<int, 1>({0}));
			CHECK(nw.getStart() == array<int, 1>({0}));
			CHECK(nw.getEnd() == array<int, 1>({num_components - 1}));
			CHECK(nw.getGhostEnd() == array<int, 1>({num_components - 1}));

			View<double, 1> ne = pview.getGhostSliceOn(Corner<2>::ne(), {xi, yi});
			CHECK(ne.getGhostStart() == array<int, 1>({0}));
			CHECK(ne.getStart() == array<int, 1>({0}));
			CHECK(ne.getEnd() == array<int, 1>({num_components - 1}));
			CHECK(ne.getGhostEnd() == array<int, 1>({num_components - 1}));

			for (int c = 0; c < num_components; c++) {
				INFO("c: " << c);
				CHECK(&pview[{-xi - 1, -yi - 1, c}] == &sw[{c}]);
				CHECK(&pview[{nx + xi, -yi - 1, c}] == &se[{c}]);
				CHECK(&pview[{-xi - 1, ny + yi, c}] == &nw[{c}]);
				CHECK(&pview[{nx + xi, ny + yi, c}] == &ne[{c}]);
			}
		}
	}
}
TEST_CASE("PatchView squarebracket operator", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, num_components};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (ny + 2 * num_ghost) * (nx + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 2> v(data, strides, lengths, num_ghost);

	double *start = &v[{0, 0, 0}];

	for (int c = -1; c < num_components + 1; c++) {
		for (int yi = -num_ghost - 1; yi < ny + num_ghost + 1; yi++) {
			for (int xi = -num_ghost - 1; xi < nx + num_ghost + 1; xi++) {
				if (xi < -num_ghost || xi >= nx + num_ghost || yi < -num_ghost || yi >= ny + num_ghost || c < 0 || c >= num_components) {
					//oob coord
					if constexpr (ENABLE_DEBUG) {
						CHECK_THROWS_AS((v[{xi, yi, c}]), RuntimeError);
					}
				} else {
					CHECK(&v[{xi, yi, c}] == start + xi + yi * (nx + 2 * num_ghost) + c * (ny + 2 * num_ghost) * (nx + 2 * num_ghost));
				}
			}
		}
	}
}
TEST_CASE("PatchView squarebracket operator const", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, num_components};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (ny + 2 * num_ghost) * (nx + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 2> v(data, strides, lengths, num_ghost);

	const double *start = &v[{0, 0, 0}];
	for (int c = -1; c < num_components + 1; c++) {
		for (int yi = -num_ghost - 1; yi < ny + num_ghost + 1; yi++) {
			for (int xi = -num_ghost - 1; xi < nx + num_ghost + 1; xi++) {
				if (xi < -num_ghost || xi >= nx + num_ghost || yi < -num_ghost || yi >= ny + num_ghost || c < 0 || c >= num_components) {
					//oob coord
					if constexpr (ENABLE_DEBUG) {
						CHECK_THROWS_AS((v[{xi, yi, c}]), RuntimeError);
					}
				} else {
					CHECK(&v[{xi, yi, c}] == start + xi + yi * (nx + 2 * num_ghost) + c * (ny + 2 * num_ghost) * (nx + 2 * num_ghost));
				}
			}
		}
	}
}
TEST_CASE("PatchView set", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, num_components};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (ny + 2 * num_ghost) * (nx + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 2> v(data, strides, lengths, num_ghost);

	double value = 0;
	for (int c = -1; c < num_components + 1; c++) {
		for (int yi = -num_ghost - 1; yi < ny + num_ghost + 1; yi++) {
			for (int xi = -num_ghost - 1; xi < nx + num_ghost + 1; xi++) {
				if ((xi < -num_ghost) || (xi >= nx + num_ghost) || (yi < -num_ghost) || (yi >= ny + num_ghost) || c < 0 || c >= num_components) {
					//oob coord
					if constexpr (ENABLE_DEBUG) {
						CHECK_THROWS_AS(v.set({xi, yi, c}, value), RuntimeError);
					}
				} else {
					v.set({xi, yi, c}, value);
					CHECK(v[{xi, yi, c}] == value);
				}
				value++;
			}
		}
	}
}
TEST_CASE("PatchView set const", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, num_components};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (ny + 2 * num_ghost) * (nx + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<const double, 2> v(data, strides, lengths, num_ghost);

	double value = 0;
	for (int c = -1; c < num_components + 1; c++) {
		for (int yi = -num_ghost - 1; yi < ny + num_ghost + 1; yi++) {
			for (int xi = -num_ghost - 1; xi < nx + num_ghost + 1; xi++) {
				if (xi >= 0 && xi < nx && yi >= 0 && yi < ny) {
					//intertior coord
					if constexpr (ENABLE_DEBUG) {
						CHECK_THROWS_AS(v.set({xi, yi, c}, value), RuntimeError);
					}
				} else if (xi < -num_ghost || xi >= nx + num_ghost || yi < -num_ghost || yi >= ny + num_ghost || c < 0 || c >= num_components) {
					//oob coord
					if constexpr (ENABLE_DEBUG) {
						CHECK_THROWS_AS(v.set({xi, yi, c}, value), RuntimeError);
					}
				} else {
					v.set({xi, yi, c}, value);
					CHECK(v[{xi, yi, c}] == value);
				}
				value++;
			}
		}
	}
}
TEST_CASE("PatchView implicit conversion to const type", "[View]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, num_components};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (ny + 2 * num_ghost) * (nx + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 2>       v(data, strides, lengths, num_ghost);
	PatchView<const double, 2> vc = v;

	CHECK(vc.getGhostStart() == v.getGhostStart());
	CHECK(vc.getStart() == v.getStart());
	CHECK(vc.getEnd() == v.getEnd());
	CHECK(vc.getGhostEnd() == v.getGhostEnd());
	CHECK(vc.getStrides() == v.getStrides());
	CHECK(&vc[vc.getGhostStart()] == &v[v.getGhostStart()]);
}
TEST_CASE("PatchView<double,2> getComponentView", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 3> lengths = {nx, ny, num_components};
	array<int, 3> strides = {1, nx + 2 * num_ghost, (ny + 2 * num_ghost) * (nx + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 2> pview(data, strides, lengths, num_ghost);

	for (int c = 0; c < num_components; c++) {
		INFO("c: " << c);
		ComponentView<double, 2> slice = pview.getComponentView(c);
		CHECK(slice.getGhostStart() == array<int, 2>({-num_ghost, -num_ghost}));
		CHECK(slice.getStart() == array<int, 2>({0, 0}));
		CHECK(slice.getEnd() == array<int, 2>({nx - 1, ny - 1}));
		CHECK(slice.getGhostEnd() == array<int, 2>({nx - 1 + num_ghost, ny - 1 + num_ghost}));
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pview[{xi, yi, c}] == &slice[{xi, yi}]);
			}
		}
	}
}
TEST_CASE("PatchView<double,3> getComponentView", "[PatchView]")
{
	auto nx             = GENERATE(2, 3);
	auto ny             = GENERATE(2, 3);
	auto nz             = GENERATE(2, 3);
	auto num_components = GENERATE(1, 2);
	auto num_ghost      = GENERATE(0, 1, 2);

	array<int, 4> lengths = {nx, ny, nz, num_components};
	array<int, 4> strides = {1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost), (nx + 2 * num_ghost) * (ny + 2 * num_ghost) * (nz + 2 * num_ghost)};
	int           size    = 1;
	for (size_t i = 0; i < 3; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	size *= num_components;
	double data[size];

	PatchView<double, 3> pview(data, strides, lengths, num_ghost);

	for (int c = 0; c < num_components; c++) {
		INFO("c: " << c);
		ComponentView<double, 3> slice = pview.getComponentView(c);
		CHECK(slice.getGhostStart() == array<int, 3>({-num_ghost, -num_ghost, -num_ghost}));
		CHECK(slice.getStart() == array<int, 3>({0, 0, 0}));
		CHECK(slice.getEnd() == array<int, 3>({nx - 1, ny - 1, nz - 1}));
		CHECK(slice.getGhostEnd() == array<int, 3>({nx - 1 + num_ghost, ny - 1 + num_ghost, nz - 1 + num_ghost}));

		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
					INFO("xi: " << xi);
					CHECK(&pview[{xi, yi, zi, c}] == &slice[{xi, yi, zi}]);
				}
			}
		}
	}
}