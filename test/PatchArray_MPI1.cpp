#include <ThunderEgg/PatchArray.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
using namespace std;
using namespace ThunderEgg;
TEST_CASE("ComponentArray getLengths", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	CHECK(pa.getLengths()[0] == nx);
	CHECK(pa.getLengths()[1] == ny);
}
TEST_CASE("ComponentArray getStart", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	CHECK(pa.getStart()[0] == 0);
	CHECK(pa.getStart()[1] == 0);
}
TEST_CASE("ComponentArray getEnd", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	CHECK(pa.getEnd()[0] == nx - 1);
	CHECK(pa.getEnd()[1] == ny - 1);
}
TEST_CASE("ComponentArray getGhostStart", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	CHECK(pa.getGhostStart()[0] == -num_ghost);
	CHECK(pa.getGhostStart()[1] == -num_ghost);
}
TEST_CASE("ComponentArray getGhostEnd", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	CHECK(pa.getGhostEnd()[0] == nx - 1 + num_ghost);
	CHECK(pa.getGhostEnd()[1] == ny - 1 + num_ghost);
}
TEST_CASE("ComponentArray getNumGhostCells", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	CHECK(pa.getNumGhostCells() == num_ghost);
}
TEST_CASE("ComponentArray getStrides", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	CHECK(pa.getStrides()[0] == 1);
	CHECK(pa.getStrides()[1] == nx + 2 * num_ghost);
}
TEST_CASE("ComponentArray squarebracket operator", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	double *start = &pa[{0, 0}];
	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			CHECK(&pa[{xi, yi}] == start + xi + yi * (nx + 2 * num_ghost));
		}
	}
	CHECK(pa.getStrides()[0] == 1);
	CHECK(pa.getStrides()[1] == nx + 2 * num_ghost);
}
TEST_CASE("ComponentArray squarebracket operator const", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	const ComponentArray<2> pa({nx, ny}, num_ghost);

	const double *start = &pa[{0, 0}];
	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			CHECK(&pa[{xi, yi}] == start + xi + yi * (nx + 2 * num_ghost));
		}
	}
	CHECK(pa.getStrides()[0] == 1);
	CHECK(pa.getStrides()[1] == nx + 2 * num_ghost);
}
TEST_CASE("ComponentArray default is zero", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			CHECK(pa[{xi, yi}] == 0.0);
		}
	}
	CHECK(pa.getStrides()[0] == 1);
	CHECK(pa.getStrides()[1] == nx + 2 * num_ghost);
}
TEST_CASE("ComponentArray<2> getSliceOn<1>", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		View<1> slice = pa.getSliceOn(Side<2>::west(), {xi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			CHECK(&pa[{xi, yi}] == &slice[{yi}]);
		}
	}

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		View<1> slice = pa.getSliceOn(Side<2>::east(), {xi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			CHECK(&pa[{nx - 1 - xi, yi}] == &slice[{yi}]);
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		View<1> slice = pa.getSliceOn(Side<2>::south(), {yi});
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&pa[{xi, yi}] == &slice[{xi}]);
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		View<1> slice = pa.getSliceOn(Side<2>::north(), {yi});
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&pa[{xi, ny - 1 - yi}] == &slice[{xi}]);
		}
	}
}
TEST_CASE("ComponentArray<2> getSliceOn<1> const", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	const ComponentArray<2> pa({nx, ny}, num_ghost);

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		ConstView<1> slice = pa.getSliceOn(Side<2>::west(), {xi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			CHECK(&pa[{xi, yi}] == &slice[{yi}]);
		}
	}

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		ConstView<1> slice = pa.getSliceOn(Side<2>::east(), {xi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			CHECK(&pa[{nx - 1 - xi, yi}] == &slice[{yi}]);
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		ConstView<1> slice = pa.getSliceOn(Side<2>::south(), {yi});
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&pa[{xi, yi}] == &slice[{xi}]);
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		ConstView<1> slice = pa.getSliceOn(Side<2>::north(), {yi});
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&pa[{xi, ny - 1 - yi}] == &slice[{xi}]);
		}
	}
}
TEST_CASE("ComponentArray<3> getSliceOn<2>", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<3> pa({nx, ny, nz}, num_ghost);

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		View<2> slice = pa.getSliceOn(Side<3>::west(), {xi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{xi, yi, zi}] == &slice[{yi, zi}]);
			}
		}
	}

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		View<2> slice = pa.getSliceOn(Side<3>::east(), {xi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{nx - 1 - xi, yi, zi}] == &slice[{yi, zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		View<2> slice = pa.getSliceOn(Side<3>::south(), {yi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, zi}] == &slice[{xi, zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		View<2> slice = pa.getSliceOn(Side<3>::north(), {yi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, ny - 1 - yi, zi}] == &slice[{xi, zi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		View<2> slice = pa.getSliceOn(Side<3>::bottom(), {zi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, zi}] == &slice[{xi, yi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		View<2> slice = pa.getSliceOn(Side<3>::top(), {zi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, nz - 1 - zi}] == &slice[{xi, yi}]);
			}
		}
	}
}
TEST_CASE("ComponentArray<3> getSliceOn<2> const", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	const ComponentArray<3> pa({nx, ny, nz}, num_ghost);

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		ConstView<2> slice = pa.getSliceOn(Side<3>::west(), {xi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{xi, yi, zi}] == &slice[{yi, zi}]);
			}
		}
	}

	for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
		INFO("xi: " << xi);
		ConstView<2> slice = pa.getSliceOn(Side<3>::east(), {xi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{nx - 1 - xi, yi, zi}] == &slice[{yi, zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		ConstView<2> slice = pa.getSliceOn(Side<3>::south(), {yi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, zi}] == &slice[{xi, zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		ConstView<2> slice = pa.getSliceOn(Side<3>::north(), {yi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, ny - 1 - yi, zi}] == &slice[{xi, zi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		ConstView<2> slice = pa.getSliceOn(Side<3>::bottom(), {zi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, zi}] == &slice[{xi, yi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		ConstView<2> slice = pa.getSliceOn(Side<3>::top(), {zi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, nz - 1 - zi}] == &slice[{xi, yi}]);
			}
		}
	}
}
TEST_CASE("ComponentArray<3> getSliceOn<1>", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<3> pa({nx, ny, nz}, num_ghost);

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			View<1> slice = pa.getSliceOn(Edge::bs(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, zi}] == &slice[{xi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			View<1> slice = pa.getSliceOn(Edge::tn(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, ny - 1 - yi, nz - 1 - zi}] == &slice[{xi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			View<1> slice = pa.getSliceOn(Edge::bn(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, ny - 1 - yi, zi}] == &slice[{xi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			View<1> slice = pa.getSliceOn(Edge::ts(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, nz - 1 - zi}] == &slice[{xi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<1> slice = pa.getSliceOn(Edge::bw(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{xi, yi, zi}] == &slice[{yi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<1> slice = pa.getSliceOn(Edge::te(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{nx - 1 - xi, yi, nz - 1 - zi}] == &slice[{yi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<1> slice = pa.getSliceOn(Edge::be(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{nx - 1 - xi, yi, zi}] == &slice[{yi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<1> slice = pa.getSliceOn(Edge::tw(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{xi, yi, nz - 1 - zi}] == &slice[{yi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<1> slice = pa.getSliceOn(Edge::sw(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&pa[{xi, yi, zi}] == &slice[{zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<1> slice = pa.getSliceOn(Edge::ne(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&pa[{nx - 1 - xi, ny - 1 - yi, zi}] == &slice[{zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<1> slice = pa.getSliceOn(Edge::se(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&pa[{nx - 1 - xi, yi, zi}] == &slice[{zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			View<1> slice = pa.getSliceOn(Edge::nw(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&pa[{xi, ny - 1 - yi, zi}] == &slice[{zi}]);
			}
		}
	}
}
TEST_CASE("ComponentArray<3> getSliceOn<1> const", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	const ComponentArray<3> pa({nx, ny, nz}, num_ghost);

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			ConstView<1> slice = pa.getSliceOn(Edge::bs(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, zi}] == &slice[{xi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			ConstView<1> slice = pa.getSliceOn(Edge::tn(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, ny - 1 - yi, nz - 1 - zi}] == &slice[{xi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			ConstView<1> slice = pa.getSliceOn(Edge::bn(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, ny - 1 - yi, zi}] == &slice[{xi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			ConstView<1> slice = pa.getSliceOn(Edge::ts(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, nz - 1 - zi}] == &slice[{xi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			ConstView<1> slice = pa.getSliceOn(Edge::bw(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{xi, yi, zi}] == &slice[{yi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			ConstView<1> slice = pa.getSliceOn(Edge::te(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{nx - 1 - xi, yi, nz - 1 - zi}] == &slice[{yi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			ConstView<1> slice = pa.getSliceOn(Edge::be(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{nx - 1 - xi, yi, zi}] == &slice[{yi}]);
			}
		}
	}

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			ConstView<1> slice = pa.getSliceOn(Edge::tw(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{xi, yi, nz - 1 - zi}] == &slice[{yi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			ConstView<1> slice = pa.getSliceOn(Edge::sw(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&pa[{xi, yi, zi}] == &slice[{zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			ConstView<1> slice = pa.getSliceOn(Edge::ne(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&pa[{nx - 1 - xi, ny - 1 - yi, zi}] == &slice[{zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			ConstView<1> slice = pa.getSliceOn(Edge::se(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&pa[{nx - 1 - xi, yi, zi}] == &slice[{zi}]);
			}
		}
	}

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			ConstView<1> slice = pa.getSliceOn(Edge::nw(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&pa[{xi, ny - 1 - yi, zi}] == &slice[{zi}]);
			}
		}
	}
}
TEST_CASE("ComponentArray<3> getSliceOn<0>", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<3> pa({nx, ny, nz}, num_ghost);

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, zi}] == &(pa.getSliceOn(Corner<3>::bsw(), {xi, yi, zi}))[{}]);
				CHECK(&pa[{nx - 1 - xi, yi, zi}] == &(pa.getSliceOn(Corner<3>::bse(), {xi, yi, zi})[{}]));
				CHECK(&pa[{xi, ny - 1 - yi, zi}] == &(pa.getSliceOn(Corner<3>::bnw(), {xi, yi, zi})[{}]));
				CHECK(&pa[{nx - 1 - xi, ny - 1 - yi, zi}] == &(pa.getSliceOn(Corner<3>::bne(), {xi, yi, zi})[{}]));
				CHECK(&pa[{xi, yi, nz - 1 - zi}] == &(pa.getSliceOn(Corner<3>::tsw(), {xi, yi, zi})[{}]));
				CHECK(&pa[{nx - 1 - xi, yi, nz - 1 - zi}] == &(pa.getSliceOn(Corner<3>::tse(), {xi, yi, zi})[{}]));
				CHECK(&pa[{xi, ny - 1 - yi, nz - 1 - zi}] == &(pa.getSliceOn(Corner<3>::tnw(), {xi, yi, zi})[{}]));
				CHECK(&pa[{nx - 1 - xi, ny - 1 - yi, nz - 1 - zi}] == &(pa.getSliceOn(Corner<3>::tne(), {xi, yi, zi})[{}]));
			}
		}
	}
}
TEST_CASE("ComponentArray<3> getSliceOn<0> const", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	const ComponentArray<3> pa({nx, ny, nz}, num_ghost);

	for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
		INFO("zi: " << zi);
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, zi}] == &(pa.getSliceOn(Corner<3>::bsw(), {xi, yi, zi}))[{}]);
				CHECK(&pa[{nx - 1 - xi, yi, zi}] == &(pa.getSliceOn(Corner<3>::bse(), {xi, yi, zi})[{}]));
				CHECK(&pa[{xi, ny - 1 - yi, zi}] == &(pa.getSliceOn(Corner<3>::bnw(), {xi, yi, zi})[{}]));
				CHECK(&pa[{nx - 1 - xi, ny - 1 - yi, zi}] == &(pa.getSliceOn(Corner<3>::bne(), {xi, yi, zi})[{}]));
				CHECK(&pa[{xi, yi, nz - 1 - zi}] == &(pa.getSliceOn(Corner<3>::tsw(), {xi, yi, zi})[{}]));
				CHECK(&pa[{nx - 1 - xi, yi, nz - 1 - zi}] == &(pa.getSliceOn(Corner<3>::tse(), {xi, yi, zi})[{}]));
				CHECK(&pa[{xi, ny - 1 - yi, nz - 1 - zi}] == &(pa.getSliceOn(Corner<3>::tnw(), {xi, yi, zi})[{}]));
				CHECK(&pa[{nx - 1 - xi, ny - 1 - yi, nz - 1 - zi}] == &(pa.getSliceOn(Corner<3>::tne(), {xi, yi, zi})[{}]));
			}
		}
	}
}
TEST_CASE("ComponentArray<2> getSliceOn<0>", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&pa[{xi, yi}] == &(pa.getSliceOn(Corner<2>::sw(), {xi, yi})[{}]));
			CHECK(&pa[{nx - 1 - xi, yi}] == &(pa.getSliceOn(Corner<2>::se(), {xi, yi})[{}]));
			CHECK(&pa[{xi, ny - 1 - yi}] == &(pa.getSliceOn(Corner<2>::nw(), {xi, yi})[{}]));
			CHECK(&pa[{nx - 1 - xi, ny - 1 - yi}] == &(pa.getSliceOn(Corner<2>::ne(), {xi, yi})[{}]));
		}
	}
}
TEST_CASE("ComponentArray<2> getSliceOn<0> const", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	const ComponentArray<2> pa({nx, ny}, num_ghost);

	for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
		INFO("yi: " << yi);
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&pa[{xi, yi}] == &(pa.getSliceOn(Corner<2>::sw(), {xi, yi})[{}]));
			CHECK(&pa[{nx - 1 - xi, yi}] == &(pa.getSliceOn(Corner<2>::se(), {xi, yi})[{}]));
			CHECK(&pa[{xi, ny - 1 - yi}] == &(pa.getSliceOn(Corner<2>::nw(), {xi, yi})[{}]));
			CHECK(&pa[{nx - 1 - xi, ny - 1 - yi}] == &(pa.getSliceOn(Corner<2>::ne(), {xi, yi})[{}]));
		}
	}
}

TEST_CASE("ComponentArray<2> getGhostSliceOn<1>", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	for (unsigned char xi = 0; xi < num_ghost; xi++) {
		INFO("xi: " << xi);
		View<1> slice_w = pa.getGhostSliceOn(Side<2>::west(), {xi});
		View<1> slice_e = pa.getGhostSliceOn(Side<2>::east(), {xi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			CHECK(&pa[{-1 - xi, yi}] == &slice_w[{yi}]);
			CHECK(&pa[{nx + xi, yi}] == &slice_e[{yi}]);
		}
	}

	for (unsigned char yi = 0; yi < num_ghost; yi++) {
		INFO("yi: " << yi);
		View<1> slice_s = pa.getGhostSliceOn(Side<2>::south(), {yi});
		View<1> slice_n = pa.getGhostSliceOn(Side<2>::north(), {yi});
		for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&pa[{xi, -1 - yi}] == &slice_s[{xi}]);
			CHECK(&pa[{xi, ny + yi}] == &slice_n[{xi}]);
		}
	}
}
TEST_CASE("ComponentArray<3> getGhostSliceOn<2>", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	const ComponentArray<3> pa({nx, ny, nz}, num_ghost);

	for (unsigned char xi = 0; xi < num_ghost; xi++) {
		INFO("xi: " << xi);
		View<2> slice_w = pa.getGhostSliceOn(Side<3>::west(), {xi});
		View<2> slice_e = pa.getGhostSliceOn(Side<3>::east(), {xi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{-1 - xi, yi, zi}] == &slice_w[{yi, zi}]);
				CHECK(&pa[{nx + xi, yi, zi}] == &slice_e[{yi, zi}]);
			}
		}
	}

	for (unsigned char yi = 0; yi < num_ghost; yi++) {
		INFO("yi: " << yi);
		View<2> slice_s = pa.getGhostSliceOn(Side<3>::south(), {yi});
		View<2> slice_n = pa.getGhostSliceOn(Side<3>::north(), {yi});
		for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
			INFO("zi: " << zi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, -1 - yi, zi}] == &slice_s[{xi, zi}]);
				CHECK(&pa[{xi, ny + yi, zi}] == &slice_n[{xi, zi}]);
			}
		}
	}

	for (unsigned char zi = 0; zi < num_ghost; zi++) {
		INFO("zi: " << zi);
		View<2> slice_b = pa.getGhostSliceOn(Side<3>::bottom(), {zi});
		View<2> slice_t = pa.getGhostSliceOn(Side<3>::top(), {zi});
		for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
			INFO("yi: " << yi);
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, yi, -1 - zi}] == &slice_b[{xi, yi}]);
				CHECK(&pa[{xi, yi, nz + zi}] == &slice_t[{xi, yi}]);
			}
		}
	}
}
TEST_CASE("ComponentArray<3> getGhostSliceOn<1>", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<3> pa({nx, ny, nz}, num_ghost);

	for (unsigned char zi = 0; zi < num_ghost; zi++) {
		INFO("zi: " << zi);
		for (unsigned char yi = 0; yi < num_ghost; yi++) {
			INFO("yi: " << yi);
			View<1> slice_bs = pa.getGhostSliceOn(Edge::bs(), {yi, zi});
			View<1> slice_tn = pa.getGhostSliceOn(Edge::tn(), {yi, zi});
			View<1> slice_bn = pa.getGhostSliceOn(Edge::bn(), {yi, zi});
			View<1> slice_ts = pa.getGhostSliceOn(Edge::ts(), {yi, zi});
			for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{xi, -1 - yi, -1 - zi}] == &slice_bs[{xi}]);
				CHECK(&pa[{xi, ny + yi, nz + zi}] == &slice_tn[{xi}]);
				CHECK(&pa[{xi, ny + yi, -1 - zi}] == &slice_bn[{xi}]);
				CHECK(&pa[{xi, -1 - yi, nz + zi}] == &slice_ts[{xi}]);
			}
		}
	}

	for (unsigned char zi = 0; zi < num_ghost; zi++) {
		INFO("zi: " << zi);
		for (unsigned char xi = 0; xi < num_ghost; xi++) {
			INFO("xi: " << xi);
			View<1> slice_bw = pa.getGhostSliceOn(Edge::bw(), {xi, zi});
			View<1> slice_te = pa.getGhostSliceOn(Edge::te(), {xi, zi});
			View<1> slice_be = pa.getGhostSliceOn(Edge::be(), {xi, zi});
			View<1> slice_tw = pa.getGhostSliceOn(Edge::tw(), {xi, zi});
			for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
				INFO("yi: " << yi);
				CHECK(&pa[{-1 - xi, yi, -1 - zi}] == &slice_bw[{yi}]);
				CHECK(&pa[{nx + xi, yi, nz + zi}] == &slice_te[{yi}]);
				CHECK(&pa[{nx + xi, yi, -1 - zi}] == &slice_be[{yi}]);
				CHECK(&pa[{-1 - xi, yi, nz + zi}] == &slice_tw[{yi}]);
			}
		}
	}

	for (unsigned char yi = 0; yi < num_ghost; yi++) {
		INFO("yi: " << yi);
		for (unsigned char xi = 0; xi < num_ghost; xi++) {
			INFO("xi: " << xi);
			View<1> slice_sw = pa.getGhostSliceOn(Edge::sw(), {xi, yi});
			View<1> slice_ne = pa.getGhostSliceOn(Edge::ne(), {xi, yi});
			View<1> slice_se = pa.getGhostSliceOn(Edge::se(), {xi, yi});
			View<1> slice_nw = pa.getGhostSliceOn(Edge::nw(), {xi, yi});
			for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
				INFO("zi: " << zi);
				CHECK(&pa[{-1 - xi, -1 - yi, zi}] == &slice_sw[{zi}]);
				CHECK(&pa[{nx + xi, ny + yi, zi}] == &slice_ne[{zi}]);
				CHECK(&pa[{nx + xi, -1 - yi, zi}] == &slice_se[{zi}]);
				CHECK(&pa[{-1 - xi, ny + yi, zi}] == &slice_nw[{zi}]);
			}
		}
	}
}
TEST_CASE("ComponentArray<3> getGhostSliceOn<0>", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto nz        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<3> pa({nx, ny, nz}, num_ghost);

	for (unsigned char zi = 0; zi < num_ghost; zi++) {
		INFO("zi: " << zi);
		for (unsigned char yi = 0; yi < num_ghost; yi++) {
			INFO("yi: " << yi);
			for (unsigned char xi = 0; xi < num_ghost; xi++) {
				INFO("xi: " << xi);
				CHECK(&pa[{-1 - xi, -1 - yi, -1 - zi}] == &(pa.getGhostSliceOn(Corner<3>::bsw(), {xi, yi, zi}))[{}]);
				CHECK(&pa[{nx + xi, -1 - yi, -1 - zi}] == &(pa.getGhostSliceOn(Corner<3>::bse(), {xi, yi, zi})[{}]));
				CHECK(&pa[{-1 - xi, ny + yi, -1 - zi}] == &(pa.getGhostSliceOn(Corner<3>::bnw(), {xi, yi, zi})[{}]));
				CHECK(&pa[{nx + xi, ny + yi, -1 - zi}] == &(pa.getGhostSliceOn(Corner<3>::bne(), {xi, yi, zi})[{}]));
				CHECK(&pa[{-1 - xi, -1 - yi, nz + zi}] == &(pa.getGhostSliceOn(Corner<3>::tsw(), {xi, yi, zi})[{}]));
				CHECK(&pa[{nx + xi, -1 - yi, nz + zi}] == &(pa.getGhostSliceOn(Corner<3>::tse(), {xi, yi, zi})[{}]));
				CHECK(&pa[{-1 - xi, ny + yi, nz + zi}] == &(pa.getGhostSliceOn(Corner<3>::tnw(), {xi, yi, zi})[{}]));
				CHECK(&pa[{nx + xi, ny + yi, nz + zi}] == &(pa.getGhostSliceOn(Corner<3>::tne(), {xi, yi, zi})[{}]));
			}
		}
	}
}
TEST_CASE("ComponentArray<2> getGhostSliceOn<0>", "[ComponentArray]")
{
	auto nx        = GENERATE(2, 3);
	auto ny        = GENERATE(2, 3);
	auto num_ghost = GENERATE(0, 1, 2);

	ComponentArray<2> pa({nx, ny}, num_ghost);

	for (unsigned char yi = 0; yi < num_ghost; yi++) {
		INFO("yi: " << yi);
		for (unsigned char xi = 0; xi < num_ghost; xi++) {
			INFO("xi: " << xi);
			CHECK(&pa[{-xi - 1, -yi - 1}] == &(pa.getGhostSliceOn(Corner<2>::sw(), {xi, yi})[{}]));
			CHECK(&pa[{nx + xi, -yi - 1}] == &(pa.getGhostSliceOn(Corner<2>::se(), {xi, yi})[{}]));
			CHECK(&pa[{-xi - 1, ny + yi}] == &(pa.getGhostSliceOn(Corner<2>::nw(), {xi, yi})[{}]));
			CHECK(&pa[{nx + xi, ny + yi}] == &(pa.getGhostSliceOn(Corner<2>::ne(), {xi, yi})[{}]));
		}
	}
}