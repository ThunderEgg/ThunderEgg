#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/ValVector.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;
using namespace Catch;

TEST_CASE("DomainTools::GetRealCoord 1D", "[getRealCoord][DomainTools]")
{
	std::shared_ptr<PatchInfo<1>> pinfo(new PatchInfo<1>());

	auto nx      = GENERATE(1, 2, 3, 10, 16, 17);
	auto startx  = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	auto lengthx = GENERATE(1.0, 0.23, 3.0);

	pinfo->ns       = {nx};
	pinfo->spacings = {lengthx / nx};
	pinfo->starts   = {startx};

	for (int coordx = -1; coordx <= nx; coordx++) {
		array<int, 1>    coord = {coordx};
		array<double, 1> expected;
		if (coordx == -1) {
			expected[0] = startx;
		} else if (coordx == nx) {
			expected[0] = startx + lengthx;
		} else {
			expected[0] = startx + (0.5 + coordx) * lengthx / nx;
		}

		array<double, 1> result;
		DomainTools::GetRealCoord<1>(pinfo, coord, result);

		CHECK(result[0] + 100 == Approx(expected[0] + 100));
	}
}
TEST_CASE("DomainTools::GetRealCoord 2D", "[getRealCoord][DomainTools]")
{
	std::shared_ptr<PatchInfo<2>> pinfo(new PatchInfo<2>());

	auto nx      = 2;
	auto ny      = GENERATE(1, 2, 3);
	auto startx  = 0.0;
	auto starty  = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	auto lengthx = 1.0;
	auto lengthy = GENERATE(1.0, 0.23, 3.0);

	pinfo->ns       = {nx, ny};
	pinfo->spacings = {lengthx / nx, lengthy / ny};
	pinfo->starts   = {startx, starty};

	for (int coordy = -1; coordy <= ny; coordy++) {
		for (int coordx = -1; coordx <= nx; coordx++) {
			array<double, 2> expected;
			if (coordx == -1) {
				expected[0] = startx;
			} else if (coordx == nx) {
				expected[0] = startx + lengthx;
			} else {
				expected[0] = startx + (0.5 + coordx) * lengthx / nx;
			}
			if (coordy == -1) {
				expected[1] = starty;
			} else if (coordy == ny) {
				expected[1] = starty + lengthy;
			} else {
				expected[1] = starty + (0.5 + coordy) * lengthy / ny;
			}

			array<double, 2> result;
			array<int, 2>    coord = {coordx, coordy};
			DomainTools::GetRealCoord<2>(pinfo, coord, result);

			CHECK(result[0] + 100 == Approx(expected[0] + 100));
			CHECK(result[1] + 100 == Approx(expected[1] + 100));
		}
	}
}
TEST_CASE("DomainTools::GetRealCoord 3D", "[getRealCoord][DomainTools]")
{
	std::shared_ptr<PatchInfo<3>> pinfo(new PatchInfo<3>());

	auto nx      = 2;
	auto ny      = 2;
	auto nz      = GENERATE(1, 2, 3);
	auto startx  = 0.0;
	auto starty  = 0.0;
	auto startz  = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	auto lengthx = 1.0;
	auto lengthy = 1.0;
	auto lengthz = GENERATE(1.0, 0.23, 3.0);

	pinfo->ns       = {nx, ny, nz};
	pinfo->spacings = {lengthx / nx, lengthy / ny, lengthz / nz};
	pinfo->starts   = {startx, starty, startz};

	for (int coordz = -1; coordz <= nz; coordz++) {
		for (int coordy = -1; coordy <= ny; coordy++) {
			for (int coordx = -1; coordx <= nx; coordx++) {
				array<double, 3> expected;
				if (coordx == -1) {
					expected[0] = startx;
				} else if (coordx == nx) {
					expected[0] = startx + lengthx;
				} else {
					expected[0] = startx + (0.5 + coordx) * lengthx / nx;
				}
				if (coordy == -1) {
					expected[1] = starty;
				} else if (coordy == ny) {
					expected[1] = starty + lengthy;
				} else {
					expected[1] = starty + (0.5 + coordy) * lengthy / ny;
				}
				if (coordz == -1) {
					expected[2] = startz;
				} else if (coordz == nz) {
					expected[2] = startz + lengthz;
				} else {
					expected[2] = startz + (0.5 + coordz) * lengthz / nz;
				}

				array<double, 3> result;
				array<int, 3>    coord = {coordx, coordy, coordz};
				DomainTools::GetRealCoord<3>(pinfo, coord, result);

				CHECK(result[0] + 100 == Approx(expected[0] + 100));
				CHECK(result[1] + 100 == Approx(expected[1] + 100));
				CHECK(result[2] + 100 == Approx(expected[2] + 100));
			}
		}
	}
}
TEST_CASE("DomainTools::getRealCoordGhost 1D", "[getRealCoordGhost][DomainTools]")
{
	std::shared_ptr<PatchInfo<1>> pinfo(new PatchInfo<1>());

	auto nx      = GENERATE(1, 2, 3, 10, 16, 17);
	auto startx  = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	auto lengthx = GENERATE(1.0, 0.23, 3.0);

	pinfo->ns       = {nx};
	pinfo->spacings = {lengthx / nx};
	pinfo->starts   = {startx};

	for (int coordx = -1; coordx <= nx; coordx++) {
		array<int, 1>    coord = {coordx};
		array<double, 1> expected;
		expected[0] = startx + (0.5 + coordx) * lengthx / nx;

		array<double, 1> result;
		DomainTools::GetRealCoordGhost<1>(pinfo, coord, result);

		CHECK(result[0] + 100 == Approx(expected[0] + 100));
	}
}
TEST_CASE("DomainTools::getRealCoordGhost 2D", "[getRealCoordGhost][DomainTools]")
{
	std::shared_ptr<PatchInfo<2>> pinfo(new PatchInfo<2>());

	auto nx      = 2;
	auto ny      = GENERATE(1, 2, 3);
	auto startx  = 0.0;
	auto starty  = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	auto lengthx = 1.0;
	auto lengthy = GENERATE(1.0, 0.23, 3.0);

	pinfo->ns       = {nx, ny};
	pinfo->spacings = {lengthx / nx, lengthy / ny};
	pinfo->starts   = {startx, starty};

	for (int coordy = -1; coordy <= ny; coordy++) {
		for (int coordx = -1; coordx <= nx; coordx++) {
			array<double, 2> expected;
			expected[0] = startx + (0.5 + coordx) * lengthx / nx;
			expected[1] = starty + (0.5 + coordy) * lengthy / ny;

			array<double, 2> result;
			array<int, 2>    coord = {coordx, coordy};
			DomainTools::GetRealCoordGhost<2>(pinfo, coord, result);

			CHECK(result[0] + 100 == Approx(expected[0] + 100));
			CHECK(result[1] + 100 == Approx(expected[1] + 100));
		}
	}
}
TEST_CASE("DomainTools::getRealCoordGhost 3D", "[getRealCoordGhost][DomainTools]")
{
	std::shared_ptr<PatchInfo<3>> pinfo(new PatchInfo<3>());

	auto nx      = 2;
	auto ny      = 2;
	auto nz      = GENERATE(1, 2, 3);
	auto startx  = 0.0;
	auto starty  = 0.0;
	auto startz  = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	auto lengthx = 1.0;
	auto lengthy = 1.0;
	auto lengthz = GENERATE(1.0, 0.23, 3.0);

	pinfo->ns       = {nx, ny, nz};
	pinfo->spacings = {lengthx / nx, lengthy / ny, lengthz / nz};
	pinfo->starts   = {startx, starty, startz};

	for (int coordz = -1; coordz <= nz; coordz++) {
		for (int coordy = -1; coordy <= ny; coordy++) {
			for (int coordx = -1; coordx <= nx; coordx++) {
				array<double, 3> expected;
				expected[0] = startx + (0.5 + coordx) * lengthx / nx;
				expected[1] = starty + (0.5 + coordy) * lengthy / ny;
				expected[2] = startz + (0.5 + coordz) * lengthz / nz;

				array<double, 3> result;
				array<int, 3>    coord = {coordx, coordy, coordz};
				DomainTools::GetRealCoordGhost<3>(pinfo, coord, result);

				CHECK(result[0] + 100 == Approx(expected[0] + 100));
				CHECK(result[1] + 100 == Approx(expected[1] + 100));
				CHECK(result[2] + 100 == Approx(expected[2] + 100));
			}
		}
	}
}
TEST_CASE("DomainTools::GetRealCoordBound 1D", "[DomainTools]")
{
	std::shared_ptr<PatchInfo<1>> pinfo(new PatchInfo<1>());

	auto nx      = GENERATE(1, 2, 3, 10, 16, 17);
	auto startx  = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	auto lengthx = GENERATE(1.0, 0.23, 3.0);

	pinfo->ns       = {nx};
	pinfo->spacings = {lengthx / nx};
	pinfo->starts   = {startx};

	std::array<int, 0>    coord;
	std::array<double, 1> result;
	DomainTools::GetRealCoordBound<1>(pinfo, coord, Side<1>::west(), result);
	CHECK(result[0] + 100 == Approx(startx + 100));
	DomainTools::GetRealCoordBound<1>(pinfo, coord, Side<1>::east(), result);
	CHECK(result[0] + 100 == Approx(startx + lengthx + 100));
}
TEST_CASE("DomainTools::GetRealCoordBound 2D", "[getRealCoord][DomainTools]")
{
	std::shared_ptr<PatchInfo<2>> pinfo(new PatchInfo<2>());

	auto nx      = 2;
	auto ny      = GENERATE(1, 2, 3);
	auto startx  = 0.0;
	auto starty  = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	auto lengthx = 1.0;
	auto lengthy = GENERATE(1.0, 0.23, 3.0);

	pinfo->ns       = {nx, ny};
	pinfo->spacings = {lengthx / nx, lengthy / ny};
	pinfo->starts   = {startx, starty};

	for (int coordx = -1; coordx <= nx; coordx++) {
		array<double, 2> expected;
		if (coordx == -1) {
			expected[0] = startx;
		} else if (coordx == nx) {
			expected[0] = startx + lengthx;
		} else {
			expected[0] = startx + (0.5 + coordx) * lengthx / nx;
		}

		array<double, 2> result;
		array<int, 1>    coord = {coordx};
		DomainTools::GetRealCoordBound<2>(pinfo, coord, Side<2>::south(), result);

		CHECK(result[0] + 100 == Approx(expected[0] + 100));
		CHECK(result[1] + 100 == Approx(starty + 100));

		DomainTools::GetRealCoordBound<2>(pinfo, coord, Side<2>::north(), result);
		CHECK(result[0] + 100 == Approx(expected[0] + 100));
		CHECK(result[1] + 100 == Approx(starty + lengthy + 100));
	}
	for (int coordy = -1; coordy <= ny; coordy++) {
		array<double, 2> expected;
		if (coordy == -1) {
			expected[1] = starty;
		} else if (coordy == ny) {
			expected[1] = starty + lengthy;
		} else {
			expected[1] = starty + (0.5 + coordy) * lengthy / ny;
		}

		array<double, 2> result;
		array<int, 1>    coord = {coordy};
		DomainTools::GetRealCoordBound<2>(pinfo, coord, Side<2>::west(), result);

		CHECK(result[0] + 100 == Approx(startx + 100));
		CHECK(result[1] + 100 == Approx(expected[1] + 100));

		DomainTools::GetRealCoordBound<2>(pinfo, coord, Side<2>::east(), result);

		CHECK(result[0] + 100 == Approx(startx + lengthx + 100));
		CHECK(result[1] + 100 == Approx(expected[1] + 100));
	}
}
TEST_CASE("DomainTools::setValues 1D g=x", "[DomainTools]")
{
	auto f = [](const std::array<double, 1> coord) { return coord[0]; };

	vector<shared_ptr<PatchInfo<1>>> pinfos(1);

	auto nx        = GENERATE(1, 2, 10, 13);
	auto startx    = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	auto spacingx  = GENERATE(0.01, 1.0, 3.14);
	auto num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<1>());
	pinfos[0]->id              = 0;
	pinfos[0]->ns              = {nx};
	pinfos[0]->spacings        = {spacingx};
	pinfos[0]->starts          = {startx};
	pinfos[0]->num_ghost_cells = num_ghost;
	shared_ptr<Domain<1>> d(new Domain<1>(0, {nx}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<1>> vec(
	new ValVector<1>(MPI_COMM_WORLD, pinfos[0]->ns, num_ghost, 1, 1));

	DomainTools::SetValues<1>(d, vec, f);
	auto ld = vec->getLocalData(0, 0);
	nested_loop<1>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 1> coord) {
		if (coord[0] < 0 || coord[0] >= nx) {
			CHECK(ld[coord] + 100 == Approx(0.0 + 100));
		} else {
			std::array<double, 1> real_coord;
			DomainTools::GetRealCoord<1>(pinfos[0], coord, real_coord);
			CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
		}
	});
}
TEST_CASE("DomainTools::setValues 1D f=x**2", "[DomainTools]")
{
	auto f = [](const std::array<double, 1> coord) { return coord[0] * coord[0]; };

	vector<shared_ptr<PatchInfo<1>>> pinfos(1);

	auto nx        = GENERATE(1, 2, 10, 13);
	auto startx    = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	auto spacingx  = GENERATE(0.01, 1.0, 3.14);
	auto num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<1>());
	pinfos[0]->id              = 0;
	pinfos[0]->ns              = {nx};
	pinfos[0]->spacings        = {spacingx};
	pinfos[0]->starts          = {startx};
	pinfos[0]->num_ghost_cells = num_ghost;
	shared_ptr<Domain<1>> d(new Domain<1>(0, {nx}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<1>> vec(
	new ValVector<1>(MPI_COMM_WORLD, pinfos[0]->ns, num_ghost, 1, 1));

	DomainTools::SetValues<1>(d, vec, f);
	auto ld = vec->getLocalData(0, 0);
	nested_loop<1>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 1> coord) {
		if (coord[0] < 0 || coord[0] >= nx) {
			CHECK(ld[coord] + 100 == Approx(0.0 + 100));
		} else {
			std::array<double, 1> real_coord;
			DomainTools::GetRealCoord<1>(pinfos[0], coord, real_coord);
			CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
		}
	});
}
TEST_CASE("DomainTools::setValues 2D f=x+y", "[DomainTools]")
{
	auto f = [](const std::array<double, 2> coord) { return coord[0] + coord[1]; };

	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    nx        = 3;
	auto   ny        = GENERATE(1, 2, 10, 13);
	double startx    = 0.0;
	auto   starty    = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	double spacingx  = 0.1;
	auto   spacingy  = GENERATE(0.01, 1.0, 3.14);
	auto   num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id              = 0;
	pinfos[0]->ns              = {nx, ny};
	pinfos[0]->spacings        = {spacingx, spacingy};
	pinfos[0]->starts          = {startx, starty};
	pinfos[0]->num_ghost_cells = num_ghost;
	shared_ptr<Domain<2>> d(new Domain<2>(0, {nx, ny}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<2>> vec(
	new ValVector<2>(MPI_COMM_WORLD, pinfos[0]->ns, num_ghost, 1, 1));

	DomainTools::SetValues<2>(d, vec, f);
	auto ld = vec->getLocalData(0, 0);
	nested_loop<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
		if (coord[0] < 0 || coord[0] >= nx || coord[1] < 0 || coord[1] >= ny) {
			CHECK(ld[coord] + 100 == Approx(0.0 + 100));
		} else {
			std::array<double, 2> real_coord;
			DomainTools::GetRealCoord<2>(pinfos[0], coord, real_coord);
			CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
		}
	});
}
TEST_CASE("DomainTools::setValues 2D f=x+y,g=x*y", "[DomainTools]")
{
	auto f = [](const std::array<double, 2> coord) { return coord[0] + coord[1]; };
	auto g = [](const std::array<double, 2> coord) { return coord[0] * coord[1]; };

	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    nx        = 3;
	auto   ny        = GENERATE(1, 2, 10, 13);
	double startx    = 0.0;
	auto   starty    = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	double spacingx  = 0.1;
	auto   spacingy  = GENERATE(0.01, 1.0, 3.14);
	auto   num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id              = 0;
	pinfos[0]->ns              = {nx, ny};
	pinfos[0]->spacings        = {spacingx, spacingy};
	pinfos[0]->starts          = {startx, starty};
	pinfos[0]->num_ghost_cells = num_ghost;
	shared_ptr<Domain<2>> d(new Domain<2>(0, {nx, ny}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<2>> vec(
	new ValVector<2>(MPI_COMM_WORLD, pinfos[0]->ns, num_ghost, 2, 1));

	DomainTools::SetValues<2>(d, vec, f, g);
	auto ld = vec->getLocalData(0, 0);
	nested_loop<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
		if (coord[0] < 0 || coord[0] >= nx || coord[1] < 0 || coord[1] >= ny) {
			CHECK(ld[coord] + 100 == Approx(0.0 + 100));
		} else {
			std::array<double, 2> real_coord;
			DomainTools::GetRealCoord<2>(pinfos[0], coord, real_coord);
			CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
		}
	});
	auto ld2 = vec->getLocalData(1, 0);
	nested_loop<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
		if (coord[0] < 0 || coord[0] >= nx || coord[1] < 0 || coord[1] >= ny) {
			CHECK(ld2[coord] + 100 == Approx(0.0 + 100));
		} else {
			std::array<double, 2> real_coord;
			DomainTools::GetRealCoord<2>(pinfos[0], coord, real_coord);
			CHECK(ld2[coord] + 100 == Approx(g(real_coord) + 100));
		}
	});
}
TEST_CASE("DomainTools::setValues 2D f=x*y", "[DomainTools]")
{
	auto f = [](const std::array<double, 2> coord) { return coord[0] * coord[1]; };

	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    nx        = 3;
	auto   ny        = GENERATE(1, 2, 10, 13);
	double startx    = 0.0;
	auto   starty    = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	double spacingx  = 0.1;
	auto   spacingy  = GENERATE(0.01, 1.0, 3.14);
	auto   num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id              = 0;
	pinfos[0]->ns              = {nx, ny};
	pinfos[0]->spacings        = {spacingx, spacingy};
	pinfos[0]->starts          = {startx, starty};
	pinfos[0]->num_ghost_cells = num_ghost;
	shared_ptr<Domain<2>> d(new Domain<2>(0, {nx, ny}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<2>> vec(
	new ValVector<2>(MPI_COMM_WORLD, pinfos[0]->ns, num_ghost, 1, 1));

	DomainTools::SetValues<2>(d, vec, f);
	auto ld = vec->getLocalData(0, 0);
	nested_loop<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
		if (coord[0] < 0 || coord[0] >= nx || coord[1] < 0 || coord[1] >= ny) {
			CHECK(ld[coord] + 100 == Approx(0.0 + 100));
		} else {
			std::array<double, 2> real_coord;
			DomainTools::GetRealCoord<2>(pinfos[0], coord, real_coord);
			CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
		}
	});
}
TEST_CASE("DomainTools::setValuesWithGhost 1D f=x", "[DomainTools]")
{
	auto f = [](const std::array<double, 1> coord) { return coord[0]; };

	vector<shared_ptr<PatchInfo<1>>> pinfos(1);

	auto nx        = GENERATE(1, 2, 10, 13);
	auto startx    = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	auto spacingx  = GENERATE(0.01, 1.0, 3.14);
	auto num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<1>());
	pinfos[0]->id              = 0;
	pinfos[0]->ns              = {nx};
	pinfos[0]->spacings        = {spacingx};
	pinfos[0]->starts          = {startx};
	pinfos[0]->num_ghost_cells = num_ghost;
	shared_ptr<Domain<1>> d(new Domain<1>(0, {nx}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<1>> vec(
	new ValVector<1>(MPI_COMM_WORLD, pinfos[0]->ns, num_ghost, 1, 1));

	DomainTools::SetValuesWithGhost<1>(d, vec, f);
	auto ld = vec->getLocalData(0, 0);
	nested_loop<1>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 1> coord) {
		std::array<double, 1> real_coord;
		DomainTools::GetRealCoordGhost<1>(pinfos[0], coord, real_coord);
		CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
	});
}
TEST_CASE("DomainTools::setValuesWithGhost 1D f=x**2", "[DomainTools]")
{
	auto f = [](const std::array<double, 1> coord) { return coord[0] * coord[0]; };

	vector<shared_ptr<PatchInfo<1>>> pinfos(1);

	auto nx        = GENERATE(1, 2, 10, 13);
	auto startx    = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	auto spacingx  = GENERATE(0.01, 1.0, 3.14);
	auto num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<1>());
	pinfos[0]->id              = 0;
	pinfos[0]->ns              = {nx};
	pinfos[0]->spacings        = {spacingx};
	pinfos[0]->starts          = {startx};
	pinfos[0]->num_ghost_cells = num_ghost;
	shared_ptr<Domain<1>> d(new Domain<1>(0, {nx}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<1>> vec(
	new ValVector<1>(MPI_COMM_WORLD, pinfos[0]->ns, num_ghost, 1, 1));

	DomainTools::SetValuesWithGhost<1>(d, vec, f);
	auto ld = vec->getLocalData(0, 0);
	nested_loop<1>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 1> coord) {
		std::array<double, 1> real_coord;
		DomainTools::GetRealCoordGhost<1>(pinfos[0], coord, real_coord);
		CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
	});
}
TEST_CASE("DomainTools::setValuesWithGhost 2D f=x+y", "[DomainTools]")
{
	auto f = [](const std::array<double, 2> coord) { return coord[0] + coord[1]; };

	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    nx        = 3;
	auto   ny        = GENERATE(1, 2, 10, 13);
	double startx    = 0.0;
	auto   starty    = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	double spacingx  = 0.1;
	auto   spacingy  = GENERATE(0.01, 1.0, 3.14);
	auto   num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id              = 0;
	pinfos[0]->ns              = {nx, ny};
	pinfos[0]->spacings        = {spacingx, spacingy};
	pinfos[0]->starts          = {startx, starty};
	pinfos[0]->num_ghost_cells = num_ghost;
	shared_ptr<Domain<2>> d(new Domain<2>(0, {nx, ny}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<2>> vec(
	new ValVector<2>(MPI_COMM_WORLD, pinfos[0]->ns, num_ghost, 1, 1));

	DomainTools::SetValuesWithGhost<2>(d, vec, f);
	auto ld = vec->getLocalData(0, 0);
	nested_loop<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
		std::array<double, 2> real_coord;
		DomainTools::GetRealCoordGhost<2>(pinfos[0], coord, real_coord);
		CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
	});
}
TEST_CASE("DomainTools::setValuesWithGhost 2D f=x+y,g=x*y", "[DomainTools]")
{
	auto f = [](const std::array<double, 2> coord) { return coord[0] + coord[1]; };
	auto g = [](const std::array<double, 2> coord) { return coord[0] * coord[1]; };

	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    nx        = 3;
	auto   ny        = GENERATE(1, 2, 10, 13);
	double startx    = 0.0;
	auto   starty    = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	double spacingx  = 0.1;
	auto   spacingy  = GENERATE(0.01, 1.0, 3.14);
	auto   num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id              = 0;
	pinfos[0]->ns              = {nx, ny};
	pinfos[0]->spacings        = {spacingx, spacingy};
	pinfos[0]->starts          = {startx, starty};
	pinfos[0]->num_ghost_cells = num_ghost;
	shared_ptr<Domain<2>> d(new Domain<2>(0, {nx, ny}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<2>> vec(
	new ValVector<2>(MPI_COMM_WORLD, pinfos[0]->ns, num_ghost, 2, 1));

	DomainTools::SetValuesWithGhost<2>(d, vec, f, g);
	auto ld = vec->getLocalData(0, 0);
	nested_loop<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
		std::array<double, 2> real_coord;
		DomainTools::GetRealCoordGhost<2>(pinfos[0], coord, real_coord);
		CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
	});
	auto ld2 = vec->getLocalData(1, 0);
	nested_loop<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
		std::array<double, 2> real_coord;
		DomainTools::GetRealCoordGhost<2>(pinfos[0], coord, real_coord);
		CHECK(ld2[coord] + 100 == Approx(g(real_coord) + 100));
	});
}
TEST_CASE("DomainTools::setValuesWithGhost 2D throws when too many functions are given",
          "[DomainTools]")
{
	auto f = [](const std::array<double, 2> coord) { return coord[0] + coord[1]; };
	auto g = [](const std::array<double, 2> coord) { return coord[0] * coord[1]; };

	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    nx        = 3;
	auto   ny        = GENERATE(1, 2, 10, 13);
	double startx    = 0.0;
	auto   starty    = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	double spacingx  = 0.1;
	auto   spacingy  = GENERATE(0.01, 1.0, 3.14);
	auto   num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id              = 0;
	pinfos[0]->ns              = {nx, ny};
	pinfos[0]->spacings        = {spacingx, spacingy};
	pinfos[0]->starts          = {startx, starty};
	pinfos[0]->num_ghost_cells = num_ghost;
	shared_ptr<Domain<2>> d(new Domain<2>(0, {nx, ny}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<2>> vec(
	new ValVector<2>(MPI_COMM_WORLD, pinfos[0]->ns, num_ghost, 1, 1));

	CHECK_THROWS_AS(DomainTools::SetValuesWithGhost<2>(d, vec, f, g), RuntimeError);
}
TEST_CASE("DomainTools::setValues 2D throws when too many functions are given", "[DomainTools]")
{
	auto f = [](const std::array<double, 2> coord) { return coord[0] + coord[1]; };
	auto g = [](const std::array<double, 2> coord) { return coord[0] * coord[1]; };

	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    nx        = 3;
	auto   ny        = GENERATE(1, 2, 10, 13);
	double startx    = 0.0;
	auto   starty    = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	double spacingx  = 0.1;
	auto   spacingy  = GENERATE(0.01, 1.0, 3.14);
	auto   num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id              = 0;
	pinfos[0]->ns              = {nx, ny};
	pinfos[0]->spacings        = {spacingx, spacingy};
	pinfos[0]->starts          = {startx, starty};
	pinfos[0]->num_ghost_cells = num_ghost;
	shared_ptr<Domain<2>> d(new Domain<2>(0, {nx, ny}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<2>> vec(
	new ValVector<2>(MPI_COMM_WORLD, pinfos[0]->ns, num_ghost, 1, 1));

	CHECK_THROWS_AS(DomainTools::SetValues<2>(d, vec, f, g), RuntimeError);
}
TEST_CASE("DomainTools::setValuesWithGhost 2D f=x*y", "[DomainTools]")
{
	auto f = [](const std::array<double, 2> coord) { return coord[0] * coord[1]; };

	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    nx        = 3;
	auto   ny        = GENERATE(1, 2, 10, 13);
	double startx    = 0.0;
	auto   starty    = GENERATE(0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0);
	double spacingx  = 0.1;
	auto   spacingy  = GENERATE(0.01, 1.0, 3.14);
	auto   num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id              = 0;
	pinfos[0]->ns              = {nx, ny};
	pinfos[0]->spacings        = {spacingx, spacingy};
	pinfos[0]->starts          = {startx, starty};
	pinfos[0]->num_ghost_cells = num_ghost;
	shared_ptr<Domain<2>> d(new Domain<2>(0, {nx, ny}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<2>> vec(
	new ValVector<2>(MPI_COMM_WORLD, pinfos[0]->ns, num_ghost, 1, 1));

	DomainTools::SetValuesWithGhost<2>(d, vec, f);
	auto ld = vec->getLocalData(0, 0);
	nested_loop<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
		std::array<double, 2> real_coord;
		DomainTools::GetRealCoordGhost<2>(pinfos[0], coord, real_coord);
		CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
	});
}
// TODO set boundary vector