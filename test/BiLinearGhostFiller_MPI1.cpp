#include "utils/DomainReader.h"
#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/ValVectorGenerator.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;
using namespace Catch;

constexpr auto single_mesh_file  = "mesh_inputs/2d_uniform_2x2_mpi1.json";
constexpr auto refined_mesh_file = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json";
constexpr auto cross_mesh_file   = "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json";
TEST_CASE("exchange uniform 2D quad BiLinearGhostFiller", "[BiLinearGhostFiller]")
{
	auto   nx        = GENERATE(1, 2, 10, 13);
	auto   ny        = GENERATE(1, 2, 10, 13);
	double startx    = 0;
	double starty    = 0;
	double lengthx   = 1;
	double lengthy   = 1;
	int    num_ghost = 1;

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 1 + ((x + 0.3) * y);
	};

	/*
	 *  2 | 3
	 *  --+---
	 *  0 | 1
	 */
	vector<PatchInfo<2>> pinfos(4);

	pinfos[0].id              = 0;
	pinfos[0].ns              = {nx, ny};
	pinfos[0].spacings        = {lengthx / nx, lengthy / ny};
	pinfos[0].starts          = {startx, starty};
	pinfos[0].num_ghost_cells = num_ghost;
	pinfos[0].nbr_info[Side<2>::east().getIndex()].reset(new NormalNbrInfo<2>(1));
	pinfos[0].nbr_info[Side<2>::north().getIndex()].reset(new NormalNbrInfo<2>(2));

	pinfos[1].id              = 1;
	pinfos[1].ns              = {nx, ny};
	pinfos[1].spacings        = {lengthx / nx, lengthy / ny};
	pinfos[1].starts          = {startx + lengthx, starty};
	pinfos[1].num_ghost_cells = num_ghost;
	pinfos[1].nbr_info[Side<2>::west().getIndex()].reset(new NormalNbrInfo<2>(0));
	pinfos[1].nbr_info[Side<2>::north().getIndex()].reset(new NormalNbrInfo<2>(3));

	pinfos[2].id              = 2;
	pinfos[2].ns              = {nx, ny};
	pinfos[2].spacings        = {lengthx / nx, lengthy / ny};
	pinfos[2].starts          = {startx, starty + lengthy};
	pinfos[2].num_ghost_cells = num_ghost;
	pinfos[2].nbr_info[Side<2>::east().getIndex()].reset(new NormalNbrInfo<2>(3));
	pinfos[2].nbr_info[Side<2>::south().getIndex()].reset(new NormalNbrInfo<2>(0));

	pinfos[3].id              = 3;
	pinfos[3].ns              = {nx, ny};
	pinfos[3].spacings        = {lengthx / nx, lengthy / ny};
	pinfos[3].starts          = {startx + lengthx, starty + lengthy};
	pinfos[3].num_ghost_cells = num_ghost;
	pinfos[3].nbr_info[Side<2>::west().getIndex()].reset(new NormalNbrInfo<2>(2));
	pinfos[3].nbr_info[Side<2>::south().getIndex()].reset(new NormalNbrInfo<2>(1));

	shared_ptr<Domain<2>> d(new Domain<2>(0, {nx, ny}, num_ghost, pinfos.begin(), pinfos.end()));

	shared_ptr<ValVector<2>> vec(
	new ValVector<2>(MPI_COMM_WORLD, pinfos[0].ns, num_ghost, 1, 4));

	DomainTools::SetValues<2>(d, vec, f);

	BiLinearGhostFiller blgf(d);
	blgf.fillGhost(vec);

	// patch 1
	{
		// check that center values weren't modified
		auto patch_1 = vec->getLocalData(0, d->getPatchInfoVector()[0].local_index);
		nested_loop<2>(patch_1.getStart(), patch_1.getEnd(), [&](const std::array<int, 2> coord) {
			std::array<double, 2> real_coord;
			DomainTools::GetRealCoord<2>(pinfos[0], coord, real_coord);
			CHECK(patch_1[coord] == f(real_coord));
		});
		// check that west values are not modified
		{
			auto west_ghost = patch_1.getGhostSliceOnSide(Side<2>::west(), 1);
			nested_loop<1>(west_ghost.getStart(), west_ghost.getEnd(),
			               [&](const std::array<int, 1> coord) { CHECK(west_ghost[coord] == 0); });
		}
		// check that east values are correct
		{
			auto east_ghost = patch_1.getGhostSliceOnSide(Side<2>::east(), 1);
			nested_loop<1>(
			east_ghost.getStart(), east_ghost.getEnd(), [&](const std::array<int, 1> coord) {
				std::array<double, 2> real_coord;
				DomainTools::GetRealCoordGhost<2>(pinfos[0], {nx, coord[0]}, real_coord);
				CHECK(east_ghost[coord] == Approx(f(real_coord)));
			});
		}
		// check that south values are not modified
		{
			auto south_ghost = patch_1.getGhostSliceOnSide(Side<2>::south(), 1);
			nested_loop<1>(south_ghost.getStart(), south_ghost.getEnd(),
			               [&](const std::array<int, 1> coord) { CHECK(south_ghost[coord] == 0); });
		}
		// check that north values are correct
		{
			auto north_ghost = patch_1.getGhostSliceOnSide(Side<2>::north(), 1);
			nested_loop<1>(
			north_ghost.getStart(), north_ghost.getEnd(), [&](const std::array<int, 1> coord) {
				std::array<double, 2> real_coord;
				DomainTools::GetRealCoordGhost<2>(pinfos[0], {coord[0], ny}, real_coord);
				CHECK(north_ghost[coord] == Approx(f(real_coord)));
			});
		}
	}
	// patch 2
	{
		// check that center values weren't modified
		auto patch_2 = vec->getLocalData(0, d->getPatchInfoVector()[1].local_index);
		nested_loop<2>(patch_2.getStart(), patch_2.getEnd(), [&](const std::array<int, 2> coord) {
			std::array<double, 2> real_coord;
			DomainTools::GetRealCoord<2>(pinfos[1], coord, real_coord);
			CHECK(patch_2[coord] == f(real_coord));
		});
		// check that west values correct
		{
			auto west_ghost = patch_2.getGhostSliceOnSide(Side<2>::west(), 1);
			nested_loop<1>(
			west_ghost.getStart(), west_ghost.getEnd(), [&](const std::array<int, 1> coord) {
				std::array<double, 2> real_coord;
				DomainTools::GetRealCoordGhost<2>(pinfos[1], {-1, coord[0]}, real_coord);
				CHECK(west_ghost[coord] == Approx(f(real_coord)));
			});
		}
		// check that east values are not modified
		{
			auto east_ghost = patch_2.getGhostSliceOnSide(Side<2>::east(), 1);
			nested_loop<1>(east_ghost.getStart(), east_ghost.getEnd(),
			               [&](const std::array<int, 1> coord) { CHECK(east_ghost[coord] == 0); });
		}
		// check that south values are not modified
		{
			auto south_ghost = patch_2.getGhostSliceOnSide(Side<2>::south(), 1);
			nested_loop<1>(south_ghost.getStart(), south_ghost.getEnd(),
			               [&](const std::array<int, 1> coord) { CHECK(south_ghost[coord] == 0); });
		}
		// check that north values are correct
		{
			auto north_ghost = patch_2.getGhostSliceOnSide(Side<2>::north(), 1);
			nested_loop<1>(
			north_ghost.getStart(), north_ghost.getEnd(), [&](const std::array<int, 1> coord) {
				std::array<double, 2> real_coord;
				DomainTools::GetRealCoordGhost<2>(pinfos[1], {coord[0], ny}, real_coord);
				CHECK(north_ghost[coord] == Approx(f(real_coord)));
			});
		}
	}
	// patch 3
	{
		// check that center values weren't modified
		auto patch_3 = vec->getLocalData(0, d->getPatchInfoVector()[2].local_index);
		nested_loop<2>(patch_3.getStart(), patch_3.getEnd(), [&](const std::array<int, 2> coord) {
			std::array<double, 2> real_coord;
			DomainTools::GetRealCoord<2>(pinfos[2], coord, real_coord);
			CHECK(patch_3[coord] == f(real_coord));
		});
		// check that west values are not modified
		{
			auto west_ghost = patch_3.getGhostSliceOnSide(Side<2>::west(), 1);
			nested_loop<1>(west_ghost.getStart(), west_ghost.getEnd(),
			               [&](const std::array<int, 1> coord) { CHECK(west_ghost[coord] == 0); });
		}
		// check that east values are correct
		{
			auto east_ghost = patch_3.getGhostSliceOnSide(Side<2>::east(), 1);
			nested_loop<1>(
			east_ghost.getStart(), east_ghost.getEnd(), [&](const std::array<int, 1> coord) {
				std::array<double, 2> real_coord;
				DomainTools::GetRealCoordGhost<2>(pinfos[2], {nx, coord[0]}, real_coord);
				CHECK(east_ghost[coord] == Approx(f(real_coord)));
			});
		}
		// check that south values are not modified
		{
			auto south_ghost = patch_3.getGhostSliceOnSide(Side<2>::south(), 1);
			nested_loop<1>(
			south_ghost.getStart(), south_ghost.getEnd(), [&](const std::array<int, 1> coord) {
				std::array<double, 2> real_coord;
				DomainTools::GetRealCoordGhost<2>(pinfos[2], {coord[0], -1}, real_coord);
				CHECK(south_ghost[coord] == Approx(f(real_coord)));
			});
		}
		// check that north values are correct
		{
			auto north_ghost = patch_3.getGhostSliceOnSide(Side<2>::north(), 1);
			nested_loop<1>(north_ghost.getStart(), north_ghost.getEnd(),
			               [&](const std::array<int, 1> coord) { CHECK(north_ghost[coord] == 0); });
		}
	}
	// patch 4
	{
		// check that center values weren't modified
		auto patch_4 = vec->getLocalData(0, d->getPatchInfoVector()[3].local_index);
		nested_loop<2>(patch_4.getStart(), patch_4.getEnd(), [&](const std::array<int, 2> coord) {
			std::array<double, 2> real_coord;
			DomainTools::GetRealCoord<2>(pinfos[3], coord, real_coord);
			CHECK(patch_4[coord] == f(real_coord));
		});
		// check that west values correct
		{
			auto west_ghost = patch_4.getGhostSliceOnSide(Side<2>::west(), 1);
			nested_loop<1>(
			west_ghost.getStart(), west_ghost.getEnd(), [&](const std::array<int, 1> coord) {
				std::array<double, 2> real_coord;
				DomainTools::GetRealCoordGhost<2>(pinfos[3], {-1, coord[0]}, real_coord);
				CHECK(west_ghost[coord] == Approx(f(real_coord)));
			});
		}
		// check that east values are not modified
		{
			auto east_ghost = patch_4.getGhostSliceOnSide(Side<2>::east(), 1);
			nested_loop<1>(east_ghost.getStart(), east_ghost.getEnd(),
			               [&](const std::array<int, 1> coord) { CHECK(east_ghost[coord] == 0); });
		}
		// check that south values are correct
		{
			auto south_ghost = patch_4.getGhostSliceOnSide(Side<2>::south(), 1);
			nested_loop<1>(
			south_ghost.getStart(), south_ghost.getEnd(), [&](const std::array<int, 1> coord) {
				std::array<double, 2> real_coord;
				DomainTools::GetRealCoordGhost<2>(pinfos[3], {coord[0], -1}, real_coord);
				CHECK(south_ghost[coord] == Approx(f(real_coord)));
			});
		}
		// check that north values are not modified
		{
			auto north_ghost = patch_4.getGhostSliceOnSide(Side<2>::north(), 1);
			nested_loop<1>(north_ghost.getStart(), north_ghost.getEnd(),
			               [&](const std::array<int, 1> coord) { CHECK(north_ghost[coord] == 0); });
		}
	}
}
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller", "[BiLinearGhostFiller]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(2, 10);
	auto ny        = GENERATE(2, 10);
	int  num_ghost = 1;

	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d = domain_reader.getFinerDomain();

	shared_ptr<ValVector<2>> vec      = ValVector<2>::GetNewVector(d, 1);
	shared_ptr<ValVector<2>> expected = ValVector<2>::GetNewVector(d, 1);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 1 + ((x * 0.3) + y);
	};

	DomainTools::SetValues<2>(d, vec, f);
	DomainTools::SetValuesWithGhost<2>(d, expected, f);

	BiLinearGhostFiller blgf(d);
	blgf.fillGhost(vec);

	for (auto pinfo : d->getPatchInfoVector()) {
		INFO("Patch: " << pinfo.id);
		INFO("x:     " << pinfo.starts[0]);
		INFO("y:     " << pinfo.starts[1]);
		INFO("nx:    " << pinfo.ns[0]);
		INFO("ny:    " << pinfo.ns[1]);
		LocalData<2> vec_ld      = vec->getLocalData(0, pinfo.local_index);
		LocalData<2> expected_ld = expected->getLocalData(0, pinfo.local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			///
			REQUIRE(vec_ld[coord] == Approx(expected_ld[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			if (pinfo.hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo.getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 1> &coord) {
					               ///
					               INFO("coord:  " << coord[0]);
					               CHECK(vec_ghost[coord] == Approx(expected_ghost[coord]));
				               });
			}
		}
	}
}
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller ghost already set",
          "[BiLinearGhostFiller]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(2, 10);
	auto ny        = GENERATE(2, 10);
	int  num_ghost = 1;

	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d = domain_reader.getFinerDomain();

	shared_ptr<ValVector<2>> vec      = ValVector<2>::GetNewVector(d, 1);
	shared_ptr<ValVector<2>> expected = ValVector<2>::GetNewVector(d, 1);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[0];
		return 1 + ((x * 0.3) + y);
	};

	DomainTools::SetValuesWithGhost<2>(d, vec, f);
	DomainTools::SetValuesWithGhost<2>(d, expected, f);

	BiLinearGhostFiller blgf(d);
	blgf.fillGhost(vec);

	for (auto pinfo : d->getPatchInfoVector()) {
		INFO("Patch: " << pinfo.id);
		INFO("x:     " << pinfo.starts[0]);
		INFO("y:     " << pinfo.starts[1]);
		INFO("nx:    " << pinfo.ns[0]);
		INFO("ny:    " << pinfo.ns[1]);
		LocalData<2> vec_ld      = vec->getLocalData(0, pinfo.local_index);
		LocalData<2> expected_ld = expected->getLocalData(0, pinfo.local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			///
			REQUIRE(vec_ld[coord] == Approx(expected_ld[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			if (pinfo.hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo.getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 1> &coord) {
					               ///
					               INFO("coord:  " << coord[0]);
					               CHECK(vec_ghost[coord] == Approx(expected_ghost[coord]));
				               });
			}
		}
	}
}
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller two components", "[BiLinearGhostFiller]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(2, 10);
	auto ny        = GENERATE(2, 10);
	int  num_ghost = 1;

	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d = domain_reader.getFinerDomain();

	shared_ptr<ValVector<2>> vec      = ValVector<2>::GetNewVector(d, 2);
	shared_ptr<ValVector<2>> expected = ValVector<2>::GetNewVector(d, 2);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 1 + ((x * 0.3) + y);
	};
	auto g = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[0];
		return 99 + ((x * 7) + y * 0.1);
	};

	DomainTools::SetValues<2>(d, vec, f, g);
	DomainTools::SetValuesWithGhost<2>(d, expected, f, g);

	BiLinearGhostFiller blgf(d);
	blgf.fillGhost(vec);

	for (auto pinfo : d->getPatchInfoVector()) {
		INFO("Patch: " << pinfo.id);
		INFO("x:     " << pinfo.starts[0]);
		INFO("y:     " << pinfo.starts[1]);
		INFO("nx:    " << pinfo.ns[0]);
		INFO("ny:    " << pinfo.ns[1]);
		LocalData<2> vec_ld       = vec->getLocalData(0, pinfo.local_index);
		LocalData<2> expected_ld  = expected->getLocalData(0, pinfo.local_index);
		LocalData<2> vec_ld2      = vec->getLocalData(1, pinfo.local_index);
		LocalData<2> expected_ld2 = expected->getLocalData(1, pinfo.local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			///
			REQUIRE(vec_ld[coord] == Approx(expected_ld[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			if (pinfo.hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo.getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 1> &coord) {
					               ///
					               INFO("coord:  " << coord[0]);
					               CHECK(vec_ghost[coord] == Approx(expected_ghost[coord]));
				               });
			}
		}
		nested_loop<2>(vec_ld2.getStart(), vec_ld2.getEnd(), [&](const array<int, 2> &coord) {
			REQUIRE(vec_ld2[coord] == Approx(expected_ld2[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost      = vec_ld2.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost = expected_ld2.getGhostSliceOnSide(s, 1);
			if (pinfo.hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo.getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 1> &coord) {
					               INFO("coord:  " << coord[0]);
					               CHECK(vec_ghost[coord] == Approx(expected_ghost[coord]));
				               });
			}
		}
	}
}
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller ghost already set two components",
          "[BiLinearGhostFiller]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(2, 10);
	auto ny        = GENERATE(2, 10);
	int  num_ghost = 1;

	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d = domain_reader.getFinerDomain();

	shared_ptr<ValVector<2>> vec      = ValVector<2>::GetNewVector(d, 2);
	shared_ptr<ValVector<2>> expected = ValVector<2>::GetNewVector(d, 2);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[0];
		return 1 + ((x * 0.3) + y);
	};
	auto g = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[0];
		return 99 + ((x * 7) + y * 0.1);
	};

	DomainTools::SetValuesWithGhost<2>(d, vec, f, g);
	DomainTools::SetValuesWithGhost<2>(d, expected, f, g);

	BiLinearGhostFiller blgf(d);
	blgf.fillGhost(vec);

	for (auto pinfo : d->getPatchInfoVector()) {
		INFO("Patch: " << pinfo.id);
		INFO("x:     " << pinfo.starts[0]);
		INFO("y:     " << pinfo.starts[1]);
		INFO("nx:    " << pinfo.ns[0]);
		INFO("ny:    " << pinfo.ns[1]);
		LocalData<2> vec_ld       = vec->getLocalData(0, pinfo.local_index);
		LocalData<2> expected_ld  = expected->getLocalData(0, pinfo.local_index);
		LocalData<2> vec_ld2      = vec->getLocalData(1, pinfo.local_index);
		LocalData<2> expected_ld2 = expected->getLocalData(1, pinfo.local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			REQUIRE(vec_ld[coord] == Approx(expected_ld[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			if (pinfo.hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo.getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 1> &coord) {
					               INFO("coord:  " << coord[0]);
					               CHECK(vec_ghost[coord] == Approx(expected_ghost[coord]));
				               });
			}
		}
		nested_loop<2>(vec_ld2.getStart(), vec_ld2.getEnd(), [&](const array<int, 2> &coord) {
			REQUIRE(vec_ld2[coord] == Approx(expected_ld2[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost      = vec_ld2.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost = expected_ld2.getGhostSliceOnSide(s, 1);
			if (pinfo.hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo.getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 1> &coord) {
					               INFO("coord:  " << coord[0]);
					               CHECK(vec_ghost[coord] == Approx(expected_ghost[coord]));
				               });
			}
		}
	}
}