#include "utils/DomainReader.h"
#include <ThunderEgg/BiQuadraticGhostFiller.h>
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/ValVectorGenerator.h>

#include "catch.hpp"
using namespace std;
using namespace ThunderEgg;

constexpr auto single_mesh_file  = "mesh_inputs/2d_uniform_2x2_mpi1.json";
constexpr auto refined_mesh_file = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json";
constexpr auto cross_mesh_file   = "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json";
TEST_CASE("exchange uniform 2D quad BiQuadraticGhostFiller", "[BiQuadraticGhostFiller]")
{
	auto   nx        = GENERATE(10, 13);
	auto   ny        = GENERATE(10, 13);
	double startx    = 0;
	double starty    = 0;
	double lengthx   = 1;
	double lengthy   = 1;
	int    num_ghost = 1;

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 1 + ((x * x + 0.3 * x + 0.3) + (0.5 * y * y));
	};

	/*
	 *  3 | 4
	 *  --+---
	 *  1 | 2
	 */
	map<int, shared_ptr<PatchInfo<2>>> pinfo_map;

	pinfo_map[1].reset(new PatchInfo<2>());
	pinfo_map[1]->id              = 1;
	pinfo_map[1]->ns              = {nx, ny};
	pinfo_map[1]->spacings        = {lengthx / nx, lengthy / ny};
	pinfo_map[1]->starts          = {startx, starty};
	pinfo_map[1]->num_ghost_cells = num_ghost;
	pinfo_map[1]->nbr_info[Side<2>::east().getIndex()].reset(new NormalNbrInfo<2>(2));
	pinfo_map[1]->nbr_info[Side<2>::north().getIndex()].reset(new NormalNbrInfo<2>(3));

	pinfo_map[2].reset(new PatchInfo<2>());
	pinfo_map[2]->id              = 2;
	pinfo_map[2]->ns              = {nx, ny};
	pinfo_map[2]->spacings        = {lengthx / nx, lengthy / ny};
	pinfo_map[2]->starts          = {startx + lengthx, starty};
	pinfo_map[2]->num_ghost_cells = num_ghost;
	pinfo_map[2]->nbr_info[Side<2>::west().getIndex()].reset(new NormalNbrInfo<2>(1));
	pinfo_map[2]->nbr_info[Side<2>::north().getIndex()].reset(new NormalNbrInfo<2>(4));

	pinfo_map[3].reset(new PatchInfo<2>());
	pinfo_map[3]->id              = 3;
	pinfo_map[3]->ns              = {nx, ny};
	pinfo_map[3]->spacings        = {lengthx / nx, lengthy / ny};
	pinfo_map[3]->starts          = {startx, starty + lengthy};
	pinfo_map[3]->num_ghost_cells = num_ghost;
	pinfo_map[3]->nbr_info[Side<2>::east().getIndex()].reset(new NormalNbrInfo<2>(4));
	pinfo_map[3]->nbr_info[Side<2>::south().getIndex()].reset(new NormalNbrInfo<2>(1));

	pinfo_map[4].reset(new PatchInfo<2>());
	pinfo_map[4]->id              = 4;
	pinfo_map[4]->ns              = {nx, ny};
	pinfo_map[4]->spacings        = {lengthx / nx, lengthy / ny};
	pinfo_map[4]->starts          = {startx + lengthx, starty + lengthy};
	pinfo_map[4]->num_ghost_cells = num_ghost;
	pinfo_map[4]->nbr_info[Side<2>::west().getIndex()].reset(new NormalNbrInfo<2>(3));
	pinfo_map[4]->nbr_info[Side<2>::south().getIndex()].reset(new NormalNbrInfo<2>(2));

	shared_ptr<Domain<2>> d(new Domain<2>(pinfo_map, {nx, ny}, num_ghost));

	shared_ptr<ValVector<2>> vec(new ValVector<2>(MPI_COMM_WORLD, pinfo_map[1]->ns, num_ghost, 4));

	DomainTools<2>::setValues(d, vec, f);

	BiQuadraticGhostFiller blgf(d);
	blgf.fillGhost(vec);

	// patch 1
	{
		// check that center values weren't modified
		auto patch_1 = vec->getLocalData(d->getPatchInfoMap()[1]->local_index);
		nested_loop<2>(patch_1.getStart(), patch_1.getEnd(), [&](const std::array<int, 2> coord) {
			std::array<double, 2> real_coord;
			DomainTools<2>::getRealCoord(pinfo_map[1], coord, real_coord);
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
				DomainTools<2>::getRealCoordGhost(pinfo_map[1], {nx, coord[0]}, real_coord);
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
				DomainTools<2>::getRealCoordGhost(pinfo_map[1], {coord[0], ny}, real_coord);
				CHECK(north_ghost[coord] == Approx(f(real_coord)));
			});
		}
	}
	// patch 2
	{
		// check that center values weren't modified
		auto patch_2 = vec->getLocalData(d->getPatchInfoMap()[2]->local_index);
		nested_loop<2>(patch_2.getStart(), patch_2.getEnd(), [&](const std::array<int, 2> coord) {
			std::array<double, 2> real_coord;
			DomainTools<2>::getRealCoord(pinfo_map[2], coord, real_coord);
			CHECK(patch_2[coord] == f(real_coord));
		});
		// check that west values correct
		{
			auto west_ghost = patch_2.getGhostSliceOnSide(Side<2>::west(), 1);
			nested_loop<1>(
			west_ghost.getStart(), west_ghost.getEnd(), [&](const std::array<int, 1> coord) {
				std::array<double, 2> real_coord;
				DomainTools<2>::getRealCoordGhost(pinfo_map[2], {-1, coord[0]}, real_coord);
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
				DomainTools<2>::getRealCoordGhost(pinfo_map[2], {coord[0], ny}, real_coord);
				CHECK(north_ghost[coord] == Approx(f(real_coord)));
			});
		}
	}
	// patch 3
	{
		// check that center values weren't modified
		auto patch_1 = vec->getLocalData(d->getPatchInfoMap()[3]->local_index);
		nested_loop<2>(patch_1.getStart(), patch_1.getEnd(), [&](const std::array<int, 2> coord) {
			std::array<double, 2> real_coord;
			DomainTools<2>::getRealCoord(pinfo_map[3], coord, real_coord);
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
				DomainTools<2>::getRealCoordGhost(pinfo_map[3], {nx, coord[0]}, real_coord);
				CHECK(east_ghost[coord] == Approx(f(real_coord)));
			});
		}
		// check that south values are not modified
		{
			auto south_ghost = patch_1.getGhostSliceOnSide(Side<2>::south(), 1);
			nested_loop<1>(
			south_ghost.getStart(), south_ghost.getEnd(), [&](const std::array<int, 1> coord) {
				std::array<double, 2> real_coord;
				DomainTools<2>::getRealCoordGhost(pinfo_map[3], {coord[0], -1}, real_coord);
				CHECK(south_ghost[coord] == Approx(f(real_coord)));
			});
		}
		// check that north values are correct
		{
			auto north_ghost = patch_1.getGhostSliceOnSide(Side<2>::north(), 1);
			nested_loop<1>(north_ghost.getStart(), north_ghost.getEnd(),
			               [&](const std::array<int, 1> coord) { CHECK(north_ghost[coord] == 0); });
		}
	}
	// patch 4
	{
		// check that center values weren't modified
		auto patch_2 = vec->getLocalData(d->getPatchInfoMap()[4]->local_index);
		nested_loop<2>(patch_2.getStart(), patch_2.getEnd(), [&](const std::array<int, 2> coord) {
			std::array<double, 2> real_coord;
			DomainTools<2>::getRealCoord(pinfo_map[4], coord, real_coord);
			CHECK(patch_2[coord] == f(real_coord));
		});
		// check that west values correct
		{
			auto west_ghost = patch_2.getGhostSliceOnSide(Side<2>::west(), 1);
			nested_loop<1>(
			west_ghost.getStart(), west_ghost.getEnd(), [&](const std::array<int, 1> coord) {
				std::array<double, 2> real_coord;
				DomainTools<2>::getRealCoordGhost(pinfo_map[4], {-1, coord[0]}, real_coord);
				CHECK(west_ghost[coord] == Approx(f(real_coord)));
			});
		}
		// check that east values are not modified
		{
			auto east_ghost = patch_2.getGhostSliceOnSide(Side<2>::east(), 1);
			nested_loop<1>(east_ghost.getStart(), east_ghost.getEnd(),
			               [&](const std::array<int, 1> coord) { CHECK(east_ghost[coord] == 0); });
		}
		// check that south values are correct
		{
			auto south_ghost = patch_2.getGhostSliceOnSide(Side<2>::south(), 1);
			nested_loop<1>(
			south_ghost.getStart(), south_ghost.getEnd(), [&](const std::array<int, 1> coord) {
				std::array<double, 2> real_coord;
				DomainTools<2>::getRealCoordGhost(pinfo_map[4], {coord[0], -1}, real_coord);
				CHECK(south_ghost[coord] == Approx(f(real_coord)));
			});
		}
		// check that north values are not modified
		{
			auto north_ghost = patch_2.getGhostSliceOnSide(Side<2>::north(), 1);
			nested_loop<1>(north_ghost.getStart(), north_ghost.getEnd(),
			               [&](const std::array<int, 1> coord) { CHECK(north_ghost[coord] == 0); });
		}
	}
}
TEST_CASE("exchange various meshes 2D BiQuadraticGhostFiller", "[BiQuadraticGhostFiller]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(10, 13);
	auto ny        = GENERATE(10, 13);
	int  num_ghost = 1;

	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d = domain_reader.getFinerDomain();

	shared_ptr<ValVector<2>> vec      = ValVector<2>::GetNewVector(d);
	shared_ptr<ValVector<2>> expected = ValVector<2>::GetNewVector(d);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 1 + ((x * x + 0.3 * x + 0.3) + (0.5 * y * y));
	};

	DomainTools<2>::setValues(d, vec, f);
	DomainTools<2>::setValuesWithGhost(d, expected, f);

	BiQuadraticGhostFiller blgf(d);
	blgf.fillGhost(vec);

	for (auto pinfo : d->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		LocalData<2> vec_ld      = vec->getLocalData(pinfo->local_index);
		LocalData<2> expected_ld = expected->getLocalData(pinfo->local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			///
			REQUIRE(vec_ld[coord] == Approx(expected_ld[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			if (pinfo->hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo->getNbrType(s));
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
TEST_CASE("exchange various meshes 2D BiQuadraticGhostFiller ghost already set",
          "[BiQuadraticGhostFiller]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(10, 13);
	auto ny        = GENERATE(10, 13);
	int  num_ghost = 1;

	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d = domain_reader.getFinerDomain();

	shared_ptr<ValVector<2>> vec      = ValVector<2>::GetNewVector(d);
	shared_ptr<ValVector<2>> expected = ValVector<2>::GetNewVector(d);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[0];
		return 1 + ((x * x + 0.3 * x + 0.3) + (0.5 * y * y));
	};

	DomainTools<2>::setValuesWithGhost(d, vec, f);
	DomainTools<2>::setValuesWithGhost(d, expected, f);

	BiQuadraticGhostFiller blgf(d);
	blgf.fillGhost(vec);

	for (auto pinfo : d->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		LocalData<2> vec_ld      = vec->getLocalData(pinfo->local_index);
		LocalData<2> expected_ld = expected->getLocalData(pinfo->local_index);
		nested_loop<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			///
			REQUIRE(vec_ld[coord] == Approx(expected_ld[coord]));
		});
		for (Side<2> s : Side<2>::getValues()) {
			LocalData<1> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<1> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			if (pinfo->hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo->getNbrType(s));
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