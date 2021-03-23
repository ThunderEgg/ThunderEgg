#include "utils/DomainReader.h"
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/TriLinearGhostFiller.h>
#include <ThunderEgg/ValVectorGenerator.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

constexpr auto single_mesh_file  = "mesh_inputs/3d_uniform_2x2x2_mpi1.json";
constexpr auto refined_mesh_file = "mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json";

TEST_CASE("exchange various meshes 3D TriLinearGhostFiller", "[TriLinearGhostFiller]")
{
	auto mesh_file = GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(10, 2);
	auto ny        = GENERATE(10, 2);
	auto nz        = GENERATE(10, 2);
	int  num_ghost = 1;

	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d = domain_reader.getFinerDomain();

	shared_ptr<ValVector<3>> vec      = ValVector<3>::GetNewVector(d, 1);
	shared_ptr<ValVector<3>> expected = ValVector<3>::GetNewVector(d, 1);

	auto f = [&](const std::array<double, 3> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		double z = coord[2];
		return 1 + 0.5 * x + y + 7 * z;
	};

	DomainTools::SetValues<3>(d, vec, f);
	DomainTools::SetValuesWithGhost<3>(d, expected, f);

	TriLinearGhostFiller tlgf(d);
	tlgf.fillGhost(vec);

	for (auto pinfo : d->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("z:     " << pinfo->starts[2]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		INFO("nz:    " << pinfo->ns[2]);
		LocalData<3> vec_ld      = vec->getLocalData(0, pinfo->local_index);
		LocalData<3> expected_ld = expected->getLocalData(0, pinfo->local_index);
		nested_loop<3>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 3> &coord) {
			REQUIRE(vec_ld[coord] == Catch::Approx(expected_ld[coord]));
		});
		for (Side<3> s : Side<3>::getValues()) {
			LocalData<2> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<2> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			if (pinfo->hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo->getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 2> &coord) {
					               INFO("coord:  " << coord[0] << ", " << coord[1]);
					               CHECK(vec_ghost[coord] == Catch::Approx(expected_ghost[coord]));
				               });
			}
		}
	}
}
TEST_CASE("exchange various meshes 3D TriLinearGhostFiller two components",
          "[TriLinearGhostFiller]")
{
	auto mesh_file = GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(10, 2);
	auto ny        = GENERATE(10, 2);
	auto nz        = GENERATE(10, 2);
	int  num_ghost = 1;

	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d = domain_reader.getFinerDomain();

	shared_ptr<ValVector<3>> vec      = ValVector<3>::GetNewVector(d, 2);
	shared_ptr<ValVector<3>> expected = ValVector<3>::GetNewVector(d, 2);

	auto f = [&](const std::array<double, 3> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		double z = coord[2];
		return 1 + 0.5 * x + y + 7 * z;
	};
	auto g = [&](const std::array<double, 3> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		double z = coord[2];
		return 100 + 20 * x + 9 * y + 7 * z;
	};

	DomainTools::SetValues<3>(d, vec, f, g);
	DomainTools::SetValuesWithGhost<3>(d, expected, f, g);

	TriLinearGhostFiller tlgf(d);
	tlgf.fillGhost(vec);

	for (auto pinfo : d->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("z:     " << pinfo->starts[2]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		INFO("nz:    " << pinfo->ns[2]);
		LocalData<3> vec_ld       = vec->getLocalData(0, pinfo->local_index);
		LocalData<3> expected_ld  = expected->getLocalData(0, pinfo->local_index);
		LocalData<3> vec_ld2      = vec->getLocalData(1, pinfo->local_index);
		LocalData<3> expected_ld2 = expected->getLocalData(1, pinfo->local_index);
		nested_loop<3>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 3> &coord) {
			REQUIRE(vec_ld[coord] == Catch::Approx(expected_ld[coord]));
		});
		for (Side<3> s : Side<3>::getValues()) {
			LocalData<2> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<2> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			if (pinfo->hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo->getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 2> &coord) {
					               INFO("coord:  " << coord[0] << ", " << coord[1]);
					               CHECK(vec_ghost[coord] == Catch::Approx(expected_ghost[coord]));
				               });
			}
		}
		nested_loop<3>(vec_ld2.getStart(), vec_ld2.getEnd(), [&](const array<int, 3> &coord) {
			REQUIRE(vec_ld2[coord] == Catch::Approx(expected_ld2[coord]));
		});
		for (Side<3> s : Side<3>::getValues()) {
			LocalData<2> vec_ghost      = vec_ld2.getGhostSliceOnSide(s, 1);
			LocalData<2> expected_ghost = expected_ld2.getGhostSliceOnSide(s, 1);
			if (pinfo->hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo->getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 2> &coord) {
					               INFO("coord:  " << coord[0] << ", " << coord[1]);
					               CHECK(vec_ghost[coord] == Catch::Approx(expected_ghost[coord]));
				               });
			}
		}
	}
}
TEST_CASE("exchange various meshes 3D TriLinearGhostFiller ghost already set two components",
          "[TriLinearGhostFiller]")
{
	auto mesh_file = GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(10, 2);
	auto ny        = GENERATE(10, 2);
	auto nz        = GENERATE(10, 2);
	int  num_ghost = 1;

	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d = domain_reader.getFinerDomain();

	shared_ptr<ValVector<3>> vec      = ValVector<3>::GetNewVector(d, 2);
	shared_ptr<ValVector<3>> expected = ValVector<3>::GetNewVector(d, 2);

	auto f = [&](const std::array<double, 3> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		double z = coord[2];
		return 1 + 0.5 * x + y + 7 * z;
	};
	auto g = [&](const std::array<double, 3> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		double z = coord[2];
		return 100 + 20 * x + 9 * y + 7 * z;
	};

	DomainTools::SetValuesWithGhost<3>(d, vec, f, g);
	DomainTools::SetValuesWithGhost<3>(d, expected, f, g);

	TriLinearGhostFiller tlgf(d);
	tlgf.fillGhost(vec);

	for (auto pinfo : d->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("z:     " << pinfo->starts[2]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		INFO("nz:    " << pinfo->ns[2]);
		LocalData<3> vec_ld       = vec->getLocalData(0, pinfo->local_index);
		LocalData<3> expected_ld  = expected->getLocalData(0, pinfo->local_index);
		LocalData<3> vec_ld2      = vec->getLocalData(1, pinfo->local_index);
		LocalData<3> expected_ld2 = expected->getLocalData(1, pinfo->local_index);
		nested_loop<3>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 3> &coord) {
			REQUIRE(vec_ld[coord] == Catch::Approx(expected_ld[coord]));
		});
		for (Side<3> s : Side<3>::getValues()) {
			LocalData<2> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<2> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			if (pinfo->hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo->getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 2> &coord) {
					               INFO("coord:  " << coord[0] << ", " << coord[1]);
					               CHECK(vec_ghost[coord] == Catch::Approx(expected_ghost[coord]));
				               });
			}
		}
		nested_loop<3>(vec_ld2.getStart(), vec_ld2.getEnd(), [&](const array<int, 3> &coord) {
			REQUIRE(vec_ld2[coord] == Catch::Approx(expected_ld2[coord]));
		});
		for (Side<3> s : Side<3>::getValues()) {
			LocalData<2> vec_ghost      = vec_ld2.getGhostSliceOnSide(s, 1);
			LocalData<2> expected_ghost = expected_ld2.getGhostSliceOnSide(s, 1);
			if (pinfo->hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo->getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 2> &coord) {
					               INFO("coord:  " << coord[0] << ", " << coord[1]);
					               CHECK(vec_ghost[coord] == Catch::Approx(expected_ghost[coord]));
				               });
			}
		}
	}
}
TEST_CASE("exchange various meshes 3D TriLinearGhostFiller ghost already set",
          "[TriLinearGhostFiller]")
{
	auto mesh_file = GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(10, 2);
	auto ny        = GENERATE(10, 2);
	auto nz        = GENERATE(10, 2);
	int  num_ghost = 1;

	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d = domain_reader.getFinerDomain();

	shared_ptr<ValVector<3>> vec      = ValVector<3>::GetNewVector(d, 1);
	shared_ptr<ValVector<3>> expected = ValVector<3>::GetNewVector(d, 1);

	auto f = [&](const std::array<double, 3> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		double z = coord[2];
		return 1 + 0.5 * x + y + 7 * z;
	};

	DomainTools::SetValuesWithGhost<3>(d, vec, f);
	DomainTools::SetValuesWithGhost<3>(d, expected, f);

	TriLinearGhostFiller tlgf(d);
	tlgf.fillGhost(vec);

	for (auto pinfo : d->getPatchInfoVector()) {
		INFO("Patch: " << pinfo->id);
		INFO("x:     " << pinfo->starts[0]);
		INFO("y:     " << pinfo->starts[1]);
		INFO("z:     " << pinfo->starts[2]);
		INFO("nx:    " << pinfo->ns[0]);
		INFO("ny:    " << pinfo->ns[1]);
		INFO("nz:    " << pinfo->ns[2]);
		LocalData<3> vec_ld      = vec->getLocalData(0, pinfo->local_index);
		LocalData<3> expected_ld = expected->getLocalData(0, pinfo->local_index);
		nested_loop<3>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 3> &coord) {
			REQUIRE(vec_ld[coord] == Catch::Approx(expected_ld[coord]));
		});
		for (Side<3> s : Side<3>::getValues()) {
			LocalData<2> vec_ghost      = vec_ld.getGhostSliceOnSide(s, 1);
			LocalData<2> expected_ghost = expected_ld.getGhostSliceOnSide(s, 1);
			if (pinfo->hasNbr(s)) {
				INFO("side:      " << s);
				INFO("nbr-type:  " << pinfo->getNbrType(s));
				nested_loop<1>(vec_ghost.getStart(), vec_ghost.getEnd(),
				               [&](const array<int, 2> &coord) {
					               INFO("coord:  " << coord[0] << ", " << coord[1]);
					               CHECK(vec_ghost[coord] == Catch::Approx(expected_ghost[coord]));
				               });
			}
		}
	}
}
TEST_CASE("TriLinearGhostFiller constructor throws error with odd number of cells",
          "[TriLinearGhostFiller]")
{
	auto mesh_file = GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file);
	INFO("MESH: " << mesh_file);
	auto axis = GENERATE(0, 1, 2);
	INFO("axis: " << axis);
	int n_even    = 10;
	int n_odd     = 11;
	int num_ghost = 1;

	array<int, 3> ns;
	ns.fill(n_even);
	ns[axis] = n_odd;
	DomainReader<3>       domain_reader(mesh_file, ns, num_ghost);
	shared_ptr<Domain<3>> d = domain_reader.getFinerDomain();

	CHECK_THROWS_AS(TriLinearGhostFiller(d), RuntimeError);
}