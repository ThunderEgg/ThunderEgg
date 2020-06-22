#include "utils/DomainReader.h"
#include <Thunderegg/BiLinearGhostFiller.h>
#include <Thunderegg/DomainTools.h>
#include <Thunderegg/ValVector.h>

#include "catch.hpp"
using namespace std;
using namespace Thunderegg;

constexpr auto uniform = "mesh_inputs/2d_uniform_2x2_nw_on_1_mpi2.json";
constexpr auto refined = "mesh_inputs/2d_uniform_2x2_refined_nw_on_1_mpi2.json";
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller", "[BiLinearGhostFiller]")
{
	auto mesh_file = GENERATE(as<std::string>{}, uniform, refined);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(2, 10);
	auto ny        = GENERATE(2, 10);
	int  num_ghost = 1;

	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d = domain_reader.getFinerDomain();

	shared_ptr<ValVector<2>> vec      = ValVector<2>::GetNewVector(d);
	shared_ptr<ValVector<2>> expected = ValVector<2>::GetNewVector(d);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[1];
		return 1 + ((x * 0.3) + y);
	};

	DomainTools<2>::setValues(d, vec, f);
	DomainTools<2>::setValuesWithGhost(d, expected, f);

	BiLinearGhostFiller blgf(d);
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
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller ghost already set",
          "[BiLinearGhostFiller]")
{
	auto mesh_file = GENERATE(as<std::string>{}, uniform, refined);
	INFO("MESH: " << mesh_file);
	auto nx        = GENERATE(2, 10);
	auto ny        = GENERATE(2, 10);
	int  num_ghost = 1;

	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d = domain_reader.getFinerDomain();

	shared_ptr<ValVector<2>> vec      = ValVector<2>::GetNewVector(d);
	shared_ptr<ValVector<2>> expected = ValVector<2>::GetNewVector(d);

	auto f = [&](const std::array<double, 2> coord) -> double {
		double x = coord[0];
		double y = coord[0];
		return 1 + ((x * 0.3) + y);
	};

	DomainTools<2>::setValuesWithGhost(d, vec, f);
	DomainTools<2>::setValuesWithGhost(d, expected, f);

	BiLinearGhostFiller blgf(d);
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