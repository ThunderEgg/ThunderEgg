#include <ThunderEgg/Domain.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Domain constructors work", "[Domain]")
{
	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	auto n         = GENERATE(1, 2, 10, 13);
	auto spacing   = GENERATE(0.01, 1.0, 3.14);
	auto num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id = 0;
	pinfos[0]->ns.fill(n);
	pinfos[0]->spacings.fill(spacing);
	pinfos[0]->num_ghost_cells = num_ghost;
	Domain<2> d(0, {n, n}, num_ghost, pinfos.begin(), pinfos.end());

	// check getters
	for (int ni : d.getNs()) {
		CHECK(ni == n);
	}
	CHECK(d.getNumGlobalPatches() == 1);
	CHECK(d.getNumLocalPatches() == 1);
	CHECK(d.getNumGlobalCells() == n * n);
	CHECK(d.getNumLocalCells() == n * n);
	CHECK(d.getNumLocalCellsWithGhost() == (n + 2 * num_ghost) * (n + 2 * num_ghost));
	CHECK(d.getNumCellsInPatch() == n * n);
	CHECK(d.getNumGhostCells() == num_ghost);
	CHECK(d.volume() == Catch::Approx(spacing * spacing * n * n));
	// TODO Check intigrate
}
TEST_CASE("Domain setTimer", "[Domain]")
{
	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    n         = 10;
	double spacing   = 0.01;
	int    num_ghost = 1;

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id = 0;
	pinfos[0]->ns.fill(n);
	pinfos[0]->spacings.fill(spacing);
	pinfos[0]->num_ghost_cells = num_ghost;
	Domain<2> d(0, {n, n}, num_ghost, pinfos.begin(), pinfos.end());

	auto timer = make_shared<Timer>(MPI_COMM_WORLD);
	d.setTimer(timer);
	CHECK(d.getTimer() == timer);
	CHECK(d.hasTimer());
}
TEST_CASE("Domain setTimer adds domain to timer", "[Domain]")
{
	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    n         = 10;
	double spacing   = 0.01;
	int    num_ghost = 1;

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id = 0;
	pinfos[0]->ns.fill(n);
	pinfos[0]->spacings.fill(spacing);
	pinfos[0]->num_ghost_cells = num_ghost;
	Domain<2> d(0, {n, n}, num_ghost, pinfos.begin(), pinfos.end());

	auto timer = make_shared<Timer>(MPI_COMM_WORLD);
	d.setTimer(timer);
	// will throw exception if domain not added to timer
	timer->startDomainTiming(0, "A");
	timer->stopDomainTiming(0, "A");
}
TEST_CASE("Domain getTimer default is no timer", "[Domain]")
{
	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    n         = 10;
	double spacing   = 0.01;
	int    num_ghost = 1;

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id = 0;
	pinfos[0]->ns.fill(n);
	pinfos[0]->spacings.fill(spacing);
	pinfos[0]->num_ghost_cells = num_ghost;
	Domain<2> d(0, {n, n}, num_ghost, pinfos.begin(), pinfos.end());

	CHECK(d.getTimer() == nullptr);
	CHECK_FALSE(d.hasTimer());
}
TEST_CASE("Domain id in constructor", "[Domain]")
{
	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    n         = 10;
	double spacing   = 0.01;
	int    num_ghost = 1;

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id = 0;
	pinfos[0]->ns.fill(n);
	pinfos[0]->spacings.fill(spacing);
	pinfos[0]->num_ghost_cells = num_ghost;
	auto      id               = GENERATE(1, 2, 9);
	Domain<2> d(id, {n, n}, num_ghost, pinfos.begin(), pinfos.end());

	CHECK(d.getId() == id);
}
TEST_CASE("Domain to_json", "[Domain]")
{
	vector<shared_ptr<PatchInfo<2>>> pinfos(1);

	int    n         = 10;
	double spacing   = 0.01;
	int    num_ghost = 1;

	pinfos[0].reset(new PatchInfo<2>());
	pinfos[0]->id = 0;
	pinfos[0]->ns.fill(n);
	pinfos[0]->spacings.fill(spacing);
	pinfos[0]->num_ghost_cells = num_ghost;
	Domain<2> d(0, {n, n}, num_ghost, pinfos.begin(), pinfos.end());

	nlohmann::json j = d;
	REQUIRE(j.is_array());
	REQUIRE(j.size() == 1);
	REQUIRE(j[0]["id"] == 0);
}
