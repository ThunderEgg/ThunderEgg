#include "catch.hpp"
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/Timer.h>
#include <sstream>
using namespace std;
using namespace ThunderEgg;
static Domain<2> GetDomain()
{
	map<int, shared_ptr<PatchInfo<2>>> pinfo_map;

	int    n         = 10;
	double spacing   = 0.01;
	int    num_ghost = 1;

	pinfo_map[0].reset(new PatchInfo<2>());
	pinfo_map[0]->id = 0;
	pinfo_map[0]->ns.fill(n);
	pinfo_map[0]->spacings.fill(spacing);
	pinfo_map[0]->num_ghost_cells = num_ghost;
	Domain<2> d(pinfo_map, {n, n}, num_ghost);
	return d;
}
TEST_CASE("Two Timings Sequential Stop second before started", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.stop("A");
	REQUIRE_THROWS_AS(timer.stop("B"), RuntimeError);
}
TEST_CASE("Two Timings Sequential Stop with empty string second before started", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.stop("A");
	REQUIRE_THROWS_AS(timer.stop(""), RuntimeError);
}
TEST_CASE("Two Timings Nested Wrong Order", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.start("B");
	REQUIRE_THROWS_AS(timer.stop("A"), RuntimeError);
}
TEST_CASE("Timer ostream operator throws with unfinished timing", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.start("B");
	stringstream writer;
	REQUIRE_THROWS_AS(writer << timer, RuntimeError);
}
TEST_CASE("Timer Two DomainTimings Sequential Stop second before started", "[Timer]")
{
	Timer timer;
	timer.addDomain(0, GetDomain());
	timer.startDomainTiming(0, "A");
	timer.stopDomainTiming(0, "A");
	int          id   = GENERATE(0, 1);
	const string name = GENERATE("A", "B", "");
	REQUIRE_THROWS_AS(timer.stopDomainTiming(id, name), RuntimeError);
}
TEST_CASE("Timer DomainTimings Nested Wrong id on stop", "[Timer]")
{
	Timer timer;
	timer.addDomain(0, GetDomain());
	timer.startDomainTiming(0, "A");
	REQUIRE_THROWS_AS(timer.stopDomainTiming(1, "A"), RuntimeError);
}
TEST_CASE("Timer DomainTimings Nested Wrong name on stop", "[Timer]")
{
	Timer timer;
	timer.addDomain(0, GetDomain());
	timer.startDomainTiming(0, "A");
	REQUIRE_THROWS_AS(timer.stopDomainTiming(0, "blah"), RuntimeError);
}
TEST_CASE("Timer DomainTimings Nested Wrong name and id on stop", "[Timer]")
{
	Timer timer;
	timer.addDomain(0, GetDomain());
	timer.startDomainTiming(0, "A");
	REQUIRE_THROWS_AS(timer.stopDomainTiming(1, "blah"), RuntimeError);
}
TEST_CASE("Timer Two DomainTimings Nested Wrong Order", "[Timer]")
{
	Timer timer;
	timer.addDomain(0, GetDomain());
	timer.addDomain(1, GetDomain());
	timer.startDomainTiming(0, "A");
	timer.startDomainTiming(1, "B");
	REQUIRE_THROWS_AS(timer.stopDomainTiming(0, "A"), RuntimeError);
}
TEST_CASE("Timer addDomain twice fails", "[Timer]")
{
	Timer timer;
	timer.addDomain(0, GetDomain());
	REQUIRE_THROWS_AS(timer.addDomain(0, GetDomain()), RuntimeError);
}
TEST_CASE("Timer startDomainTiming fails without added domain", "[Timer]")
{
	Timer timer;
	REQUIRE_THROWS_AS(timer.startDomainTiming(0, "A"), RuntimeError);
}
TEST_CASE("Timer to_json empty timer", "[Timer]")
{
	Timer          timer;
	nlohmann::json j = timer;
	INFO(j.dump(4));
	REQUIRE(j == nullptr);
}
TEST_CASE("Timer to_json unassociated timing", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.stop("A");
	const nlohmann::json j = timer;
	INFO(j.dump(4));
	REQUIRE(j != nullptr);
	CHECK(j.size() == 1);
	CHECK(j["timings"].is_array());
	CHECK(j["timings"].size() == 1);
	CHECK(j["timings"][0]["min"].is_number());
	CHECK(j["timings"][0]["max"].is_number());
	CHECK(j["timings"][0]["sum"].is_number());
	CHECK(j["timings"][0]["num_calls"].is_number());
	CHECK(j["timings"][0]["name"] == "A");
}
TEST_CASE("Timer to_json two unassociated timings sequential", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.stop("A");
	timer.start("B");
	timer.stop("B");
	const nlohmann::json j = timer;
	INFO(j.dump(4));
	REQUIRE(j != nullptr);
	CHECK(j.size() == 1);
	CHECK(j["timings"].is_array());
	CHECK(j["timings"].size() == 2);
	CHECK(j["timings"][0]["min"].is_number());
	CHECK(j["timings"][0]["max"].is_number());
	CHECK(j["timings"][0]["sum"].is_number());
	CHECK(j["timings"][0]["num_calls"].is_number());
	CHECK(j["timings"][0]["name"] == "A");
	CHECK(j["timings"][1]["min"].is_number());
	CHECK(j["timings"][1]["max"].is_number());
	CHECK(j["timings"][1]["sum"].is_number());
	CHECK(j["timings"][1]["num_calls"].is_number());
	CHECK(j["timings"][1]["name"] == "B");
}
TEST_CASE("Timer to_json nested timing", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.start("B");
	timer.stop("B");
	timer.stop("A");
	const nlohmann::json j = timer;
	INFO(j.dump(4));
	REQUIRE(j != nullptr);
	CHECK(j.size() == 1);
	CHECK(j["timings"].is_array());
	CHECK(j["timings"].size() == 1);
	CHECK(j["timings"][0]["min"].is_number());
	CHECK(j["timings"][0]["max"].is_number());
	CHECK(j["timings"][0]["sum"].is_number());
	CHECK(j["timings"][0]["num_calls"].is_number());
	CHECK(j["timings"][0]["name"] == "A");
	CHECK(j["timings"][0]["timings"].is_array());
	CHECK(j["timings"][0]["timings"].size() == 1);
}
TEST_CASE("Timer to_json domain timing", "[Timer]")
{
	Timer timer;
	timer.addDomain(0, GetDomain());
	timer.startDomainTiming(0, "A");
	timer.stopDomainTiming(0, "A");
	const nlohmann::json j = timer;
	REQUIRE(j != nullptr);
	INFO(j.dump(4));
	CHECK(j.size() == 2);
	CHECK(j["domains"].is_array());
	CHECK(j["domains"].size() == 1);
	CHECK(j["domains"][0].is_array());
	CHECK(j["domains"][0].size() == 2);
	CHECK(j["timings"].is_array());
	CHECK(j["timings"].size() == 1);
	CHECK(j["timings"][0]["min"].is_number());
	CHECK(j["timings"][0]["max"].is_number());
	CHECK(j["timings"][0]["sum"].is_number());
	CHECK(j["timings"][0]["num_calls"].is_number());
	CHECK(j["timings"][0]["name"] == "A");
	CHECK(j["timings"][0]["domain_id"] == 0);
}