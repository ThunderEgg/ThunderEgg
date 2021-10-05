#include <ThunderEgg/PatchInfo.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("CoarseNbrInfo Serialization/Deserialization", "[CoarseNbrInfo]")
{
	CoarseNbrInfo<2> info;
	info.id             = 5;
	info.rank           = 1;
	info.orth_on_coarse = Orthant<2>::nw();
	// serialize and then deserialize
	char *buff = new char[info.serialize(nullptr)];
	info.serialize(buff);
	CoarseNbrInfo<2> out;
	out.deserialize(buff);
	delete[] buff;
	REQUIRE(out.id == 5);
	REQUIRE(out.rank == 1);
	REQUIRE(out.orth_on_coarse == Orthant<2>::nw());
}
TEST_CASE("CoarseNbrInfo to_json", "[CoarseNbrInfo]")
{
	CoarseNbrInfo<2> info;
	info.id             = GENERATE(1, 2, 3);
	info.rank           = GENERATE(0, 1, 2);
	info.orth_on_coarse = GENERATE(Orthant<2>::sw(), Orthant<3>::se(), Orthant<2>::nw());

	ThunderEgg::tpl::nlohmann::json j = info;

	CHECK(j["type"] == "COARSE");
	REQUIRE(j["ids"].is_array());
	CHECK(j["ids"].size() == 1);
	CHECK(j["ids"][0] == info.id);
	REQUIRE(j["ranks"].is_array());
	CHECK(j["ranks"].size() == 1);
	CHECK(j["ranks"][0] == info.rank);
	CHECK(j["orth_on_coarse"].get<Orthant<2>>() == info.orth_on_coarse);
}
TEST_CASE("CoarseNbrInfo from_json", "[CoarseNbrInfo]")
{
	int        id             = GENERATE(1, 2, 3);
	int        rank           = GENERATE(0, 1, 2);
	Orthant<2> orth_on_coarse = GENERATE(Orthant<2>::sw(), Orthant<3>::se(), Orthant<2>::nw());

	ThunderEgg::tpl::nlohmann::json j;
	j["type"]           = "COARSE";
	j["ids"]            = {id};
	j["ranks"]          = {rank};
	j["orth_on_coarse"] = orth_on_coarse;

	CoarseNbrInfo<2> info = j.get<CoarseNbrInfo<2>>();
	CHECK(info.id == id);
	CHECK(info.rank == rank);
	CHECK(info.orth_on_coarse == orth_on_coarse);
}