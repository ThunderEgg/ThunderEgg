#include <ThunderEgg/PatchInfo.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("NormalNbrInfo getNbrType works", "[NormalNbrInfo]")
{
	NbrInfo<3> *info = new NormalNbrInfo<3>();
	REQUIRE(info->getNbrType() == NbrType::Normal);
	delete info;
}

TEST_CASE("NormalNbrInfo Serialization/Deserialization", "[NormalNbrInfo]")
{
	NormalNbrInfo<3> info;
	info.id   = 5;
	info.rank = 1;
	// serialize and then deserialize
	char *buff = new char[info.serialize(nullptr)];
	info.serialize(buff);
	NormalNbrInfo<3> out;
	out.deserialize(buff);
	delete[] buff;
	REQUIRE(out.id == 5);
	REQUIRE(out.rank == 1);
}
TEST_CASE("NormalNbrInfo to_json", "[NormalNbrInfo]")
{
	NormalNbrInfo<3> info;
	info.id   = GENERATE(1, 2, 3);
	info.rank = GENERATE(0, 1, 2);

	nlohmann::json j = info;

	CHECK(j["type"] == "NORMAL");
	REQUIRE(j["ids"].is_array());
	CHECK(j["ids"].size() == 1);
	CHECK(j["ids"][0] == info.id);
	REQUIRE(j["ranks"].is_array());
	CHECK(j["ranks"].size() == 1);
	CHECK(j["ranks"][0] == info.rank);
}
TEST_CASE("NormalNbrInfo from_json", "[NormalNbrInfo]")
{
	int id   = GENERATE(1, 2, 3);
	int rank = GENERATE(0, 1, 2);

	nlohmann::json j;
	j["type"]  = "NORMAL";
	j["ids"]   = {id};
	j["ranks"] = {rank};

	NormalNbrInfo<3> info = j.get<NormalNbrInfo<3>>();
	CHECK(info.id == id);
	CHECK(info.rank == rank);
}