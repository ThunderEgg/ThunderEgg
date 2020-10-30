#include "catch.hpp"
#include <ThunderEgg/PatchInfo.h>
using namespace std;
using namespace ThunderEgg;
TEST_CASE("FineNbrInfo Serialization/Deserialization", "[FineNbrInfo]")
{
	FineNbrInfo<3> info;
	info.ids[0]   = 1;
	info.ids[1]   = 2;
	info.ids[2]   = 3;
	info.ids[3]   = 4;
	info.ranks[0] = 9;
	info.ranks[1] = 8;
	info.ranks[2] = 7;
	info.ranks[3] = 6;
	// serialize and then deserialize
	char *buff = new char[info.serialize(nullptr)];
	info.serialize(buff);
	FineNbrInfo<3> out;
	out.deserialize(buff);
	delete[] buff;
	REQUIRE(out.ids[0] == 1);
	REQUIRE(out.ids[1] == 2);
	REQUIRE(out.ids[2] == 3);
	REQUIRE(out.ids[3] == 4);
	REQUIRE(out.ranks[0] == 9);
	REQUIRE(out.ranks[1] == 8);
	REQUIRE(out.ranks[2] == 7);
	REQUIRE(out.ranks[3] == 6);
}
TEST_CASE("FineNbrInfo to_json", "[FineNbrInfo]")
{
	FineNbrInfo<3> info;
	info.ids[0]   = GENERATE(1, 2);
	info.ids[1]   = GENERATE(1, 2);
	info.ids[2]   = GENERATE(1, 2);
	info.ids[3]   = GENERATE(1, 2);
	info.ranks[0] = GENERATE(0, 1);
	info.ranks[1] = GENERATE(0, 1);
	info.ranks[2] = GENERATE(0, 1);
	info.ranks[3] = GENERATE(0, 1);

	nlohmann::json j = info;

	CHECK(j["type"] == "FINE");
	REQUIRE(j["ids"].is_array());
	REQUIRE(j["ids"].size() == 4);
	CHECK(j["ids"][0] == info.ids[0]);
	CHECK(j["ids"][1] == info.ids[1]);
	CHECK(j["ids"][2] == info.ids[2]);
	CHECK(j["ids"][3] == info.ids[3]);
	REQUIRE(j["ranks"].is_array());
	REQUIRE(j["ranks"].size() == 4);
	CHECK(j["ranks"][0] == info.ranks[0]);
	CHECK(j["ranks"][1] == info.ranks[1]);
	CHECK(j["ranks"][2] == info.ranks[2]);
	CHECK(j["ranks"][3] == info.ranks[3]);
}
TEST_CASE("FineNbrInfo from_json", "[FineNbrInfo]")
{
	int id1   = GENERATE(1, 2);
	int id2   = GENERATE(1, 2);
	int id3   = GENERATE(1, 2);
	int id4   = GENERATE(1, 2);
	int rank1 = GENERATE(0, 1);
	int rank2 = GENERATE(0, 1);
	int rank3 = GENERATE(0, 1);
	int rank4 = GENERATE(0, 1);

	nlohmann::json j;
	j["type"]  = "NORMAL";
	j["ids"]   = {id1, id2, id3, id4};
	j["ranks"] = {rank1, rank2, rank3, rank4};

	FineNbrInfo<3> info = j.get<FineNbrInfo<3>>();
	CHECK(info.ids[0] == id1);
	CHECK(info.ids[1] == id2);
	CHECK(info.ids[2] == id3);
	CHECK(info.ids[3] == id4);
	CHECK(info.ranks[0] == rank1);
	CHECK(info.ranks[1] == rank2);
	CHECK(info.ranks[2] == rank3);
	CHECK(info.ranks[3] == rank4);
}