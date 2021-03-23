#include <ThunderEgg/NbrType.h>

#include <catch2/catch_test_macros.hpp>

using namespace ThunderEgg;
TEST_CASE("NbrType to_json", "[NbrType]")
{
	nlohmann::json j;
	j["normal"] = NbrType::Normal;
	j["coarse"] = NbrType::Coarse;
	j["fine"]   = NbrType::Fine;
	CHECK(j["normal"] == "NORMAL");
	CHECK(j["coarse"] == "COARSE");
	CHECK(j["fine"] == "FINE");
}
TEST_CASE("NbrType from_json", "[NbrType]")
{
	nlohmann::json j;
	j["normal"] = "NORMAL";
	j["coarse"] = "COARSE";
	j["fine"]   = "FINE";
	CHECK(j["normal"].get<NbrType>() == NbrType::Normal);
	CHECK(j["coarse"].get<NbrType>() == NbrType::Coarse);
	CHECK(j["fine"].get<NbrType>() == NbrType::Fine);
}