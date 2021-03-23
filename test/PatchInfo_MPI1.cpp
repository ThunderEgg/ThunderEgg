#include <ThunderEgg/PatchInfo.h>

#include <catch2/catch_test_macros.hpp>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("PatchInfo Serialization/Deserialization", "[PatchInfo]")
{
	PatchInfo<3> *d_ptr = new PatchInfo<3>;
	PatchInfo<3> &d     = *d_ptr;
	d.id                = 0;
	d.nbr_info[Side<3>::north().getIndex()].reset(new NormalNbrInfo<3>(1));
	d.nbr_info[Side<3>::east().getIndex()].reset(new CoarseNbrInfo<3>(2, Orthant<2>::nw()));
	d.nbr_info[Side<3>::south().getIndex()].reset(new FineNbrInfo<3>({3, 4, 5, 6}));

	// serialize and then deserialize
	char *buff = new char[d.serialize(nullptr)];
	d.serialize(buff);
	delete d_ptr;
	PatchInfo<3> out;
	out.deserialize(buff);
	delete[] buff;

	// check that deserialized version has the same information
	REQUIRE(out.id == 0);

	REQUIRE(!out.hasNbr(Side<3>::west()));

	REQUIRE(out.hasNbr(Side<3>::east()));
	REQUIRE(out.getNbrType(Side<3>::east()) == NbrType::Coarse);
	REQUIRE(out.getCoarseNbrInfo(Side<3>::east()).id == 2);
	REQUIRE(out.getCoarseNbrInfo(Side<3>::east()).orth_on_coarse == Orthant<2>::nw());

	REQUIRE(out.hasNbr(Side<3>::south()));
	REQUIRE(out.getNbrType(Side<3>::south()) == NbrType::Fine);
	REQUIRE(out.getFineNbrInfo(Side<3>::south()).ids[0] == 3);
	REQUIRE(out.getFineNbrInfo(Side<3>::south()).ids[1] == 4);
	REQUIRE(out.getFineNbrInfo(Side<3>::south()).ids[2] == 5);
	REQUIRE(out.getFineNbrInfo(Side<3>::south()).ids[3] == 6);

	REQUIRE(out.hasNbr(Side<3>::north()));
	REQUIRE(out.getNbrType(Side<3>::north()) == NbrType::Normal);
	REQUIRE(out.getNormalNbrInfo(Side<3>::north()).id == 1);

	REQUIRE(!out.hasNbr(Side<3>::bottom()));
	REQUIRE(!out.hasNbr(Side<3>::top()));
}
TEST_CASE("PatchInfo Default Values", "[PatchInfo]")
{
	PatchInfo<3> pinfo;
	CHECK(pinfo.id == 0);
	CHECK(pinfo.local_index == 0);
	CHECK(pinfo.global_index == 0);
	CHECK(pinfo.refine_level == -1);
	CHECK(pinfo.parent_id == -1);
	CHECK(pinfo.parent_rank == -1);
	for (int child_id : pinfo.child_ids) {
		CHECK(child_id == -1);
	}
	for (int child_rank : pinfo.child_ids) {
		CHECK(child_rank == -1);
	}
	CHECK(pinfo.num_ghost_cells == 0);
	CHECK(pinfo.rank == -1);
	CHECK(pinfo.orth_on_parent == Orthant<3>::null());
	CHECK(pinfo.neumann.to_ulong() == 0);
	for (int n : pinfo.ns) {
		CHECK(n == 1);
	}
	for (double start : pinfo.starts) {
		CHECK(start == 0);
	}
	for (double spacing : pinfo.spacings) {
		CHECK(spacing == 1);
	}
	for (auto nbr_info : pinfo.nbr_info) {
		CHECK(nbr_info == nullptr);
	}
	for (int idx : pinfo.bc_local_index) {
		CHECK(idx == -1);
	}
	for (int idx : pinfo.bc_global_index) {
		CHECK(idx == -1);
	}
}
TEST_CASE("PatchInfo to_json no children", "[PatchInfo]")
{
	PatchInfo<3> d;
	d.id             = 9;
	d.rank           = 0;
	d.parent_id      = 2;
	d.parent_rank    = 3;
	d.orth_on_parent = Orthant<3>::tnw();
	d.starts         = {1, 2, 3};
	d.spacings       = {0.1, 0.2, 0.3};
	d.ns             = {10, 20, 30};
	d.nbr_info[Side<3>::north().getIndex()].reset(new NormalNbrInfo<3>(1));
	d.nbr_info[Side<3>::east().getIndex()].reset(new CoarseNbrInfo<3>(2, Orthant<2>::nw()));
	d.nbr_info[Side<3>::south().getIndex()].reset(new FineNbrInfo<3>({3, 4, 5, 6}));

	nlohmann::json j = d;

	CHECK(j["id"] == d.id);
	CHECK(j["parent_id"] == d.parent_id);
	CHECK(j["parent_rank"] == d.parent_rank);
	CHECK(j["orth_on_parent"] == "TNW");
	CHECK(j["rank"] == d.rank);
	CHECK(j["child_ids"] == nullptr);
	CHECK(j["child_ranks"] == nullptr);

	REQUIRE(j["starts"].is_array());
	REQUIRE(j["starts"].size() == 3);
	CHECK(j["starts"][0] == d.starts[0]);
	CHECK(j["starts"][1] == d.starts[1]);
	CHECK(j["starts"][2] == d.starts[2]);

	REQUIRE(j["lengths"].is_array());
	REQUIRE(j["lengths"].size() == 3);
	CHECK(j["lengths"][0] == d.spacings[0] * d.ns[0]);
	CHECK(j["lengths"][1] == d.spacings[1] * d.ns[1]);
	CHECK(j["lengths"][2] == d.spacings[2] * d.ns[2]);

	REQUIRE(j["nbrs"].is_array());
	REQUIRE(j["nbrs"].size() == 3);

	CHECK(j["nbrs"][0]["type"] == "COARSE");
	CHECK(j["nbrs"][0]["side"] == "EAST");

	CHECK(j["nbrs"][1]["type"] == "FINE");
	CHECK(j["nbrs"][1]["side"] == "SOUTH");

	CHECK(j["nbrs"][2]["type"] == "NORMAL");
	CHECK(j["nbrs"][2]["side"] == "NORTH");
}
TEST_CASE("PatchInfo to_json no children no neighbors", "[PatchInfo]")
{
	PatchInfo<3> d;
	d.id          = 9;
	d.rank        = 0;
	d.parent_id   = 2;
	d.parent_rank = 3;
	d.starts      = {1, 2, 3};
	d.spacings    = {0.1, 0.2, 0.3};
	d.ns          = {10, 20, 30};

	nlohmann::json j = d;

	CHECK(j["id"] == d.id);
	CHECK(j["parent_id"] == d.parent_id);
	CHECK(j["parent_rank"] == d.parent_rank);
	CHECK(j["rank"] == d.rank);
	CHECK(j["child_ids"] == nullptr);
	CHECK(j["child_ranks"] == nullptr);
	CHECK(j["orth_on_parent"] == nullptr);

	REQUIRE(j["starts"].is_array());
	REQUIRE(j["starts"].size() == 3);
	CHECK(j["starts"][0] == d.starts[0]);
	CHECK(j["starts"][1] == d.starts[1]);
	CHECK(j["starts"][2] == d.starts[2]);

	REQUIRE(j["lengths"].is_array());
	REQUIRE(j["lengths"].size() == 3);
	CHECK(j["lengths"][0] == d.spacings[0] * d.ns[0]);
	CHECK(j["lengths"][1] == d.spacings[1] * d.ns[1]);
	CHECK(j["lengths"][2] == d.spacings[2] * d.ns[2]);

	REQUIRE(j["nbrs"].is_array());
	REQUIRE(j["nbrs"].size() == 0);
}
TEST_CASE("PatchInfo to_json with children", "[PatchInfo]")
{
	PatchInfo<3> d;
	d.id          = 9;
	d.rank        = 0;
	d.parent_id   = 2;
	d.parent_rank = 3;
	d.starts      = {1, 2, 3};
	d.spacings    = {0.1, 0.2, 0.3};
	d.ns          = {10, 20, 30};
	d.child_ids   = {3, 4, 5, 6, 7, 8, 9, 10};
	d.child_ranks = {1, 2, 3, 4, 5, 6, 7, 8};
	d.nbr_info[Side<3>::north().getIndex()].reset(new NormalNbrInfo<3>(1));
	d.nbr_info[Side<3>::east().getIndex()].reset(new CoarseNbrInfo<3>(2, Orthant<2>::nw()));
	d.nbr_info[Side<3>::south().getIndex()].reset(new FineNbrInfo<3>({3, 4, 5, 6}));

	nlohmann::json j = d;

	CHECK(j["id"] == d.id);
	CHECK(j["parent_id"] == d.parent_id);
	CHECK(j["parent_rank"] == d.parent_rank);
	CHECK(j["rank"] == d.rank);

	REQUIRE(j["child_ids"].is_array());
	REQUIRE(j["child_ids"].size() == 8);
	CHECK(j["child_ids"][0] == d.child_ids[0]);
	CHECK(j["child_ids"][1] == d.child_ids[1]);
	CHECK(j["child_ids"][2] == d.child_ids[2]);
	CHECK(j["child_ids"][3] == d.child_ids[3]);
	CHECK(j["child_ids"][4] == d.child_ids[4]);
	CHECK(j["child_ids"][5] == d.child_ids[5]);
	CHECK(j["child_ids"][6] == d.child_ids[6]);
	CHECK(j["child_ids"][7] == d.child_ids[7]);

	REQUIRE(j["child_ranks"].is_array());
	REQUIRE(j["child_ranks"].size() == 8);
	CHECK(j["child_ranks"][0] == d.child_ranks[0]);
	CHECK(j["child_ranks"][1] == d.child_ranks[1]);
	CHECK(j["child_ranks"][2] == d.child_ranks[2]);
	CHECK(j["child_ranks"][3] == d.child_ranks[3]);
	CHECK(j["child_ranks"][4] == d.child_ranks[4]);
	CHECK(j["child_ranks"][5] == d.child_ranks[5]);
	CHECK(j["child_ranks"][6] == d.child_ranks[6]);
	CHECK(j["child_ranks"][7] == d.child_ranks[7]);

	REQUIRE(j["starts"].is_array());
	REQUIRE(j["starts"].size() == 3);
	CHECK(j["starts"][0] == d.starts[0]);
	CHECK(j["starts"][1] == d.starts[1]);
	CHECK(j["starts"][2] == d.starts[2]);

	REQUIRE(j["lengths"].is_array());
	REQUIRE(j["lengths"].size() == 3);
	CHECK(j["lengths"][0] == d.spacings[0] * d.ns[0]);
	CHECK(j["lengths"][1] == d.spacings[1] * d.ns[1]);
	CHECK(j["lengths"][2] == d.spacings[2] * d.ns[2]);

	REQUIRE(j["nbrs"].is_array());
	REQUIRE(j["nbrs"].size() == 3);

	CHECK(j["nbrs"][0]["type"] == "COARSE");
	CHECK(j["nbrs"][0]["side"] == "EAST");

	CHECK(j["nbrs"][1]["type"] == "FINE");
	CHECK(j["nbrs"][1]["side"] == "SOUTH");

	CHECK(j["nbrs"][2]["type"] == "NORMAL");
	CHECK(j["nbrs"][2]["side"] == "NORTH");
}
TEST_CASE("PatchInfo from_json no children", "[PatchInfo]")
{
	nlohmann::json j;
	j["id"]          = 9;
	j["rank"]        = 3;
	j["parent_id"]   = 2;
	j["parent_rank"] = 3;
	j["starts"]      = {1, 2, 3};
	j["lengths"]     = {10, 20, 30};
	j["nbrs"]
	= {NormalNbrInfo<3>(1), CoarseNbrInfo<3>(2, Orthant<2>::nw()), FineNbrInfo<3>({3, 4, 5, 6})};
	j["nbrs"][0]["side"] = "NORTH";
	j["nbrs"][1]["side"] = "EAST";
	j["nbrs"][2]["side"] = "SOUTH";

	PatchInfo<3> d = j.get<PatchInfo<3>>();
	CHECK(d.id == 9);
	CHECK(d.rank == 3);
	CHECK(d.parent_id == 2);
	CHECK(d.parent_rank == 3);
	CHECK(d.orth_on_parent == Orthant<3>::null());
	CHECK(d.starts[0] == 1);
	CHECK(d.starts[1] == 2);
	CHECK(d.starts[2] == 3);
	CHECK(d.spacings[0] == 10);
	CHECK(d.spacings[1] == 20);
	CHECK(d.spacings[2] == 30);
	CHECK(d.ns[0] == 1);
	CHECK(d.ns[1] == 1);
	CHECK(d.ns[2] == 1);
	CHECK_FALSE(d.hasNbr(Side<3>::west()));
	CHECK(d.hasNbr(Side<3>::east()));
	CHECK(d.getNbrType(Side<3>::east()) == NbrType::Coarse);
	CHECK(d.hasNbr(Side<3>::south()));
	CHECK(d.getNbrType(Side<3>::south()) == NbrType::Fine);
	CHECK(d.hasNbr(Side<3>::north()));
	CHECK(d.getNbrType(Side<3>::north()) == NbrType::Normal);
	CHECK_FALSE(d.hasNbr(Side<3>::bottom()));
	CHECK_FALSE(d.hasNbr(Side<3>::top()));
}
TEST_CASE("PatchInfo from_json with children", "[PatchInfo]")
{
	nlohmann::json j;
	j["id"]             = 9;
	j["rank"]           = 3;
	j["parent_id"]      = 2;
	j["parent_rank"]    = 3;
	j["orth_on_parent"] = "TNW";
	j["starts"]         = {1, 2, 3};
	j["lengths"]        = {10, 20, 30};
	j["child_ids"]      = {1, 2, 3, 4, 5, 6, 7, 8};
	j["child_ranks"]    = {0, 1, 2, 3, 4, 5, 6, 7};
	j["nbrs"]
	= {NormalNbrInfo<3>(1), CoarseNbrInfo<3>(2, Orthant<2>::nw()), FineNbrInfo<3>({3, 4, 5, 6})};
	j["nbrs"][0]["side"] = "NORTH";
	j["nbrs"][1]["side"] = "EAST";
	j["nbrs"][2]["side"] = "SOUTH";

	PatchInfo<3> d = j.get<PatchInfo<3>>();
	CHECK(d.id == 9);
	CHECK(d.rank == 3);
	CHECK(d.parent_id == 2);
	CHECK(d.parent_rank == 3);
	CHECK(d.orth_on_parent == Orthant<3>::tnw());
	CHECK(d.starts[0] == 1);
	CHECK(d.starts[1] == 2);
	CHECK(d.starts[2] == 3);
	CHECK(d.spacings[0] == 10);
	CHECK(d.spacings[1] == 20);
	CHECK(d.spacings[2] == 30);
	CHECK(d.ns[0] == 1);
	CHECK(d.ns[1] == 1);
	CHECK(d.ns[2] == 1);
	CHECK(d.child_ids[0] == 1);
	CHECK(d.child_ids[1] == 2);
	CHECK(d.child_ids[2] == 3);
	CHECK(d.child_ids[3] == 4);
	CHECK(d.child_ids[4] == 5);
	CHECK(d.child_ids[5] == 6);
	CHECK(d.child_ids[6] == 7);
	CHECK(d.child_ids[7] == 8);
	CHECK(d.child_ranks[0] == 0);
	CHECK(d.child_ranks[1] == 1);
	CHECK(d.child_ranks[2] == 2);
	CHECK(d.child_ranks[3] == 3);
	CHECK(d.child_ranks[4] == 4);
	CHECK(d.child_ranks[5] == 5);
	CHECK(d.child_ranks[6] == 6);
	CHECK(d.child_ranks[7] == 7);
	CHECK_FALSE(d.hasNbr(Side<3>::west()));
	CHECK(d.hasNbr(Side<3>::east()));
	CHECK(d.getNbrType(Side<3>::east()) == NbrType::Coarse);
	CHECK(d.hasNbr(Side<3>::south()));
	CHECK(d.getNbrType(Side<3>::south()) == NbrType::Fine);
	CHECK(d.hasNbr(Side<3>::north()));
	CHECK(d.getNbrType(Side<3>::north()) == NbrType::Normal);
	CHECK_FALSE(d.hasNbr(Side<3>::bottom()));
	CHECK_FALSE(d.hasNbr(Side<3>::top()));
}