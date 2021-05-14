#include <ThunderEgg/PatchInfo.h>

#include <catch2/catch_test_macros.hpp>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("PatchInfo getNbrIds NormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<3> *nbr_info = new NormalNbrInfo<3>(2);
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int> ids = pinfo.getNbrIds();

	REQUIRE(ids.size() == 1);
	CHECK(ids[0] == 2);
}
TEST_CASE("PatchInfo getNbrRanks NormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<3> *nbr_info = new NormalNbrInfo<3>(2);
	nbr_info->rank             = 3;
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int> ranks = pinfo.getNbrRanks();

	REQUIRE(ranks.size() == 1);
	CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists NormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<3> *nbr_info = new NormalNbrInfo<3>(2);
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[2] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_index == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist NormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<3> *nbr_info = new NormalNbrInfo<3>(2);
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[59] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_index == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists NormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<3> *nbr_info = new NormalNbrInfo<3>(2);
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_global_index_map;
	id_to_global_index_map[2] = 30;
	pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

	CHECK(nbr_info->global_index == 30);
}
TEST_CASE("PatchInfo getNbrIds CoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<3> *nbr_info = new CoarseNbrInfo<3>(2, Orthant<2>::sw());
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int> ids = pinfo.getNbrIds();

	REQUIRE(ids.size() == 1);
	CHECK(ids[0] == 2);
}
TEST_CASE("PatchInfo getNbrRanks CoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<3> *nbr_info = new CoarseNbrInfo<3>(2, Orthant<2>::sw());
	nbr_info->rank             = 3;
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int> ranks = pinfo.getNbrRanks();

	REQUIRE(ranks.size() == 1);
	CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists CoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<3> *nbr_info = new CoarseNbrInfo<3>(2, Orthant<2>::sw());
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[2] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_index == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist CoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<3> *nbr_info = new CoarseNbrInfo<3>(2, Orthant<2>::sw());
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[59] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_index == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists CoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<3> *nbr_info = new CoarseNbrInfo<3>(2, Orthant<2>::sw());
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_global_index_map;
	id_to_global_index_map[2] = 30;
	pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

	CHECK(nbr_info->global_index == 30);
}
TEST_CASE("PatchInfo getNbrIds FineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<3> *nbr_info = new FineNbrInfo<3>({1, 2, 3, 4});
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int> ids = pinfo.getNbrIds();

	REQUIRE(ids.size() == 4);
	CHECK(ids[0] == 1);
	CHECK(ids[1] == 2);
	CHECK(ids[2] == 3);
	CHECK(ids[3] == 4);
}
TEST_CASE("PatchInfo getNbrRanks FineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<3> *nbr_info = new FineNbrInfo<3>({1, 2, 3, 4});
	nbr_info->ranks          = {3, 4, 6, 7};
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int> ranks = pinfo.getNbrRanks();

	REQUIRE(ranks.size() == 4);
	CHECK(ranks[0] == 3);
	CHECK(ranks[1] == 4);
	CHECK(ranks[2] == 6);
	CHECK(ranks[3] == 7);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists FineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<3> *nbr_info = new FineNbrInfo<3>({2, 3, 4, 5});
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[2] = 30;
	id_to_local_index_map[3] = 31;
	id_to_local_index_map[4] = 32;
	id_to_local_index_map[5] = 33;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_indexes[0] == 30);
	CHECK(nbr_info->local_indexes[1] == 31);
	CHECK(nbr_info->local_indexes[2] == 32);
	CHECK(nbr_info->local_indexes[3] == 33);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist FineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<3> *nbr_info = new FineNbrInfo<3>({2, 3, 4, 5});
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[2]  = 30;
	id_to_local_index_map[3]  = 31;
	id_to_local_index_map[4]  = 32;
	id_to_local_index_map[59] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_indexes[0] == 30);
	CHECK(nbr_info->local_indexes[1] == 31);
	CHECK(nbr_info->local_indexes[2] == 32);
	CHECK(nbr_info->local_indexes[3] == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists FineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<3> *nbr_info = new FineNbrInfo<3>({2, 3, 4, 5});
	pinfo.nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_global_index_map;
	id_to_global_index_map[2] = 30;
	id_to_global_index_map[3] = 31;
	id_to_global_index_map[4] = 32;
	id_to_global_index_map[5] = 33;
	pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

	CHECK(nbr_info->global_indexes[0] == 30);
	CHECK(nbr_info->global_indexes[1] == 31);
	CHECK(nbr_info->global_indexes[2] == 32);
	CHECK(nbr_info->global_indexes[3] == 33);
}
TEST_CASE("PatchInfo getNbrIds EdgeNormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<2> *nbr_info = new NormalNbrInfo<2>(2);
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int> ids = pinfo.getNbrIds();

	REQUIRE(ids.size() == 1);
	CHECK(ids[0] == 2);
}
TEST_CASE("PatchInfo getNbrRanks EdgeNormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<2> *nbr_info = new NormalNbrInfo<2>(2);
	nbr_info->rank             = 3;
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int> ranks = pinfo.getNbrRanks();

	REQUIRE(ranks.size() == 1);
	CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists EdgeNormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<2> *nbr_info = new NormalNbrInfo<2>(2);
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[2] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_index == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist EdgeNormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<2> *nbr_info = new NormalNbrInfo<2>(2);
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[59] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_index == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists EdgeNormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<2> *nbr_info = new NormalNbrInfo<2>(2);
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_global_index_map;
	id_to_global_index_map[2] = 30;
	pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

	CHECK(nbr_info->global_index == 30);
}
TEST_CASE("PatchInfo getNbrIds EdgeCoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<2> *nbr_info = new CoarseNbrInfo<2>(2, Orthant<1>::lower());
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int> ids = pinfo.getNbrIds();

	REQUIRE(ids.size() == 1);
	CHECK(ids[0] == 2);
}
TEST_CASE("PatchInfo getNbrRanks EdgeCoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<2> *nbr_info = new CoarseNbrInfo<2>(2, Orthant<1>::lower());
	nbr_info->rank             = 3;
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int> ranks = pinfo.getNbrRanks();

	REQUIRE(ranks.size() == 1);
	CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists EdgeCoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<2> *nbr_info = new CoarseNbrInfo<2>(2, Orthant<1>::lower());
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[2] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_index == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist EdgeCoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<2> *nbr_info = new CoarseNbrInfo<2>(2, Orthant<1>::lower());
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[59] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_index == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists EdgeCoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<2> *nbr_info = new CoarseNbrInfo<2>(2, Orthant<1>::lower());
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_global_index_map;
	id_to_global_index_map[2] = 30;
	pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

	CHECK(nbr_info->global_index == 30);
}
TEST_CASE("PatchInfo getNbrIds EdgeFineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<2> *nbr_info = new FineNbrInfo<2>({1, 2});
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int> ids = pinfo.getNbrIds();

	REQUIRE(ids.size() == 2);
	CHECK(ids[0] == 1);
	CHECK(ids[1] == 2);
}
TEST_CASE("PatchInfo getNbrRanks EdgeFineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<2> *nbr_info = new FineNbrInfo<2>({1, 2});
	nbr_info->ranks          = {3, 4};
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int> ranks = pinfo.getNbrRanks();

	REQUIRE(ranks.size() == 2);
	CHECK(ranks[0] == 3);
	CHECK(ranks[1] == 4);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists EdgeFineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<2> *nbr_info = new FineNbrInfo<2>({2, 3});
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[2] = 30;
	id_to_local_index_map[3] = 31;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_indexes[0] == 30);
	CHECK(nbr_info->local_indexes[1] == 31);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist EdgeFineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<2> *nbr_info = new FineNbrInfo<2>({2, 3});
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[2]  = 30;
	id_to_local_index_map[59] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_indexes[0] == 30);
	CHECK(nbr_info->local_indexes[1] == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists EdgeFineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<2> *nbr_info = new FineNbrInfo<2>({2, 3});
	pinfo.edge_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_global_index_map;
	id_to_global_index_map[2] = 30;
	id_to_global_index_map[3] = 31;
	pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

	CHECK(nbr_info->global_indexes[0] == 30);
	CHECK(nbr_info->global_indexes[1] == 31);
}
TEST_CASE("PatchInfo getNbrIds CornerNormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<1> *nbr_info = new NormalNbrInfo<1>(2);
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int> ids = pinfo.getNbrIds();

	REQUIRE(ids.size() == 1);
	CHECK(ids[0] == 2);
}
TEST_CASE("PatchInfo getNbrRanks CornerNormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<1> *nbr_info = new NormalNbrInfo<1>(2);
	nbr_info->rank             = 3;
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int> ranks = pinfo.getNbrRanks();

	REQUIRE(ranks.size() == 1);
	CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists CornerNormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<1> *nbr_info = new NormalNbrInfo<1>(2);
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[2] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_index == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist CornerNormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<1> *nbr_info = new NormalNbrInfo<1>(2);
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[59] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_index == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists CornerNormalNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	NormalNbrInfo<1> *nbr_info = new NormalNbrInfo<1>(2);
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_global_index_map;
	id_to_global_index_map[2] = 30;
	pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

	CHECK(nbr_info->global_index == 30);
}
TEST_CASE("PatchInfo getNbrIds CornerCoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<1> *nbr_info = new CoarseNbrInfo<1>(2, Orthant<0>::null());
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int> ids = pinfo.getNbrIds();

	REQUIRE(ids.size() == 1);
	CHECK(ids[0] == 2);
}
TEST_CASE("PatchInfo getNbrRanks CornerCoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<1> *nbr_info = new CoarseNbrInfo<1>(2, Orthant<0>::null());
	nbr_info->rank             = 3;
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int> ranks = pinfo.getNbrRanks();

	REQUIRE(ranks.size() == 1);
	CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists CornerCoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<1> *nbr_info = new CoarseNbrInfo<1>(2, Orthant<0>::null());
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[2] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_index == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist CornerCoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<1> *nbr_info = new CoarseNbrInfo<1>(2, Orthant<0>::null());
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[59] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_index == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists CornerCoarseNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>      pinfo;
	CoarseNbrInfo<1> *nbr_info = new CoarseNbrInfo<1>(2, Orthant<0>::null());
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_global_index_map;
	id_to_global_index_map[2] = 30;
	pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

	CHECK(nbr_info->global_index == 30);
}
TEST_CASE("PatchInfo getNbrIds CornerFineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<1> *nbr_info = new FineNbrInfo<1>({1});
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int> ids = pinfo.getNbrIds();

	REQUIRE(ids.size() == 1);
	CHECK(ids[0] == 1);
}
TEST_CASE("PatchInfo getNbrRanks CornerFineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<1> *nbr_info = new FineNbrInfo<1>({1});
	nbr_info->ranks          = {3};
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int> ranks = pinfo.getNbrRanks();

	REQUIRE(ranks.size() == 1);
	CHECK(ranks[0] == 3);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes exists CornerFineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<1> *nbr_info = new FineNbrInfo<1>({2});
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[2] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_indexes[0] == 30);
}
TEST_CASE("PatchInfo setNeighborLocalIndexes does not exist CornerFineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<1> *nbr_info = new FineNbrInfo<1>({2});
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_local_index_map;
	id_to_local_index_map[59] = 30;
	pinfo.setNeighborLocalIndexes(id_to_local_index_map);

	CHECK(nbr_info->local_indexes[0] == -1);
}
TEST_CASE("PatchInfo setNeighborGlobalIndexes exists CornerFineNbrInfo", "[PatchInfo]")
{
	PatchInfo<3>    pinfo;
	FineNbrInfo<1> *nbr_info = new FineNbrInfo<1>({2});
	pinfo.corner_nbr_info[0].reset(nbr_info);
	deque<int>    ranks = pinfo.getNbrRanks();
	map<int, int> id_to_global_index_map;
	id_to_global_index_map[2] = 30;
	pinfo.setNeighborGlobalIndexes(id_to_global_index_map);

	CHECK(nbr_info->global_indexes[0] == 30);
}
TEST_CASE("PatchInfo Serialization/Deserialization", "[PatchInfo]")
{
	PatchInfo<3> *d_ptr = new PatchInfo<3>;
	PatchInfo<3> &d     = *d_ptr;
	d.id                = 0;
	d.nbr_info[Side<3>::north().getIndex()].reset(new NormalNbrInfo<3>(1));
	d.nbr_info[Side<3>::east().getIndex()].reset(new CoarseNbrInfo<3>(2, Orthant<2>::nw()));
	d.nbr_info[Side<3>::south().getIndex()].reset(new FineNbrInfo<3>({3, 4, 5, 6}));
	d.corner_nbr_info[Corner<3>::bsw().getIndex()].reset(new NormalNbrInfo<1>(1));
	d.corner_nbr_info[Corner<3>::tse().getIndex()].reset(new CoarseNbrInfo<1>(2, Orthant<0>::null()));
	d.corner_nbr_info[Corner<3>::bnw().getIndex()].reset(new FineNbrInfo<1>({1}));
	d.edge_nbr_info[Edge::sw().getIndex()].reset(new NormalNbrInfo<2>(1));
	d.edge_nbr_info[Edge::bn().getIndex()].reset(new CoarseNbrInfo<2>(2, Orthant<1>::lower()));
	d.edge_nbr_info[Edge::tw().getIndex()].reset(new FineNbrInfo<2>({1, 2}));

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

	// Corners

	REQUIRE(out.hasNbr(Corner<3>::bsw()));
	REQUIRE(out.getNbrType(Corner<3>::bsw()) == NbrType::Normal);
	REQUIRE(out.getNormalNbrInfo(Corner<3>::bsw()).id == 1);

	REQUIRE(!out.hasNbr(Corner<3>::bse()));

	REQUIRE(out.hasNbr(Corner<3>::bnw()));
	REQUIRE(out.getNbrType(Corner<3>::bnw()) == NbrType::Fine);
	REQUIRE(out.getFineNbrInfo(Corner<3>::bnw()).ids[0] == 1);

	REQUIRE(!out.hasNbr(Corner<3>::bne()));
	REQUIRE(!out.hasNbr(Corner<3>::tsw()));

	REQUIRE(out.hasNbr(Corner<3>::tse()));
	REQUIRE(out.getNbrType(Corner<3>::tse()) == NbrType::Coarse);
	REQUIRE(out.getCoarseNbrInfo(Corner<3>::tse()).id == 2);
	REQUIRE(out.getCoarseNbrInfo(Corner<3>::tse()).orth_on_coarse == Orthant<0>::null());

	REQUIRE(!out.hasNbr(Corner<3>::tnw()));
	REQUIRE(!out.hasNbr(Corner<3>::tne()));

	// Edges

	REQUIRE(!out.hasNbr(Edge::bs()));
	REQUIRE(!out.hasNbr(Edge::tn()));

	REQUIRE(out.hasNbr(Edge::bn()));
	REQUIRE(out.getNbrType(Edge::bn()) == NbrType::Coarse);
	REQUIRE(out.getCoarseNbrInfo(Edge::bn()).id == 2);
	REQUIRE(out.getCoarseNbrInfo(Edge::bn()).orth_on_coarse == Orthant<1>::lower());

	REQUIRE(!out.hasNbr(Edge::ts()));
	REQUIRE(!out.hasNbr(Edge::bw()));
	REQUIRE(!out.hasNbr(Edge::te()));
	REQUIRE(!out.hasNbr(Edge::be()));

	REQUIRE(out.hasNbr(Edge::tw()));
	REQUIRE(out.getNbrType(Edge::tw()) == NbrType::Fine);
	REQUIRE(out.getFineNbrInfo(Edge::tw()).ids[0] == 1);
	REQUIRE(out.getFineNbrInfo(Edge::tw()).ids[1] == 2);

	REQUIRE(out.hasNbr(Edge::sw()));
	REQUIRE(out.getNbrType(Edge::sw()) == NbrType::Normal);
	REQUIRE(out.getNormalNbrInfo(Edge::sw()).id == 1);

	REQUIRE(!out.hasNbr(Edge::ne()));
	REQUIRE(!out.hasNbr(Edge::se()));
	REQUIRE(!out.hasNbr(Edge::nw()));
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
	for (int n : pinfo.ns) {
		CHECK(n == 1);
	}
	for (double start : pinfo.starts) {
		CHECK(start == 0);
	}
	for (double spacing : pinfo.spacings) {
		CHECK(spacing == 1);
	}
	for (auto &nbr_info : pinfo.nbr_info) {
		CHECK(nbr_info == nullptr);
	}
	for (auto &nbr_info : pinfo.corner_nbr_info) {
		CHECK(nbr_info == nullptr);
	}
	for (auto &nbr_info : pinfo.edge_nbr_info) {
		CHECK(nbr_info == nullptr);
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
	d.corner_nbr_info[Corner<3>::bsw().getIndex()].reset(new NormalNbrInfo<1>(1));
	d.corner_nbr_info[Corner<3>::tse().getIndex()].reset(new CoarseNbrInfo<1>(2, Orthant<0>(0)));
	d.corner_nbr_info[Corner<3>::bnw().getIndex()].reset(new FineNbrInfo<1>({1}));
	d.edge_nbr_info[Edge::sw().getIndex()].reset(new NormalNbrInfo<2>(1));
	d.edge_nbr_info[Edge::bn().getIndex()].reset(new CoarseNbrInfo<2>(2, Orthant<1>::lower()));
	d.edge_nbr_info[Edge::tw().getIndex()].reset(new FineNbrInfo<2>({1, 2}));

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

	REQUIRE(j["corner_nbrs"].is_array());
	REQUIRE(j["corner_nbrs"].size() == 3);

	CHECK(j["corner_nbrs"][0]["type"] == "NORMAL");
	CHECK(j["corner_nbrs"][0]["corner"] == "BSW");

	CHECK(j["corner_nbrs"][1]["type"] == "FINE");
	CHECK(j["corner_nbrs"][1]["corner"] == "BNW");

	CHECK(j["corner_nbrs"][2]["type"] == "COARSE");
	CHECK(j["corner_nbrs"][2]["corner"] == "TSE");

	REQUIRE(j["edge_nbrs"].is_array());
	REQUIRE(j["edge_nbrs"].size() == 3);

	CHECK(j["edge_nbrs"][0]["type"] == "COARSE");
	CHECK(j["edge_nbrs"][0]["edge"] == "BN");

	CHECK(j["edge_nbrs"][1]["type"] == "FINE");
	CHECK(j["edge_nbrs"][1]["edge"] == "TW");

	CHECK(j["edge_nbrs"][2]["type"] == "NORMAL");
	CHECK(j["edge_nbrs"][2]["edge"] == "SW");
}
TEST_CASE("PatchInfo to_json no children no neighbors", "[PatchInfo]")
{
	PatchInfo<3> d;
	d.id           = 9;
	d.rank         = 0;
	d.parent_id    = 2;
	d.parent_rank  = 3;
	d.refine_level = 329;
	d.starts       = {1, 2, 3};
	d.spacings     = {0.1, 0.2, 0.3};
	d.ns           = {10, 20, 30};

	nlohmann::json j = d;

	CHECK(j["id"] == d.id);
	CHECK(j["parent_id"] == d.parent_id);
	CHECK(j["parent_rank"] == d.parent_rank);
	CHECK(j["rank"] == d.rank);
	CHECK(j["refine_level"] == 329);
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
	d.id           = 9;
	d.rank         = 0;
	d.parent_id    = 2;
	d.parent_rank  = 3;
	d.refine_level = 329;
	d.starts       = {1, 2, 3};
	d.spacings     = {0.1, 0.2, 0.3};
	d.ns           = {10, 20, 30};
	d.child_ids    = {3, 4, 5, 6, 7, 8, 9, 10};
	d.child_ranks  = {1, 2, 3, 4, 5, 6, 7, 8};
	d.nbr_info[Side<3>::north().getIndex()].reset(new NormalNbrInfo<3>(1));
	d.nbr_info[Side<3>::east().getIndex()].reset(new CoarseNbrInfo<3>(2, Orthant<2>::nw()));
	d.nbr_info[Side<3>::south().getIndex()].reset(new FineNbrInfo<3>({3, 4, 5, 6}));

	nlohmann::json j = d;

	CHECK(j["id"] == d.id);
	CHECK(j["parent_id"] == d.parent_id);
	CHECK(j["parent_rank"] == d.parent_rank);
	CHECK(j["rank"] == d.rank);
	CHECK(j["refine_level"] == 329);

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
	j["id"]                       = 9;
	j["rank"]                     = 3;
	j["refine_level"]             = 329;
	j["parent_id"]                = 2;
	j["parent_rank"]              = 3;
	j["starts"]                   = {1, 2, 3};
	j["lengths"]                  = {10, 20, 30};
	j["nbrs"]                     = {NormalNbrInfo<3>(1), CoarseNbrInfo<3>(2, Orthant<2>::nw()), FineNbrInfo<3>({3, 4, 5, 6})};
	j["nbrs"][0]["side"]          = "NORTH";
	j["nbrs"][1]["side"]          = "EAST";
	j["nbrs"][2]["side"]          = "SOUTH";
	j["corner_nbrs"]              = {NormalNbrInfo<1>(1), CoarseNbrInfo<1>(2, Orthant<0>(0)), FineNbrInfo<1>({1})};
	j["corner_nbrs"][0]["corner"] = "BSW";
	j["corner_nbrs"][1]["corner"] = "TSE";
	j["corner_nbrs"][2]["corner"] = "BNW";
	j["edge_nbrs"]                = {NormalNbrInfo<2>(1), CoarseNbrInfo<2>(2, Orthant<1>::lower()), FineNbrInfo<2>({1, 2})};
	j["edge_nbrs"][0]["edge"]     = "SW";
	j["edge_nbrs"][1]["edge"]     = "BN";
	j["edge_nbrs"][2]["edge"]     = "TW";

	PatchInfo<3>
	d = j.get<PatchInfo<3>>();
	CHECK(d.id == 9);
	CHECK(d.rank == 3);
	CHECK(d.refine_level == 329);
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

	CHECK(d.hasNbr(Corner<3>::bsw()));
	CHECK(d.getNbrType(Corner<3>::bsw()) == NbrType::Normal);
	CHECK_FALSE(d.hasNbr(Corner<3>::bse()));
	CHECK(d.hasNbr(Corner<3>::bnw()));
	CHECK(d.getNbrType(Corner<3>::bnw()) == NbrType::Fine);
	CHECK_FALSE(d.hasNbr(Corner<3>::bne()));
	CHECK_FALSE(d.hasNbr(Corner<3>::tsw()));
	CHECK(d.hasNbr(Corner<3>::tse()));
	CHECK(d.getNbrType(Corner<3>::tse()) == NbrType::Coarse);
	CHECK_FALSE(d.hasNbr(Corner<3>::tnw()));
	CHECK_FALSE(d.hasNbr(Corner<3>::tne()));

	CHECK_FALSE(d.hasNbr(Edge::bs()));
	CHECK_FALSE(d.hasNbr(Edge::tn()));
	CHECK(d.hasNbr(Edge::bn()));
	CHECK(d.getNbrType(Edge::bn()) == NbrType::Coarse);
	CHECK_FALSE(d.hasNbr(Edge::ts()));
	CHECK_FALSE(d.hasNbr(Edge::bw()));
	CHECK_FALSE(d.hasNbr(Edge::te()));
	CHECK_FALSE(d.hasNbr(Edge::be()));
	CHECK(d.hasNbr(Edge::tw()));
	CHECK(d.getNbrType(Edge::tw()) == NbrType::Fine);
	CHECK(d.hasNbr(Edge::sw()));
	CHECK(d.getNbrType(Edge::sw()) == NbrType::Normal);
	CHECK_FALSE(d.hasNbr(Edge::ne()));
	CHECK_FALSE(d.hasNbr(Edge::se()));
	CHECK_FALSE(d.hasNbr(Edge::nw()));
}
TEST_CASE("PatchInfo from_json with children", "[PatchInfo]")
{
	nlohmann::json j;
	j["id"]             = 9;
	j["rank"]           = 3;
	j["refine_level"]   = 329;
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
	CHECK(d.refine_level == 329);
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
TEST_CASE("PatchInfo copy constructor", "[PatchInfo]")
{
	PatchInfo<3> d;
	d.id              = 9;
	d.local_index     = 10;
	d.global_index    = 10;
	d.rank            = 0;
	d.parent_id       = 2;
	d.parent_rank     = 3;
	d.num_ghost_cells = 239;
	d.refine_level    = 329;
	d.starts          = {1, 2, 3};
	d.spacings        = {0.1, 0.2, 0.3};
	d.ns              = {10, 20, 30};
	d.child_ids       = {3, 4, 5, 6, 7, 8, 9, 10};
	d.child_ranks     = {1, 2, 3, 4, 5, 6, 7, 8};
	d.nbr_info[Side<3>::north().getIndex()].reset(new NormalNbrInfo<3>(1));
	d.nbr_info[Side<3>::east().getIndex()].reset(new CoarseNbrInfo<3>(2, Orthant<2>::nw()));
	d.nbr_info[Side<3>::south().getIndex()].reset(new FineNbrInfo<3>({3, 4, 5, 6}));
	d.corner_nbr_info[Corner<3>::bsw().getIndex()].reset(new NormalNbrInfo<1>(1));
	d.corner_nbr_info[Corner<3>::tse().getIndex()].reset(new CoarseNbrInfo<1>(2, Orthant<0>(0)));
	d.corner_nbr_info[Corner<3>::bnw().getIndex()].reset(new FineNbrInfo<1>({1}));
	d.edge_nbr_info[Edge::sw().getIndex()].reset(new NormalNbrInfo<2>(1));
	d.edge_nbr_info[Edge::bn().getIndex()].reset(new CoarseNbrInfo<2>(2, Orthant<1>::lower()));
	d.edge_nbr_info[Edge::tw().getIndex()].reset(new FineNbrInfo<2>({1, 2}));

	PatchInfo<3> d2(d);

	CHECK(d.id == d2.id);
	CHECK(d.local_index == d2.global_index);
	CHECK(d.rank == d2.rank);
	CHECK(d.parent_id == d2.parent_id);
	CHECK(d.parent_rank == d2.parent_rank);
	CHECK(d.num_ghost_cells == d2.num_ghost_cells);
	CHECK(d.refine_level == d2.refine_level);
	CHECK(d.starts == d2.starts);
	CHECK(d.spacings == d2.spacings);
	CHECK(d.ns == d2.ns);
	CHECK(d.child_ids == d2.child_ids);
	CHECK(d.child_ranks == d2.child_ranks);

	for (Side<3> s : Side<3>::getValues()) {
		REQUIRE(d.hasNbr(s) == d2.hasNbr(s));
		if (d.hasNbr(s)) {
			CHECK(d.nbr_info[s.getIndex()] != d2.nbr_info[s.getIndex()]);
			switch (d.getNbrType(s)) {
				case NbrType::Normal:
					CHECK(d.getNormalNbrInfo(s).id == d2.getNormalNbrInfo(s).id);
					break;
				case NbrType::Fine:
					CHECK(d.getFineNbrInfo(s).ids[0] == d2.getFineNbrInfo(s).ids[0]);
					break;
				case NbrType::Coarse:
					CHECK(d.getCoarseNbrInfo(s).id == d2.getCoarseNbrInfo(s).id);
					break;
			}
		}
	}
	for (Corner<3> c : Corner<3>::getValues()) {
		REQUIRE(d.hasNbr(c) == d2.hasNbr(c));
		if (d.hasNbr(c)) {
			CHECK(d.corner_nbr_info[c.getIndex()] != d2.corner_nbr_info[c.getIndex()]);
			switch (d.getNbrType(c)) {
				case NbrType::Normal:
					CHECK(d.getNormalNbrInfo(c).id == d2.getNormalNbrInfo(c).id);
					break;
				case NbrType::Fine:
					CHECK(d.getFineNbrInfo(c).ids[0] == d2.getFineNbrInfo(c).ids[0]);
					break;
				case NbrType::Coarse:
					CHECK(d.getCoarseNbrInfo(c).id == d2.getCoarseNbrInfo(c).id);
					break;
			}
		}
	}
	for (Edge c : Edge::getValues()) {
		REQUIRE(d.hasNbr(c) == d2.hasNbr(c));
		if (d.hasNbr(c)) {
			CHECK(d.edge_nbr_info[c.getIndex()] != d2.edge_nbr_info[c.getIndex()]);
			switch (d.getNbrType(c)) {
				case NbrType::Normal:
					CHECK(d.getNormalNbrInfo(c).id == d2.getNormalNbrInfo(c).id);
					break;
				case NbrType::Fine:
					CHECK(d.getFineNbrInfo(c).ids[0] == d2.getFineNbrInfo(c).ids[0]);
					break;
				case NbrType::Coarse:
					CHECK(d.getCoarseNbrInfo(c).id == d2.getCoarseNbrInfo(c).id);
					break;
			}
		}
	}
}
TEST_CASE("PatchInfo copy assignment", "[PatchInfo]")
{
	PatchInfo<3> d;
	d.id              = 9;
	d.local_index     = 10;
	d.global_index    = 10;
	d.rank            = 0;
	d.parent_id       = 2;
	d.parent_rank     = 3;
	d.num_ghost_cells = 239;
	d.refine_level    = 329;
	d.starts          = {1, 2, 3};
	d.spacings        = {0.1, 0.2, 0.3};
	d.ns              = {10, 20, 30};
	d.child_ids       = {3, 4, 5, 6, 7, 8, 9, 10};
	d.child_ranks     = {1, 2, 3, 4, 5, 6, 7, 8};
	d.nbr_info[Side<3>::north().getIndex()].reset(new NormalNbrInfo<3>(1));
	d.nbr_info[Side<3>::east().getIndex()].reset(new CoarseNbrInfo<3>(2, Orthant<2>::nw()));
	d.nbr_info[Side<3>::south().getIndex()].reset(new FineNbrInfo<3>({3, 4, 5, 6}));
	d.corner_nbr_info[Corner<3>::bsw().getIndex()].reset(new NormalNbrInfo<1>(1));
	d.corner_nbr_info[Corner<3>::tse().getIndex()].reset(new CoarseNbrInfo<1>(2, Orthant<0>(0)));
	d.corner_nbr_info[Corner<3>::bnw().getIndex()].reset(new FineNbrInfo<1>({1}));
	d.edge_nbr_info[Edge::sw().getIndex()].reset(new NormalNbrInfo<2>(1));
	d.edge_nbr_info[Edge::bn().getIndex()].reset(new CoarseNbrInfo<2>(2, Orthant<1>::lower()));
	d.edge_nbr_info[Edge::tw().getIndex()].reset(new FineNbrInfo<2>({1, 2}));

	PatchInfo<3> d2;
	d2 = d;

	CHECK(d.id == d2.id);
	CHECK(d.local_index == d2.global_index);
	CHECK(d.rank == d2.rank);
	CHECK(d.parent_id == d2.parent_id);
	CHECK(d.parent_rank == d2.parent_rank);
	CHECK(d.num_ghost_cells == d2.num_ghost_cells);
	CHECK(d.refine_level == d2.refine_level);
	CHECK(d.starts == d2.starts);
	CHECK(d.spacings == d2.spacings);
	CHECK(d.ns == d2.ns);
	CHECK(d.child_ids == d2.child_ids);
	CHECK(d.child_ranks == d2.child_ranks);

	for (Side<3> s : Side<3>::getValues()) {
		REQUIRE(d.hasNbr(s) == d2.hasNbr(s));
		if (d.hasNbr(s)) {
			CHECK(d.nbr_info[s.getIndex()] != d2.nbr_info[s.getIndex()]);
			switch (d.getNbrType(s)) {
				case NbrType::Normal:
					CHECK(d.getNormalNbrInfo(s).id == d2.getNormalNbrInfo(s).id);
					break;
				case NbrType::Fine:
					CHECK(d.getFineNbrInfo(s).ids[0] == d2.getFineNbrInfo(s).ids[0]);
					break;
				case NbrType::Coarse:
					CHECK(d.getCoarseNbrInfo(s).id == d2.getCoarseNbrInfo(s).id);
					break;
			}
		}
	}
	for (Corner<3> c : Corner<3>::getValues()) {
		REQUIRE(d.hasNbr(c) == d2.hasNbr(c));
		if (d.hasNbr(c)) {
			CHECK(d.corner_nbr_info[c.getIndex()] != d2.corner_nbr_info[c.getIndex()]);
			switch (d.getNbrType(c)) {
				case NbrType::Normal:
					CHECK(d.getNormalNbrInfo(c).id == d2.getNormalNbrInfo(c).id);
					break;
				case NbrType::Fine:
					CHECK(d.getFineNbrInfo(c).ids[0] == d2.getFineNbrInfo(c).ids[0]);
					break;
				case NbrType::Coarse:
					CHECK(d.getCoarseNbrInfo(c).id == d2.getCoarseNbrInfo(c).id);
					break;
			}
		}
	}
	for (Edge c : Edge::getValues()) {
		REQUIRE(d.hasNbr(c) == d2.hasNbr(c));
		if (d.hasNbr(c)) {
			CHECK(d.edge_nbr_info[c.getIndex()] != d2.edge_nbr_info[c.getIndex()]);
			switch (d.getNbrType(c)) {
				case NbrType::Normal:
					CHECK(d.getNormalNbrInfo(c).id == d2.getNormalNbrInfo(c).id);
					break;
				case NbrType::Fine:
					CHECK(d.getFineNbrInfo(c).ids[0] == d2.getFineNbrInfo(c).ids[0]);
					break;
				case NbrType::Coarse:
					CHECK(d.getCoarseNbrInfo(c).id == d2.getCoarseNbrInfo(c).id);
					break;
			}
		}
	}
}