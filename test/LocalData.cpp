#include "catch.hpp"
#include <ThunderEgg/Vector.h>
using namespace std;
using namespace ThunderEgg;
TEST_CASE("LocalData constructor", "[DomainTools]")
{
	auto nx        = GENERATE(1, 2, 10, 13);
	auto ny        = GENERATE(1, 2, 10, 13);
	auto num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	array<int, 2> lengths = {nx, ny};
	array<int, 2> strides = {1, nx + 2 * num_ghost};
	int           size    = 1;
	for (size_t i = 0; i < 2; i++) {
		size *= (lengths[i] + 2 * num_ghost);
	}
	vector<double> vec(size);
	iota(vec.begin(), vec.end(), 0);

	LocalData<2> ld(vec.data() + num_ghost * strides[0] + num_ghost * strides[1], strides, lengths,
	                num_ghost);

	CHECK(ld.getNumGhostCells() == num_ghost);
	CHECK(ld.getPtr() == vec.data() + num_ghost * strides[0] + num_ghost * strides[1]);
	for (int i = 0; i < 2; i++) {
		CHECK(ld.getLengths()[i] == lengths[i]);
		CHECK(ld.getStrides()[i] == strides[i]);
		CHECK(ld.getStart()[i] == 0);
		CHECK(ld.getEnd()[i] == lengths[i] - 1);
		CHECK(ld.getGhostStart()[i] == -num_ghost);
		CHECK(ld.getGhostEnd()[i] == lengths[i] - 1 + num_ghost);
	}
	CHECK(ld.getPtr(ld.getGhostStart()) == vec.data());
}