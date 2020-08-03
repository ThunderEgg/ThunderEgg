/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#include "catch.hpp"
#include <ThunderEgg/Domain.h>
using namespace std;
using namespace ThunderEgg;
TEST_CASE("Domain construisctors work", "[Domain]")
{
	map<int, shared_ptr<PatchInfo<2>>> pinfo_map;

	auto n         = GENERATE(1, 2, 10, 13);
	auto spacing   = GENERATE(0.01, 1.0, 3.14);
	auto num_ghost = GENERATE(0, 1, 2, 3, 4, 5);

	pinfo_map[0].reset(new PatchInfo<2>());
	pinfo_map[0]->id = 0;
	pinfo_map[0]->ns.fill(n);
	pinfo_map[0]->spacings.fill(spacing);
	pinfo_map[0]->num_ghost_cells = num_ghost;
	Domain<2> d(pinfo_map, {n, n}, num_ghost);

	// check getters
	for (int ni : d.getNs()) {
		CHECK(ni == n);
	}
	CHECK(d.getNumGlobalPatches() == 1);
	CHECK(d.getNumLocalPatches() == 1);
	CHECK(d.getNumGlobalCells() == n * n);
	CHECK(d.getNumLocalCells() == n * n);
	CHECK(d.getNumLocalCellsWithGhost() == (n + 2 * num_ghost) * (n + 2 * num_ghost));
	CHECK(d.getNumLocalBCCells() == 4 * n);
	CHECK(d.getNumCellsInPatch() == n * n);
	CHECK(d.getNumGhostCells() == num_ghost);
	CHECK(d.volume() == Approx(spacing * spacing * n * n));
	// TODO Check intigrate
}