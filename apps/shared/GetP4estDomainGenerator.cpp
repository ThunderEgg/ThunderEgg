/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2021      Scott Aiton
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
#include "GetP4estDomainGenerator.h"
#include <ThunderEgg/tpl/json.hpp>
#include <fstream>
#include <p4est_connectivity.h>
#include <p4est_extended.h>
using namespace ThunderEgg::tpl;
int refine_fn(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
	const nlohmann::json &json = *(nlohmann::json *) p4est->user_pointer;
	// double x = p4est_root_length
	nlohmann::json node = json["forest"]["nodes"][std::to_string(0)];
	for (int i = 1; i <= quadrant->level; i++) {
		if (node["child_ids"][0] != -1) {
			int quad = (quadrant->x / P4EST_QUADRANT_LEN(i)) & 0x1;
			quad |= ((quadrant->y / P4EST_QUADRANT_LEN(i)) & 0x1) << 1;
			node = json["forest"]["nodes"][std::to_string(node["child_ids"][quad].get<int>())];
		} else {
			break;
		}
	}
	return node["child_ids"][0] != -1;
}

ThunderEgg::P4estDomainGenerator
getP4estDomainGenerator(p4est_connectivity *conn, const std::string &mesh_file, std::array<int, 2> ns,
                        int                                                   num_ghost_cells,
                        const ThunderEgg::P4estDomainGenerator::BlockMapFunc &bmf)
{
	nlohmann::json json;
	std::ifstream  input_stream(mesh_file);
	if (!input_stream.good()) {
		throw "could not open file";
	}
	input_stream >> json;
	input_stream.close();

	p4est *p4est = p4est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, &json);
	p4est_refine(p4est, true, refine_fn, nullptr);
	p4est_partition(p4est, true, nullptr);
	ThunderEgg::P4estDomainGenerator dom_gen(p4est, ns, num_ghost_cells, bmf);
	p4est_destroy(p4est);
	return dom_gen;
}