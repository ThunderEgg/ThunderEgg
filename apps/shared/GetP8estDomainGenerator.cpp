#include "GetP8estDomainGenerator.h"
#include <ThunderEgg/tpl/json.hpp>
#include <fstream>
#include <p8est_connectivity.h>
#include <p8est_extended.h>
int refine_fn(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *quadrant)
{
	const ThunderEgg::tpl::nlohmann::json &json = *(ThunderEgg::tpl::nlohmann::json *) p8est->user_pointer;
	// double x = p4est_root_length
	ThunderEgg::tpl::nlohmann::json node = json["forest"]["nodes"][std::to_string(0)];
	for (int i = 1; i <= quadrant->level; i++) {
		if (node["child_ids"][0] != -1) {
			int quad = (quadrant->x / P8EST_QUADRANT_LEN(i)) & 0x1;
			quad |= ((quadrant->y / P8EST_QUADRANT_LEN(i)) & 0x1) << 1;
			quad |= ((quadrant->z / P8EST_QUADRANT_LEN(i)) & 0x1) << 2;
			node = json["forest"]["nodes"][std::to_string(node["child_ids"][quad].get<int>())];
		} else {
			break;
		}
	}
	return node["child_ids"][0] != -1;
}

ThunderEgg::P8estDomainGenerator
getP8estDomainGenerator(p8est_connectivity *conn, const std::string &mesh_file, std::array<int, 3> ns,
                        int                                                   num_ghost_cells,
                        const ThunderEgg::P8estDomainGenerator::BlockMapFunc &bmf)
{
	ThunderEgg::tpl::nlohmann::json json;
	std::ifstream                   input_stream(mesh_file);
	if (!input_stream.good()) {
		throw "could not open file";
	}
	input_stream >> json;
	input_stream.close();

	p8est *p8est = p8est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, &json);
	p8est_refine(p8est, true, refine_fn, nullptr);
	p8est_partition(p8est, true, nullptr);
	ThunderEgg::P8estDomainGenerator dom_gen(p8est, ns, num_ghost_cells, bmf);
	p8est_destroy(p8est);
	return dom_gen;
}