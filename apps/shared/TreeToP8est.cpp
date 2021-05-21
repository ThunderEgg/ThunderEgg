#include "TreeToP8est.h"
#include <p8est_connectivity.h>
#include <p8est_extended.h>
using namespace ThunderEgg::Experimental;
int refine_fn(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *quadrant)
{
	Tree<3> &t = *(Tree<3> *) p8est->user_pointer;
	// double x = p4est_root_length
	Node<3> node = t.nodes[t.root];
	for (int i = 1; i <= quadrant->level; i++) {
		if (node.hasChildren()) {
			int quad = (quadrant->x / P8EST_QUADRANT_LEN(i)) & 0x1;
			quad |= ((quadrant->y / P8EST_QUADRANT_LEN(i)) & 0x1) << 1;
			quad |= ((quadrant->z / P8EST_QUADRANT_LEN(i)) & 0x1) << 2;
			node = t.nodes[node.child_id[quad]];
		} else {
			break;
		}
	}
	return node.hasChildren();
}
TreeToP8est::TreeToP8est(Tree<3> t)
{
	conn  = p8est_connectivity_new_unitcube();
	p8est = p8est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, &t);
	p8est_refine(p8est, true, refine_fn, nullptr);
	p8est_partition(p8est, true, nullptr);
	ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
	mesh  = p8est_mesh_new_ext(p8est, ghost, 1, 1, P8EST_CONNECT_FULL);
}
