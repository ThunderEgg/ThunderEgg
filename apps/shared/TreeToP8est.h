#ifndef TREETOP8EST
#define TREETOP8EST
#include <ThunderEgg/Experimental/OctTree.h>
#include <p8est.h>
#include <p8est_mesh.h>
class TreeToP8est
{
	public:
	p8est_connectivity_t *conn;
	p8est_t *             p8est;
	p8est_ghost_t *       ghost;
	p8est_mesh_t *        mesh;
	TreeToP8est(ThunderEgg::Experimental::Tree<3> t);
};
#endif
