#include <ThunderEgg/P4estDomainGenerator.h>
ThunderEgg::P4estDomainGenerator
getP4estDomainGenerator(p4est_connectivity_t *conn, const std::string &mesh_file, std::array<int, 2> ns,
                        int                                                   num_ghost_cells,
                        const ThunderEgg::P4estDomainGenerator::BlockMapFunc &bmf);
