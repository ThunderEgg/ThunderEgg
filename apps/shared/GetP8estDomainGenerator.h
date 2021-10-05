#include <ThunderEgg/P8estDomainGenerator.h>
ThunderEgg::P8estDomainGenerator
getP8estDomainGenerator(p8est_connectivity_t *conn, const std::string &mesh_file, std::array<int, 3> ns,
                        int                                                   num_ghost_cells,
                        const ThunderEgg::P8estDomainGenerator::BlockMapFunc &bmf);
