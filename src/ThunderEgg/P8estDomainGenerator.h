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

#ifndef THUNDEREGG_P8ESTDOMAINGENERATOR_H
#define THUNDEREGG_P8ESTDOMAINGENERATOR_H
/**
 * @file
 *
 * @brief P8estDomainGenerator class
 */
#include <ThunderEgg/DomainGenerator.h>
#include <functional>
#include <p8est.h>
namespace ThunderEgg
{
/**
 * @brief Generates Domain objects form a given p4est object
 */
class P8estDomainGenerator : public DomainGenerator<3>
{
	public:
	/**
	 * @brief Maps coordinate in a block to a coordinate in the domain.
	 *
	 * Each block is treated as a unit square. The input wil be the block number and a coordinate
	 * withing that unit square.
	 *
	 * @param block_no the block number
	 * @param unit_x the x coordinate in the block
	 * @param unit_y the y coordinate in the block
	 * @param unit_z the z coordinate in the block
	 * @param x the resulting x coordinate of the mapping function
	 * @param y the resulting y coordinate of the mapping function
	 * @param z the resulting z coordinate of the mapping function
	 */
	using BlockMapFunc = std::function<void(int block_no, double unit_x, double unit_y, double unit_z, double &x, double &y, double &z)>;

	private:
	/**
	 * @brief copy of p8est tree
	 */
	p8est_t *my_p8est;
	/**
	 * @brief List of the domain patches
	 *
	 * Finesr domain is stored in back
	 */
	std::list<std::vector<PatchInfo<3>>> domain_patches;

	/**
	 * @brief The dimensions of each patch
	 */
	std::array<int, 3> ns;
	/**
	 * @brief the number of ghost cells on each side of the patch
	 */
	int num_ghost_cells;
	/**
	 * @brief The current level that has been generated.
	 *
	 * Will start with num_levels-1
	 */
	int curr_level;
	/**
	 * @brief The block Mapping function being used.
	 */
	BlockMapFunc bmf;

	/**
	 * @brief The id to set the next domain with
	 */
	int id = 0;

	/**
	 * @brief Get a new coarser level and add it to the end of domain_list
	 */
	void extractLevel();

	/**
	 * @brief coaren the p8est tree
	 */
	void coarsenTree();

	/**
	 * @brief create the patchinfo objects using the cureent leaves of the tree
	 */
	void createPatchInfos();

	/**
	 * @brief Add NbrInfo objects to PatchInfo objects
	 */
	void linkNeighbors();

	/**
	 * @brief update the parent_rank for the previous domain
	 */
	void updateParentRanksOfPreviousDomain();

	public:
	/**
	 * @brief Construct a new P8estDomainGenerator object
	 *
	 * @param p8est the p8est object
	 * @param ns the number of cells in each direction
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 * @param bmf the function used to map the blocks to the domain
	 */
	P8estDomainGenerator(p8est_t *p8est, const std::array<int, 3> &ns, int num_ghost_cells, const BlockMapFunc &bmf);
	P8estDomainGenerator(P8estDomainGenerator &);
	P8estDomainGenerator(P8estDomainGenerator &&) = default;
	P8estDomainGenerator &operator                =(const P8estDomainGenerator &);
	P8estDomainGenerator &operator=(P8estDomainGenerator &&) = default;
	~P8estDomainGenerator();
	Domain<3> getFinestDomain();
	bool      hasCoarserDomain();
	Domain<3> getCoarserDomain();
};
} // namespace ThunderEgg
#endif
