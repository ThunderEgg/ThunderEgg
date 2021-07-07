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

#ifndef THUNDEREGG_P4ESTDOMAINGENERATOR_H
#define THUNDEREGG_P4ESTDOMAINGENERATOR_H
#include <ThunderEgg/DomainGenerator.h>
#include <functional>
#include <list>
#include <p4est_extended.h>
namespace ThunderEgg
{
/**
 * @brief Generates Domain objects form a given p4est object
 */
class P4estDomainGenerator : public DomainGenerator<2>
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
	 * @param x the resulting x coordinate of the mapping function
	 * @param y the resulting y coordinate of the mapping function
	 */
	using BlockMapFunc = std::function<void(int block_no, double unit_x, double unit_y, double &x, double &y)>;

	private:
	/**
	 * @brief copy of p4est tree
	 */
	p4est_t *my_p4est;
	/**
	 * @brief List of the domain patches
	 *
	 * Finesr domain is stored in back
	 */
	std::list<std::vector<PatchInfo<2>>> domain_patches;

	/**
	 * @brief The dimensions of each patch
	 */
	std::array<int, 2> ns;
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
	 * @brief Construct a new P4estDomainGenerator object
	 *
	 * @param p4est the p4est object
	 * @param ns the number of cells in each direction
	 * @param num_ghost_cells the number of ghost cells on each side of the patch
	 * @param bmf the function used to map the blocks to the domain
	 */
	P4estDomainGenerator(p4est_t *p4est, const std::array<int, 2> &ns, int num_ghost_cells, const BlockMapFunc &bmf);
	~P4estDomainGenerator();
	P4estDomainGenerator(const P4estDomainGenerator &) = delete;
	P4estDomainGenerator(P4estDomainGenerator &&)      = delete;
	P4estDomainGenerator &operator=(const P4estDomainGenerator &) = delete;
	P4estDomainGenerator &operator=(P4estDomainGenerator &&) = delete;
	Domain<2>             getFinestDomain();
	bool                  hasCoarserDomain();
	Domain<2>             getCoarserDomain();
};
} // namespace ThunderEgg
#endif
