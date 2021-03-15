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

#ifndef THUNDEREGG_GMG_CYCLEOPTS_H
#define THUNDEREGG_GMG_CYCLEOPTS_H
#include <string>
namespace ThunderEgg
{
namespace GMG
{
struct CycleOpts {
	/**
	 * @brief The max number of levels in GMG cycle. 0 means no limit.
	 */
	int max_levels = 0;
	/**
	 * @brief Lowest level is guaranteed to have at least this number of patches per processor.
	 */
	double patches_per_proc = 0;
	/**
	 * @brief Number of sweeps on down cycle
	 */
	int pre_sweeps = 1;
	/**
	 * @brief Number of sweeps on up cycle
	 */
	int post_sweeps = 1;
	/**
	 * @brief Number of sweeps inbetween up and down
	 */
	int mid_sweeps = 1;
	/**
	 * @brief Number of sweeps on coarse level
	 */
	int coarse_sweeps = 1;
	/**
	 * @brief Cycle type
	 */
	std::string cycle_type = "V";
};
} // namespace GMG
} // namespace ThunderEgg
#endif