/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
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

#ifndef THUNDEREGG_GMG_INTERLEVELCOMM_H
#define THUNDEREGG_GMG_INTERLEVELCOMM_H

#include <Thunderegg/Domain.h>

namespace Thunderegg
{
namespace GMG
{
/**
 * @brief Exception that the InterLevelComm class trows
 */
struct InterLevelCommException : std::runtime_error {
	InterLevelCommException(std::string message) : std::runtime_error(message){};
};
/**
 * @brief Creates a mapping from fine to coarse levels.
 */
template <size_t D> class InterLevelComm
{
	private:
	public:
	/**
	 * @brief Create a new InterLevelComm object.
	 *
	 * @param coarse_domain the coarser DomainCollection.
	 * @param fine_domain the finer DomainCollection.
	 */
	InterLevelComm(std::shared_ptr<const Domain<D>> coarser_domain,
	               std::shared_ptr<const Domain<D>> finer_domain)
	{
	}

	/**
	 * @brief Allocate a new vector for ghost patch values
	 *
	 * @return the newly allocated vector.
	 */
	std::shared_ptr<Vector<D>> getNewGhostVector() const {}

	const std::vector<std::pair<int, std::shared_ptr<const PatchInfo<D>>>> &
	getPatchesWithLocalParent() const
	{
	}
	const std::vector<std::pair<int, std::shared_ptr<const PatchInfo<D>>>> &
	getPatchesWithGhostParent() const
	{
	}

	/**
	 * @brief Get a PETSc VecScatter object that will scatter values from the coarse vector to the
	 * processes that need them for interpolation / restriction.
	 *
	 * @return the VecScatter object
	 */
	void sendGhostPatchesStart(std::shared_ptr<Vector<D>>       vector,
	                           std::shared_ptr<const Vector<D>> ghost_vector)
	{
	}
	void sendGhostPatchesFinish(std::shared_ptr<Vector<D>>       vector,
	                            std::shared_ptr<const Vector<D>> ghost_vector)
	{
	}
	void getGhostPatchesStart(std::shared_ptr<const Vector<D>> vector,
	                          std::shared_ptr<Vector<D>>       ghost_vector)
	{
	}
	void getGhostPatchesFinish(std::shared_ptr<const Vector<D>> vector,
	                           std::shared_ptr<Vector<D>>       ghost_vector)
	{
	}
}; // namespace GMG
} // namespace GMG
} // namespace Thunderegg
#endif