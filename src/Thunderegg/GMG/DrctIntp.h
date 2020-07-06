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

#ifndef THUNDEREGG_GMG_DRCTINTP_H
#define THUNDEREGG_GMG_DRCTINTP_H

#include <Thunderegg/Domain.h>
#include <Thunderegg/GMG/InterLevelComm.h>
#include <Thunderegg/GMG/MPIInterpolator.h>
#include <memory>

namespace Thunderegg
{
namespace GMG
{
/**
 * @brief Simple class that directly places values from coarse cell into the corresponding fine
 * cells.
 */
template <size_t D> class DrctIntp : public MPIInterpolator<D>
{
	private:
	void interpolatePatches(
	const std::vector<std::pair<int, std::shared_ptr<const PatchInfo<D>>>> &patches,
	std::shared_ptr<const Vector<D>>                                        coarser_vector,
	std::shared_ptr<Vector<D>> finer_vector) const override
	{
		for (auto pair : patches) {
			auto         pinfo             = pair.second;
			LocalData<D> coarse_local_data = coarser_vector->getLocalData(pair.first);
			LocalData<D> fine_data         = finer_vector->getLocalData(pinfo->local_index);

			if (pinfo->hasCoarseParent()) {
				Orthant<D>         orth = pinfo->orth_on_parent;
				std::array<int, D> starts;
				for (size_t i = 0; i < D; i++) {
					starts[i]
					= orth.isOnSide(Side<D>(2 * i)) ? 0 : coarse_local_data.getLengths()[i];
				}

				nested_loop<D>(fine_data.getStart(), fine_data.getEnd(),
				               [&](const std::array<int, D> &coord) {
					               std::array<int, D> coarse_coord;
					               for (size_t x = 0; x < D; x++) {
						               coarse_coord[x] = (coord[x] + starts[x]) / 2;
					               }
					               fine_data[coord] += coarse_local_data[coarse_coord];
				               });
			} else {
				nested_loop<D>(fine_data.getStart(), fine_data.getEnd(),
				               [&](const std::array<int, D> &coord) {
					               fine_data[coord] += coarse_local_data[coord];
				               });
			}
		}
	}

	public:
	/**
	 * @brief Create new DrctIntp object.
	 *
	 * @param ilc the communcation package for the two levels.
	 */
	DrctIntp(std::shared_ptr<InterLevelComm<D>> ilc_in) : MPIInterpolator<D>(ilc_in) {}

	class Generator
	{
		public:
		Generator() {}
		std::shared_ptr<const DrctIntp<D>> operator()(std::shared_ptr<const Level<D>> level)
		{
			auto ilc = std::make_shared<InterLevelComm<D>>(level->getDomain(),
			                                               level->getFiner()->getDomain());
			return std::make_shared<DrctIntp<D>>(ilc);
		}
	};
};

} // namespace GMG
} // namespace Thunderegg
#endif