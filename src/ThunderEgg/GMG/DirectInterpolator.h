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

#ifndef THUNDEREGG_GMG_DRCTINTP_H
#define THUNDEREGG_GMG_DRCTINTP_H

#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GMG/InterLevelComm.h>
#include <ThunderEgg/GMG/MPIInterpolator.h>
#include <memory>

namespace ThunderEgg
{
namespace GMG
{
/**
 * @brief Simple class that directly places values from coarse cell into the corresponding fine
 * cells.
 */
template <int D> class DirectInterpolator : public MPIInterpolator<D>
{
	public:
	/**
	 * @brief Create new DirectInterpolator object.
	 *
	 * @param coarse_domain the coarser Domain
	 * @param fine_domain the finer Domain
	 * @param num_components the number of components in each cell
	 */
	DirectInterpolator(std::shared_ptr<Domain<D>> coarse_domain,
	                   std::shared_ptr<Domain<D>> fine_domain, int num_components)
	: MPIInterpolator<D>(
	  std::make_shared<InterLevelComm<D>>(coarse_domain, num_components, fine_domain))
	{
	}
	void interpolatePatches(
	const std::vector<std::pair<int, std::shared_ptr<const PatchInfo<D>>>> &patches,
	std::shared_ptr<const Vector<D>>                                        coarser_vector,
	std::shared_ptr<Vector<D>> finer_vector) const override
	{
		for (auto pair : patches) {
			auto pinfo              = pair.second;
			auto coarse_local_datas = coarser_vector->getLocalDatas(pair.first);
			auto fine_datas         = finer_vector->getLocalDatas(pinfo->local_index);

			if (pinfo->hasCoarseParent()) {
				Orthant<D>         orth = pinfo->orth_on_parent;
				std::array<int, D> starts;
				for (size_t i = 0; i < D; i++) {
					starts[i]
					= orth.isOnSide(Side<D>(2 * i)) ? 0 : coarse_local_datas[0].getLengths()[i];
				}

				for (size_t c = 0; c < fine_datas.size(); c++) {
					nested_loop<D>(fine_datas[c].getStart(), fine_datas[c].getEnd(),
					               [&](const std::array<int, D> &coord) {
						               std::array<int, D> coarse_coord;
						               for (size_t x = 0; x < D; x++) {
							               coarse_coord[x] = (coord[x] + starts[x]) / 2;
						               }
						               fine_datas[c][coord] += coarse_local_datas[c][coarse_coord];
					               });
				}
			} else {
				for (size_t c = 0; c < fine_datas.size(); c++) {
					nested_loop<D>(fine_datas[c].getStart(), fine_datas[c].getEnd(),
					               [&](const std::array<int, D> &coord) {
						               fine_datas[c][coord] += coarse_local_datas[c][coord];
					               });
				}
			}
		}
	}
};
extern template class DirectInterpolator<2>;
extern template class DirectInterpolator<3>;
} // namespace GMG
} // namespace ThunderEgg
#endif