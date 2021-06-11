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
	DirectInterpolator(std::shared_ptr<Domain<D>> coarse_domain, std::shared_ptr<Domain<D>> fine_domain, int num_components)
	: MPIInterpolator<D>(std::make_shared<InterLevelComm<D>>(coarse_domain, num_components, fine_domain))
	{
	}
	void interpolatePatches(const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>> &patches,
	                        std::shared_ptr<const Vector<D>>                                               coarser_vector,
	                        std::shared_ptr<Vector<D>>                                                     finer_vector) const override
	{
		for (auto pair : patches) {
			const PatchInfo<D> &       pinfo       = pair.second.get();
			PatchView<const double, D> coarse_view = coarser_vector->getPatchView(pair.first);
			PatchView<double, D>       fine_view   = finer_vector->getPatchView(pinfo.local_index);

			if (pinfo.hasCoarseParent()) {
				Orthant<D>         orth = pinfo.orth_on_parent;
				std::array<int, D> starts;
				for (size_t i = 0; i < D; i++) {
					starts[i] = orth.isOnSide(Side<D>(2 * i)) ? 0 : coarse_view.getEnd()[i] + 1;
				}

				loop_over_interior_indexes<D + 1>(fine_view, [&](const std::array<int, D + 1> &coord) {
					std::array<int, D + 1> coarse_coord;
					for (size_t x = 0; x < D; x++) {
						coarse_coord[x] = (coord[x] + starts[x]) / 2;
					}
					coarse_coord[D] = coord[D];
					fine_view[coord] += coarse_view[coarse_coord];
				});
			} else {
				loop_over_interior_indexes<D + 1>(fine_view, [&](const std::array<int, D + 1> &coord) { fine_view[coord] += coarse_view[coord]; });
			}
		}
	}
};
extern template class DirectInterpolator<2>;
extern template class DirectInterpolator<3>;
} // namespace GMG
} // namespace ThunderEgg
#endif