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

#ifndef THUNDEREGG_GMG_LINEARRESTRICTOR_H
#define THUNDEREGG_GMG_LINEARRESTRICTOR_H
#include <Thunderegg/GMG/MPIRestrictor.h>
#include <memory>
namespace Thunderegg
{
namespace GMG
{
/**
 * @brief Restrictor that averages the corresponding fine cells into each coarse cell.
 */
template <size_t D> class LinearRestrictor : public MPIRestrictor<D>
{
	private:
	void
	restrictPatches(const std::vector<std::pair<int, std::shared_ptr<const PatchInfo<D>>>> &patches,
	                std::shared_ptr<const Vector<D>> finer_vector,
	                std::shared_ptr<Vector<D>>       coarser_vector) const override
	{
		for (const auto &pair : patches) {
			std::shared_ptr<const PatchInfo<D>> pinfo = pair.second;
			LocalData<D> coarse_local_data            = coarser_vector->getLocalData(pair.first);
			LocalData<D> fine_data = finer_vector->getLocalData(pinfo->local_index);

			if (pinfo->hasCoarseParent()) {
				Orthant<D> orth = pinfo->orth_on_parent;
				// get starting index in coarser patch
				std::array<int, D> starts;
				for (size_t i = 0; i < D; i++) {
					starts[i] = orth.isOnSide(2 * i) ? 0 : coarse_local_data.getLengths()[i];
				}

				// interpolate interior values
				nested_loop<D>(fine_data.getStart(), fine_data.getEnd(),
				               [&](const std::array<int, D> &coord) {
					               std::array<int, D> coarse_coord;
					               for (size_t x = 0; x < D; x++) {
						               coarse_coord[x] = (coord[x] + starts[x]) / 2;
					               }
					               coarse_local_data[coarse_coord] += fine_data[coord] / (1 << D);
				               });

				// extrapolate ghost values
				for (Side<D> s : orth.getExteriorSides()) {
					auto fine_ghost    = fine_data.getGhostSliceOnSide(s, 1);
					auto fine_interior = fine_data.getSliceOnSide(s);
					auto coarse_ghost  = coarse_local_data.getGhostSliceOnSide(s, 1);
					nested_loop<D>(fine_ghost.getStart(), fine_ghost.getEnd(),
					               [&](const std::array<int, D - 1> &coord) {
						               std::array<int, D - 1> coarse_coord;
						               for (size_t x = 0; x < s.axis(); x++) {
							               coarse_coord[x] = (coord[x] + starts[x]) / 2;
						               }
						               for (size_t x = s.axis() + 1; x < D; x++) {
							               coarse_coord[x - 1] = (coord[x - 1] + starts[x]) / 2;
						               }
						               coarse_ghost[coarse_coord]
						               += (3 * fine_ghost[coord] - fine_interior[coord]) / (1 << D);
					               });
				}
			} else {
				// just copy the values
				nested_loop<D>(fine_data.getGhostStart(), fine_data.getGhostEnd(),
				               [&](const std::array<int, D> &coord) {
					               coarse_local_data[coord] += fine_data[coord];
				               });
			}
		}
	}

	public:
	/**
	 * @brief Create new LinearRestrictor object.
	 *
	 * @param ilc the communcation package for the two levels.
	 */
	LinearRestrictor(std::shared_ptr<InterLevelComm<D>> ilc_in) : MPIRestrictor<D>(ilc_in) {}
	/**
	 * @brief Generates LinearRestrictors for GMG levels
	 */
	class Generator
	{
		public:
		/**
		 * @brief Generate a LinearRestrictor from the finer level to the coarser level
		 *
		 * @param level the finer level
		 * @return std::shared_ptr<const LinearRestrictor<D>> the restrictor
		 */
		std::shared_ptr<const LinearRestrictor<D>> operator()(std::shared_ptr<const Level<D>> level)
		{
			auto ilc = std::make_shared<InterLevelComm<D>>(level->getCoarser()->getDomain(),
			                                               level->getDomain());
			return std::make_shared<LinearRestrictor<D>>(ilc);
		}
	};
};
} // namespace GMG
} // namespace Thunderegg
// explicit instantiation
extern template class Thunderegg::GMG::LinearRestrictor<2>;
extern template class Thunderegg::GMG::LinearRestrictor<3>;
#endif