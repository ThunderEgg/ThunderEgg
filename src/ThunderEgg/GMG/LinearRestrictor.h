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

#ifndef THUNDEREGG_GMG_LINEARRESTRICTOR_H
#define THUNDEREGG_GMG_LINEARRESTRICTOR_H
#include <ThunderEgg/GMG/MPIRestrictor.h>
#include <memory>
namespace ThunderEgg
{
namespace GMG
{
/**
 * @brief Restrictor that averages the corresponding fine cells into each coarse cell.
 */
template <int D> class LinearRestrictor : public MPIRestrictor<D>
{
	private:
	/**
	 * @brief true if ghost values at boundaries should be extrapolated
	 */
	bool extrapolate_boundary_ghosts;

	/**
	 * @brief Extrapolate to the ghosts on the parent patch
	 *
	 * @param pinfo the patch
	 * @param fine_data the finer patch
	 * @param coarse_data the coarser patch
	 */
	void extrapolateBoundaries(const PatchInfo<D> &pinfo, const LocalData<D> &fine_data, LocalData<D> &coarse_data) const
	{
		Orthant<D>         orth = pinfo.orth_on_parent;
		std::array<int, D> starts;
		for (size_t i = 0; i < D; i++) {
			starts[i] = orth.isOnSide(Side<D>::LowerSideOnAxis(i)) ? 0 : coarse_data.getLengths()[i];
		}
		// extrapolate ghost values
		for (Side<D> s : pinfo.orth_on_parent.getExteriorSides()) {
			if (!pinfo.hasNbr(s)) {
				auto fine_ghost    = fine_data.getGhostSliceOnSide(s, 1);
				auto fine_interior = fine_data.getSliceOnSide(s);
				auto coarse_ghost  = coarse_data.getGhostSliceOnSide(s, 1);
				nested_loop<D - 1>(fine_ghost.getStart(), fine_ghost.getEnd(), [&](const std::array<int, D - 1> &coord) {
					std::array<int, D - 1> coarse_coord;
					for (size_t x = 0; x < s.getAxisIndex(); x++) {
						coarse_coord[x] = (coord[x] + starts[x]) / 2;
					}
					for (size_t x = s.getAxisIndex() + 1; x < D; x++) {
						coarse_coord[x - 1] = (coord[x - 1] + starts[x]) / 2;
					}
					coarse_ghost[coarse_coord] += (3 * fine_ghost[coord] - fine_interior[coord]) / (1 << D);
				});
			}
		}
	}
	/**
	 * @brief Restrict to a coarser patch
	 *
	 * @param pinfo the patch
	 * @param parent_index the index of the coarser patch
	 * @param finer_vector the finer vector
	 * @param coarser_vector the coarser vector
	 */
	void restrictToCoarserParent(const PatchInfo<D> &             pinfo,
	                             int                              parent_index,
	                             std::shared_ptr<const Vector<D>> finer_vector,
	                             std::shared_ptr<Vector<D>>       coarser_vector) const
	{
		auto coarse_local_datas = coarser_vector->getLocalDatas(parent_index);
		auto fine_datas         = finer_vector->getLocalDatas(pinfo.local_index);
		// get starting index in coarser patch
		Orthant<D>         orth = pinfo.orth_on_parent;
		std::array<int, D> starts;
		for (size_t i = 0; i < D; i++) {
			starts[i] = orth.isOnSide(Side<D>::LowerSideOnAxis(i)) ? 0 : coarse_local_datas[0].getLengths()[i];
		}

		for (size_t c = 0; c < fine_datas.size(); c++) {
			// interpolate interior values
			nested_loop<D>(fine_datas[c].getStart(), fine_datas[c].getEnd(), [&](const std::array<int, D> &coord) {
				std::array<int, D> coarse_coord;
				for (size_t x = 0; x < D; x++) {
					coarse_coord[x] = (coord[x] + starts[x]) / 2;
				}
				coarse_local_datas[c][coarse_coord] += fine_datas[c][coord] / (1 << D);
			});

			if (extrapolate_boundary_ghosts) {
				extrapolateBoundaries(pinfo, fine_datas[c], coarse_local_datas[c]);
			}
		}
	}
	/**
	 * @brief Copy to a parent patch
	 *
	 * @param pinfo the patch
	 * @param parent_index the index of the coarser patch
	 * @param finer_vector the finer vector
	 * @param coarser_vector the coarser vector
	 */

	void copyToParent(const PatchInfo<D> &             pinfo,
	                  int                              parent_index,
	                  std::shared_ptr<const Vector<D>> finer_vector,
	                  std::shared_ptr<Vector<D>>       coarser_vector) const
	{
		auto coarse_local_datas = coarser_vector->getLocalDatas(parent_index);
		auto fine_datas         = finer_vector->getLocalDatas(pinfo.local_index);
		for (size_t c = 0; c < fine_datas.size(); c++) {
			// just copy the values
			nested_loop<D>(fine_datas[c].getStart(), fine_datas[c].getEnd(), [&](const std::array<int, D> &coord) {
				coarse_local_datas[c][coord] += fine_datas[c][coord];
			});
			if (extrapolate_boundary_ghosts) {
				// copy boundary ghost values
				for (Side<D> s : Side<D>::getValues()) {
					if (!pinfo.hasNbr(s)) {
						auto fine_ghost   = fine_datas[c].getGhostSliceOnSide(s, 1);
						auto coarse_ghost = coarse_local_datas[c].getGhostSliceOnSide(s, 1);
						nested_loop<D - 1>(fine_ghost.getStart(), fine_ghost.getEnd(), [&](const std::array<int, D - 1> &coord) {
							coarse_ghost[coord] += fine_ghost[coord];
						});
					}
				}
			}
		}
	}

	public:
	/**
	 * @brief Create new LinearRestrictor object.
	 *
	 * @param fine_domain the finer Domain
	 * @param coarse_domain the coarser Domain
	 * @param num_components the number of components in each cell
	 * @param extrapolate_boundary_ghosts set to true if ghost values at the boundaries should be
	 * extrapolated
	 */
	LinearRestrictor(std::shared_ptr<Domain<D>> fine_domain,
	                 std::shared_ptr<Domain<D>> coarse_domain,
	                 int                        num_components,
	                 bool                       extrapolate_boundary_ghosts = false)
	: MPIRestrictor<D>(std::make_shared<InterLevelComm<D>>(coarse_domain, num_components, fine_domain)),
	  extrapolate_boundary_ghosts(extrapolate_boundary_ghosts)
	{
	}
	void restrictPatches(const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>> &patches,
	                     std::shared_ptr<const Vector<D>>                                               finer_vector,
	                     std::shared_ptr<Vector<D>>                                                     coarser_vector) const override
	{
		for (const auto &pair : patches) {
			if (pair.second.get().hasCoarseParent()) {
				restrictToCoarserParent(pair.second.get(), pair.first, finer_vector, coarser_vector);
			} else {
				copyToParent(pair.second.get(), pair.first, finer_vector, coarser_vector);
			}
		}
	}
};
} // namespace GMG
} // namespace ThunderEgg
// explicit instantiation
extern template class ThunderEgg::GMG::LinearRestrictor<2>;
extern template class ThunderEgg::GMG::LinearRestrictor<3>;
#endif