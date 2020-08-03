/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019-2020 ThunderEgg Developers. See AUTHORS.md file at the
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

#ifndef THUNDEREGG_GMG_MPIINTERPOLATOR_H
#define THUNDEREGG_GMG_MPIINTERPOLATOR_H
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GMG/InterLevelComm.h>
#include <ThunderEgg/GMG/Interpolator.h>
#include <ThunderEgg/GMG/Level.h>
#include <memory>
namespace ThunderEgg
{
namespace GMG
{
/**
 * @brief Interpolator that implements the necessary mpi calls, derived classes only have to
 * implement interpolatePatches method
 */
template <size_t D> class MPIInterpolator : public Interpolator<D>
{
	private:
	/**
	 * @brief The communication package for restricting between levels.
	 */
	std::shared_ptr<InterLevelComm<D>> ilc;

	/**
	 * @brief Interpolate values from coarse vector to the finer vector
	 *
	 * The idea behind this is that this function will be called twice. Once to interpolate from the
	 * local values, and once to interpolate from the ghost values. The local values will be
	 * interpolated from first, while MPI communication is happening, and the ghost values will be
	 * interpolated from last.
	 *
	 * @param patches pairs where the first value is the index in the coarse vector and the second
	 * value is a pointer to the PatchInfo object
	 * @param finer_vector the finer vector
	 * @param coarser_vector the coaser vector
	 */
	virtual void interpolatePatches(
	const std::vector<std::pair<int, std::shared_ptr<const PatchInfo<D>>>> &patches,
	std::shared_ptr<const Vector<D>>                                        coarser_vector,
	std::shared_ptr<Vector<D>>                                              finer_vector) const = 0;

	public:
	/**
	 * @brief Create new MPIInterpolator object.
	 *
	 * @param ilc the communcation package for the two levels.
	 */
	MPIInterpolator(std::shared_ptr<InterLevelComm<D>> ilc_in)
	{
		this->ilc = ilc_in;
	}
	/**
	 * @brief interpolation function
	 *
	 * @param coarse the input vector that is interpolated from
	 * @param fine the output vector that is interpolated to.
	 */
	void interpolate(std::shared_ptr<const Vector<D>> coarse, std::shared_ptr<Vector<D>> fine) const
	{
		std::shared_ptr<Vector<D>> coarse_ghost = ilc->getNewGhostVector();

		// start scatter for ghost values
		ilc->getGhostPatchesStart(coarse, coarse_ghost);

		// interpolate form local values
		interpolatePatches(ilc->getPatchesWithLocalParent(), coarse, fine);

		// finish scatter for ghost values
		ilc->getGhostPatchesFinish(coarse, coarse_ghost);

		// interpolator from ghost values
		interpolatePatches(ilc->getPatchesWithGhostParent(), coarse_ghost, fine);
	}
};
} // namespace GMG
} // namespace ThunderEgg
// explicit instantiation
extern template class ThunderEgg::GMG::MPIInterpolator<2>;
extern template class ThunderEgg::GMG::MPIInterpolator<3>;
#endif