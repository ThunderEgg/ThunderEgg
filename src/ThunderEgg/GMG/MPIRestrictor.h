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

#ifndef THUNDEREGG_GMG_MPIRESTRICTOR_H
#define THUNDEREGG_GMG_MPIRESTRICTOR_H
/**
 * @file
 *
 * @brief MPIRestrictor class
 */
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GMG/InterLevelComm.h>
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/GMG/Restrictor.h>
namespace ThunderEgg::GMG
{
/**
 * @brief Base class that makes the necessary mpi calls, derived classes only have to
 * implement restrictPatches() method
 */
template <int D> class MPIRestrictor : public Restrictor<D>
{
	private:
	/**
	 * @brief The communication package for restricting between levels.
	 */
	mutable InterLevelComm<D> ilc;

	public:
	/**
	 * @brief Create new LinearRestrictor object.
	 *
	 * @param ilc the communcation package for the two levels.
	 */
	MPIRestrictor(const Domain<D> &coarser_domain, const Domain<D> &finer_domain) : ilc(coarser_domain, finer_domain) {}
	Vector<D> restrict(const Vector<D> &fine) const override
	{
		if constexpr (ENABLE_DEBUG) {
			if (fine.getNumLocalPatches() != ilc.getFinerDomain().getNumLocalPatches()) {
				throw RuntimeError("fine vector is incorrect length. Expected Length of " + std::to_string(ilc.getFinerDomain().getNumLocalPatches())
				                   + " but vector was length " + std::to_string(fine.getNumLocalPatches()));
			}
		}
		Vector<D> coarse(ilc.getCoarserDomain(), fine.getNumComponents());
		Vector<D> coarse_ghost = ilc.getNewGhostVector(fine.getNumComponents());

		// fill in ghost values
		restrictPatches(ilc.getPatchesWithGhostParent(), fine, coarse_ghost);

		// clear values in coarse vector
		coarse.setWithGhost(0);

		// start scatter for ghost values
		ilc.sendGhostPatchesStart(coarse, coarse_ghost);

		// fill in local values
		restrictPatches(ilc.getPatchesWithLocalParent(), fine, coarse);

		// finish scatter for ghost values
		ilc.sendGhostPatchesFinish(coarse, coarse_ghost);

		return coarse;
	}
	/**
	 * @brief Restrict values into coarse vector
	 *
	 * The idea behind this is that this function will be called twice. Once to fill in the ghost
	 * values, and once to fill in the local values. The ghost values will be filled first and the
	 * local values will be fill while MPI communication is happening.
	 *
	 * @param patches pairs where the first value is the index in the coarse vector and the second
	 * value is a reference to the PatchInfo object
	 * @param finer_vector the finer vector
	 * @param coarser_vector the coarser vector
	 */
	virtual void restrictPatches(const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>> &patches,
	                             const Vector<D> &                                                              finer_vector,
	                             Vector<D> &                                                                    coarser_vector) const = 0;
};
} // namespace ThunderEgg::GMG
// explicit instantiation
extern template class ThunderEgg::GMG::MPIRestrictor<2>;
extern template class ThunderEgg::GMG::MPIRestrictor<3>;
#endif