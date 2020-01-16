/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019-2020 Thunderegg Developers. See AUTHORS.md file at the
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
#include <Thunderegg/Domain.h>
#include <Thunderegg/GMG/InterLevelComm.h>
#include <Thunderegg/GMG/Level.h>
#include <Thunderegg/GMG/Restrictor.h>
#include <memory>
namespace Thunderegg
{
namespace GMG
{
/**
 * @brief Restrictor that averages the corresponding fine cells into each coarse cell.
 */
template <size_t D> class MPIRestrictor : public Restrictor<D>
{
	private:
	/**
	 * @brief The communication package for restricting between levels.
	 */
	std::shared_ptr<InterLevelComm<D>> ilc;

	/**
	 * @brief Restrict values into coarse vector
	 *
	 * The idea behind this is that this function will be called twice. Once to fill in the ghost
	 * values, and once to fill in the local values. The ghost values will be filled first and the
	 * local values will be fill while MPI communication is happening.
	 *
	 * @param patches pairs where the first value is the index in the coarse vector and the second
	 * value is a pointer to the PatchInfo object
	 * @param finer_vector the finer vector
	 * @param coarser_vector the coaser vector
	 */
	virtual void
	restrictPatches(const std::vector<std::pair<int, std::shared_ptr<const PatchInfo<D>>>> &patches,
	                std::shared_ptr<const Vector<D>> finer_vector,
	                std::shared_ptr<Vector<D>>       coarser_vector) const = 0;

	public:
	/**
	 * @brief Create new LinearRestrictor object.
	 *
	 * @param ilc the communcation package for the two levels.
	 */
	MPIRestrictor(std::shared_ptr<InterLevelComm<D>> ilc_in)
	{
		this->ilc = ilc_in;
	}
	/**
	 * @brief restriction function
	 *
	 * @param coarse the output vector that is restricted to.
	 * @param fine the input vector that is restricted.
	 */
	void restrict(std::shared_ptr<Vector<D>> coarse, std::shared_ptr<const Vector<D>> fine) const
	{
		std::shared_ptr<Vector<D>> coarse_ghost = ilc->getNewGhostVector();

		// fill in ghost values
		restrictPatches(ilc->getPatchesWithGhostParent(), fine, coarse_ghost);

		// clear values in coarse vector
		coarse->set(0);

		// start scatter for ghost values
		ilc->sendGhostPatchesStart(coarse, coarse_ghost);

		// fill in local values
		restrictPatches(ilc->getPatchesWithLocalParent(), fine, coarse);

		// finish scatter for ghost values
		ilc->sendGhostPatchesFinish(coarse, coarse_ghost);
	}
};
} // namespace GMG
} // namespace Thunderegg
// explicit instantiation
extern template class Thunderegg::GMG::MPIRestrictor<2>;
extern template class Thunderegg::GMG::MPIRestrictor<3>;
#endif