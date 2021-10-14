/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
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
/**
 * @file
 *
 * @brief MPIInterpolator class
 */
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GMG/InterLevelComm.h>
#include <ThunderEgg/GMG/Interpolator.h>
#include <ThunderEgg/GMG/Level.h>
namespace ThunderEgg::GMG {
/**
 * @brief Base class that makes the necessary mpi calls, derived classes only have to
 * implement interpolatePatches() method
 */
template<int D>
class MPIInterpolator : public Interpolator<D>
{
private:
  /**
   * @brief The communication package for restricting between levels.
   */
  mutable InterLevelComm<D> ilc;

public:
  /**
   * @brief Create new MPIInterpolator object.
   *
   * @param ilc the communcation package for the two levels.
   */
  MPIInterpolator(const Domain<D>& coarser_domain, const Domain<D>& finer_domain)
    : ilc(coarser_domain, finer_domain)
  {}
  /**
   * @brief Interpolate values from coarse vector to the finer vector
   *
   * The idea behind this is that this function will be called twice. Once to interpolate from the
   * local values, and once to interpolate from the ghost values. The local values will be
   * interpolated from first, while MPI communication is happening, and the ghost values will be
   * interpolated from last.
   *
   * @param patches pairs where the first value is the index in the coarse vector and the second
   * value is a reference to the PatchInfo object
   * @param finer_vector the finer vector
   * @param coarser_vector the coarser vector
   */
  virtual void interpolatePatches(
    const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>>& patches,
    const Vector<D>& coarser_vector,
    Vector<D>& finer_vector) const = 0;

  /**
   * @brief interpolation function
   *
   * @param coarse the input vector that is interpolated from
   * @param fine the output vector that is interpolated to.
   */
  void interpolate(const Vector<D>& coarse, Vector<D>& fine) const
  {
    if constexpr (ENABLE_DEBUG) {
      if (coarse.getNumLocalPatches() != ilc.getCoarserDomain().getNumLocalPatches()) {
        throw RuntimeError("coarse vector is incorrect length. Expected Length of " +
                           std::to_string(ilc.getCoarserDomain().getNumLocalPatches()) +
                           " but vector was length " + std::to_string(coarse.getNumLocalPatches()));
      }
      if (fine.getNumLocalPatches() != ilc.getFinerDomain().getNumLocalPatches()) {
        throw RuntimeError("fine vector is incorrect length. Expected Length of " +
                           std::to_string(ilc.getFinerDomain().getNumLocalPatches()) +
                           " but vector was length " + std::to_string(fine.getNumLocalPatches()));
      }
    }
    Vector<D> coarse_ghost = ilc.getNewGhostVector(coarse.getNumComponents());

    // start scatter for ghost values
    ilc.getGhostPatchesStart(coarse, coarse_ghost);

    // interpolate form local values
    interpolatePatches(ilc.getPatchesWithLocalParent(), coarse, fine);

    // finish scatter for ghost values
    ilc.getGhostPatchesFinish(coarse, coarse_ghost);

    // interpolator from ghost values
    interpolatePatches(ilc.getPatchesWithGhostParent(), coarse_ghost, fine);
  }
};
} // namespace ThunderEgg::GMG
// explicit instantiation
extern template class ThunderEgg::GMG::MPIInterpolator<2>;
extern template class ThunderEgg::GMG::MPIInterpolator<3>;
#endif