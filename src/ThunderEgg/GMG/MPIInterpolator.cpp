/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020      Scott Aiton
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

#include <ThunderEgg/GMG/MPIInterpolator.h>

namespace ThunderEgg::GMG {

template<int D>
  requires is_supported_dimension<D>
class MPIInterpolator<D>::Implimentation
{
private:
  /**
   * @brief The communication package for restricting between levels.
   */
  InterLevelComm<D> ilc;

public:
  /**
   * @brief Create new MPIInterpolator object.
   *
   * @param ilc the communcation package for the two levels.
   */
  Implimentation(const Domain<D>& coarser_domain, const Domain<D>& finer_domain)
    : ilc(coarser_domain, finer_domain)
  {
  }

  /**
   * @brief interpolation function
   *
   * @param coarse the input vector that is interpolated from
   * @param fine the output vector that is interpolated to.
   */
  void
  interpolate(const MPIInterpolator<D>& interpolator,
              const Vector<D>& coarse,
              Vector<D>& fine) const
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

    InterLevelComm my_ilc = ilc;
    // start scatter for ghost values
    my_ilc.getGhostPatchesStart(coarse, coarse_ghost);

    // interpolate form local values
    interpolator.interpolatePatches(ilc.getPatchesWithLocalParent(), coarse, fine);

    // finish scatter for ghost values
    my_ilc.getGhostPatchesFinish(coarse, coarse_ghost);

    // interpolator from ghost values
    interpolator.interpolatePatches(ilc.getPatchesWithGhostParent(), coarse_ghost, fine);
  }
};

template<int D>
  requires is_supported_dimension<D>
MPIInterpolator<D>::MPIInterpolator(const Domain<D>& coarser_domain, const Domain<D>& finer_domain)
  : implimentation(new Implimentation(coarser_domain, finer_domain))
{
}

template<int D>
  requires is_supported_dimension<D>
MPIInterpolator<D>::~MPIInterpolator() = default;

template<int D>
  requires is_supported_dimension<D>
MPIInterpolator<D>::MPIInterpolator(const MPIInterpolator& other)
  : implimentation(new Implimentation(*other.implimentation))
{
}

template<int D>
  requires is_supported_dimension<D>
MPIInterpolator<D>::MPIInterpolator(MPIInterpolator&& other) = default;

template<int D>
  requires is_supported_dimension<D>
MPIInterpolator<D>&
MPIInterpolator<D>::operator=(const MPIInterpolator& other)
{
  implimentation.reset(new Implimentation(*other.implimentation));
  return *this;
}

template<int D>
  requires is_supported_dimension<D>
MPIInterpolator<D>&
MPIInterpolator<D>::operator=(MPIInterpolator&& other) = default;

template<int D>
  requires is_supported_dimension<D>
void
MPIInterpolator<D>::interpolate(const Vector<D>& coarse, Vector<D>& fine) const
{
  implimentation->interpolate(*this, coarse, fine);
}

template class MPIInterpolator<2>;
template class MPIInterpolator<3>;
} // namespace ThunderEgg::GMG