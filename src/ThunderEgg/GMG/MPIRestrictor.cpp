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

#include <ThunderEgg/GMG/MPIRestrictor.h>

namespace ThunderEgg::GMG {
template<int D>
  requires is_supported_dimension<D>
class MPIRestrictor<D>::Implimentation
{
private:
  /**
   * @brief The communication package for restricting between levels.
   */
  InterLevelComm<D> ilc;

public:
  /**
   * @brief Create new LinearRestrictor object.
   *
   * @param ilc the communcation package for the two levels.
   */
  Implimentation(const Domain<D>& coarser_domain, const Domain<D>& finer_domain)
    : ilc(coarser_domain, finer_domain)
  {
  }

  Vector<D> restrict(const MPIRestrictor<D>& restrictor, const Vector<D>& fine) const
  {
    if constexpr (ENABLE_DEBUG) {
      if (fine.getNumLocalPatches() != ilc.getFinerDomain().getNumLocalPatches()) {
        throw RuntimeError("fine vector is incorrect length. Expected Length of " +
                           std::to_string(ilc.getFinerDomain().getNumLocalPatches()) +
                           " but vector was length " + std::to_string(fine.getNumLocalPatches()));
      }
    }
    Vector<D> coarse(ilc.getCoarserDomain(), fine.getNumComponents());
    Vector<D> coarse_ghost = ilc.getNewGhostVector(fine.getNumComponents());

    // fill in ghost values
    restrictor.restrictPatches(ilc.getPatchesWithGhostParent(), fine, coarse_ghost);

    // clear values in coarse vector
    coarse.setWithGhost(0);

    InterLevelComm my_ilc = ilc;
    // start scatter for ghost values
    my_ilc.sendGhostPatchesStart(coarse, coarse_ghost);

    // fill in local values
    restrictor.restrictPatches(ilc.getPatchesWithLocalParent(), fine, coarse);

    // finish scatter for ghost values
    my_ilc.sendGhostPatchesFinish(coarse, coarse_ghost);

    return coarse;
  }
};

template<int D>
  requires is_supported_dimension<D>
MPIRestrictor<D>::MPIRestrictor(const Domain<D>& coarser_domain, const Domain<D>& finer_domain)
  : implimentation(new Implimentation(coarser_domain, finer_domain))
{
}

template<int D>
  requires is_supported_dimension<D>
MPIRestrictor<D>::~MPIRestrictor() = default;

template<int D>
  requires is_supported_dimension<D>
MPIRestrictor<D>::MPIRestrictor(const MPIRestrictor<D>& other)
  : implimentation(new Implimentation(*other.implimentation))
{
}

template<int D>
  requires is_supported_dimension<D>
MPIRestrictor<D>::MPIRestrictor(MPIRestrictor<D>&& other) = default;

template<int D>
  requires is_supported_dimension<D>
MPIRestrictor<D>&
MPIRestrictor<D>::operator=(const MPIRestrictor<D>& other)
{
  implimentation.reset(new Implimentation(*other.implimentation));
  return *this;
}

template<int D>
  requires is_supported_dimension<D>
MPIRestrictor<D>&
MPIRestrictor<D>::operator=(MPIRestrictor<D>&& other) = default;

template<int D>
  requires is_supported_dimension<D>
Vector<D> MPIRestrictor<D>::restrict(const Vector<D>& fine) const
{
  return implimentation->restrict(*this, fine);
}

template class MPIRestrictor<2>;
template class MPIRestrictor<3>;
} // namespace ThunderEgg::GMG