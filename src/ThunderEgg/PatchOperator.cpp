/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2021      Scott Aiton
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

#include <ThunderEgg/PatchOperator.h>

namespace ThunderEgg {
template<int D>
  requires is_supported_dimension<D>
PatchOperator<D>::PatchOperator(const Domain<D>& domain, const GhostFiller<D>& ghost_filler)
  : domain(domain)
  , ghost_filler(ghost_filler.clone())
{
}

template<int D>
  requires is_supported_dimension<D>
PatchOperator<D>::~PatchOperator()
{
}

template<int D>
  requires is_supported_dimension<D>
void
PatchOperator<D>::apply(const Vector<D>& u, Vector<D>& f) const
{
  if constexpr (ENABLE_DEBUG) {
    if (u.getNumLocalPatches() != domain.getNumLocalPatches()) {
      throw RuntimeError("u vector is incorrect length");
    }
    if (f.getNumLocalPatches() != domain.getNumLocalPatches()) {
      throw RuntimeError("f vector is incorrect length");
    }
  }
  f.setWithGhost(0);
  ghost_filler->fillGhost(u);
  for (const PatchInfo<D>& pinfo : domain.getPatchInfoVector()) {
    PatchView<const double, D> u_view = u.getPatchView(pinfo.local_index);
    PatchView<double, D> f_view = f.getPatchView(pinfo.local_index);
    applySinglePatch(pinfo, u_view, f_view);
  }
}

template<int D>
  requires is_supported_dimension<D>
const Domain<D>&
PatchOperator<D>::getDomain() const
{
  return domain;
}

template<int D>
  requires is_supported_dimension<D>
const GhostFiller<D>&
PatchOperator<D>::getGhostFiller() const
{
  return *ghost_filler;
}

// EXPLICIT INSTANTIATIONS

template class PatchOperator<2>;
template class PatchOperator<3>;

}
