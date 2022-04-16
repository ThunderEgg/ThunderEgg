/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2019-2021 Scott Aiton
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

#include <ThunderEgg/PatchSolver.h>

namespace ThunderEgg {

template<int D>
  requires is_supported_dimension<D>
PatchSolver<D>::PatchSolver(const Domain<D>& domain, const GhostFiller<D>& ghost_filler)
  : domain(domain)
  , ghost_filler(ghost_filler.clone())
{
}

template<int D>
  requires is_supported_dimension<D>
PatchSolver<D>::~PatchSolver()
{
}

template<int D>
  requires is_supported_dimension<D>
const Domain<D>&
PatchSolver<D>::getDomain() const
{
  return domain;
}

template<int D>
  requires is_supported_dimension<D>
const GhostFiller<D>&
PatchSolver<D>::getGhostFiller() const
{
  return *ghost_filler;
}

template<int D>
  requires is_supported_dimension<D>
void
PatchSolver<D>::apply(const Vector<D>& f, Vector<D>& u) const
{
  if constexpr (ENABLE_DEBUG) {
    if (u.getNumLocalPatches() != this->domain.getNumLocalPatches()) {
      throw RuntimeError("u vector is incorrect length");
    }
    if (f.getNumLocalPatches() != this->domain.getNumLocalPatches()) {
      throw RuntimeError("f vector is incorrect length");
    }
  }
  u.setWithGhost(0);
  if (domain.hasTimer()) {
    domain.getTimer()->startDomainTiming(domain.getId(), "Total Patch Solve");
  }
  for (const PatchInfo<D>& pinfo : domain.getPatchInfoVector()) {
    if (domain.hasTimer()) {
      domain.getTimer()->start("Single Patch Solve");
    }
    PatchView<const double, D> f_view = f.getPatchView(pinfo.local_index);
    PatchView<double, D> u_view = u.getPatchView(pinfo.local_index);
    solveSinglePatch(pinfo, f_view, u_view);
    if (domain.hasTimer()) {
      domain.getTimer()->stop("Single Patch Solve");
    }
  }
  if (domain.hasTimer()) {
    domain.getTimer()->stopDomainTiming(domain.getId(), "Total Patch Solve");
  }
}

template<int D>
  requires is_supported_dimension<D>
void
PatchSolver<D>::smooth(const Vector<D>& f, Vector<D>& u) const
{
  if constexpr (ENABLE_DEBUG) {
    if (u.getNumLocalPatches() != this->domain.getNumLocalPatches()) {
      throw RuntimeError("u vector is incorrect length");
    }
    if (f.getNumLocalPatches() != this->domain.getNumLocalPatches()) {
      throw RuntimeError("f vector is incorrect length");
    }
  }
  if (domain.hasTimer()) {
    domain.getTimer()->startDomainTiming(domain.getId(), "Total Patch Smooth");
  }
  ghost_filler->fillGhost(u);
  for (const PatchInfo<D>& pinfo : domain.getPatchInfoVector()) {
    if (domain.hasTimer()) {
      domain.getTimer()->startPatchTiming(pinfo.id, domain.getId(), "Single Patch Solve");
    }
    PatchView<const double, D> f_view = f.getPatchView(pinfo.local_index);
    PatchView<double, D> u_view = u.getPatchView(pinfo.local_index);
    solveSinglePatch(pinfo, f_view, u_view);
    if (domain.hasTimer()) {
      domain.getTimer()->stopPatchTiming(pinfo.id, domain.getId(), "Single Patch Solve");
    }
  }
  if (domain.hasTimer()) {
    domain.getTimer()->stopDomainTiming(domain.getId(), "Total Patch Smooth");
  }
}

// EXPLICIT INSTANTIATIONS

template class PatchSolver<2>;
template class PatchSolver<3>;

} // namespace ThunderEgg