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

#include "PatchView.h"
namespace ThunderEgg
{
template class PatchView<double, 1>;
template View<double, 1> PatchView<double, 1>::getSliceOn(Face<1, 0>, const std::array<int, 1> &) const;
template View<double, 1> PatchView<double, 1>::getGhostSliceOn(Face<1, 0>, const std::array<size_t, 1> &) const;

template class PatchView<double, 2>;
template View<double, 1> PatchView<double, 2>::getSliceOn(Face<2, 0>, const std::array<int, 2> &) const;
template View<double, 2> PatchView<double, 2>::getSliceOn(Face<2, 1>, const std::array<int, 1> &) const;
template View<double, 1> PatchView<double, 2>::getGhostSliceOn(Face<2, 0>, const std::array<size_t, 2> &) const;
template View<double, 2> PatchView<double, 2>::getGhostSliceOn(Face<2, 1>, const std::array<size_t, 1> &) const;

template class PatchView<double, 3>;
template View<double, 1> PatchView<double, 3>::getSliceOn(Face<3, 0>, const std::array<int, 3> &) const;
template View<double, 2> PatchView<double, 3>::getSliceOn(Face<3, 1>, const std::array<int, 2> &) const;
template View<double, 3> PatchView<double, 3>::getSliceOn(Face<3, 2>, const std::array<int, 1> &) const;
template View<double, 1> PatchView<double, 3>::getGhostSliceOn(Face<3, 0>, const std::array<size_t, 3> &) const;
template View<double, 2> PatchView<double, 3>::getGhostSliceOn(Face<3, 1>, const std::array<size_t, 2> &) const;
template View<double, 3> PatchView<double, 3>::getGhostSliceOn(Face<3, 2>, const std::array<size_t, 1> &) const;

template class PatchView<const double, 1>;
template View<const double, 1> PatchView<const double, 1>::getSliceOn(Face<1, 0>, const std::array<int, 1> &) const;
template View<double, 1>       PatchView<const double, 1>::getGhostSliceOn(Face<1, 0>, const std::array<size_t, 1> &) const;

template class PatchView<const double, 2>;
template View<const double, 1> PatchView<const double, 2>::getSliceOn(Face<2, 0>, const std::array<int, 2> &) const;
template View<const double, 2> PatchView<const double, 2>::getSliceOn(Face<2, 1>, const std::array<int, 1> &) const;
template View<double, 1>       PatchView<const double, 2>::getGhostSliceOn(Face<2, 0>, const std::array<size_t, 2> &) const;
template View<double, 2>       PatchView<const double, 2>::getGhostSliceOn(Face<2, 1>, const std::array<size_t, 1> &) const;

template class PatchView<const double, 3>;
template View<const double, 1> PatchView<const double, 3>::getSliceOn(Face<3, 0>, const std::array<int, 3> &) const;
template View<const double, 2> PatchView<const double, 3>::getSliceOn(Face<3, 1>, const std::array<int, 2> &) const;
template View<const double, 3> PatchView<const double, 3>::getSliceOn(Face<3, 2>, const std::array<int, 1> &) const;
template View<double, 1>       PatchView<const double, 3>::getGhostSliceOn(Face<3, 0>, const std::array<size_t, 3> &) const;
template View<double, 2>       PatchView<const double, 3>::getGhostSliceOn(Face<3, 1>, const std::array<size_t, 2> &) const;
template View<double, 3>       PatchView<const double, 3>::getGhostSliceOn(Face<3, 2>, const std::array<size_t, 1> &) const;
} // namespace ThunderEgg
