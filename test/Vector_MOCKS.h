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

#include <ThunderEgg/Vector.h>
namespace ThunderEgg {
/**
 * @brief Check if a coordinate is a ghost coordinate
 *
 * @tparam D the number of Cartesian Dimensions
 * @param coord the coordinate
 * @param ns the number of cells in each direction of the patch
 * @param num_ghost_cells the number of ghost cells on each side of the patch
 * @return true if it is a ghost coordinate
 * @return false if it is not a ghost coordinate
 */
template<size_t D>
bool
isGhost(const std::array<int, D>& coord, const std::array<int, D>& ns, int num_ghost_cells)
{
  for (size_t i = 0; i < D; i++) {
    if (coord[i] < 0 || coord[i] >= ns[i]) {
      return true;
    }
  }
  return false;
}
} // namespace ThunderEgg