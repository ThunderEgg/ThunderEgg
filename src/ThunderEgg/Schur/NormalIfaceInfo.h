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

#ifndef THUNDEREGG_SCHUR_NORMALIFACEINFO_H
#define THUNDEREGG_SCHUR_NORMALIFACEINFO_H
/**
 * @file
 *
 * @brief NormalIfaceInfo class
 */
#include <ThunderEgg/PatchInfo.h>
#include <ThunderEgg/Schur/IfaceInfo.h>
namespace ThunderEgg {
namespace Schur {
/**
 * @brief This represents an interface where the neighbor is at the same refinement level
 *
 * @tparam D the number of Cartesian dimensions in a patch
 */
template<int D>
class NormalIfaceInfo : public IfaceInfo<D>
{
private:
  /**
   * @brief Get the rank for the interface on a given side of the patch
   *
   * @param pinfo the patch
   * @param s the side
   * @return int the rank
   */
  static int GetRank(const PatchInfo<D>& pinfo, Side<D> s)
  {
    if (s.isLowerOnAxis()) {
      // lower axis interface belongs to neighboring rank
      auto nbr_info = pinfo.getNormalNbrInfo(s);
      return nbr_info.rank;
    } else {
      // higher axis interafce belongs to this patch's rank
      return pinfo.rank;
    }
  }

  /**
   * @brief Get the id for the interface on a given side of the patch
   *
   * @param pinfo the patch
   * @param s the side
   * @return int the id
   */
  static int GetId(const PatchInfo<D>& pinfo, Side<D> s)
  {
    if (s.isLowerOnAxis()) {
      // lower axis interface belongs to neighboring rank
      auto nbr_info = pinfo.getNormalNbrInfo(s);
      return (int)(nbr_info.id * Side<D>::number_of + s.opposite().getIndex());
    } else {
      // higher axis interafce belongs to this patch's rank
      return (int)(pinfo.id * Side<D>::number_of + s.getIndex());
    }
  }

public:
  /**
   * @brief Construct a new NormalIfaceInfo object
   *
   * @param pinfo the associated PatchInfo object
   * @param s the side of the patch that the interface is on
   */
  NormalIfaceInfo(const PatchInfo<D>& pinfo, Side<D> s)
    : IfaceInfo<D>(GetRank(pinfo, s), GetId(pinfo, s))
  {}
};
} // namespace Schur
} // namespace ThunderEgg
#endif
