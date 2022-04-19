/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2018-2021 Scott Aiton
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

#ifndef THUNDEREGG_GMG_WCYCLE_H
#define THUNDEREGG_GMG_WCYCLE_H
/**
 * @file
 *
 * @brief WCycle class
 */
#include <ThunderEgg/GMG/Cycle.h>
#include <ThunderEgg/GMG/CycleOpts.h>

namespace ThunderEgg::GMG {
/**
 * @brief Implementation of a W-cycle
 */
template<int D>
  requires is_supported_dimension<D>
class WCycle : public Cycle<D>
{
private:
  /**
   * @brief Implimentation class
   */
  class Implimentation;

  /**
   * @brief pointer to the implimentation
   */
  std::shared_ptr<const Implimentation> implimentation;

protected:
  /**
   * @brief Implements W-cycle. Pre-smooth, visit coarser level, smooth, visit coarse level, and
   * then post-smooth.
   *
   * @param level the current level that is being visited.
   */
  void
  visit(const Level<D>& level, const Vector<D>& f, Vector<D>& u) const override;

public:
  /**
   * @brief Create new W-cycle
   *
   * @param finest_level a pointer to the finest level
   */
  WCycle(const Level<D>& finest_level, const CycleOpts& opts);

  /**
   * @brief Get a clone of this WCycle
   *
   * @return WCycle<D>* a newly allocated copy
   */
  WCycle<D>*
  clone() const override;
};
extern template class WCycle<2>;
extern template class WCycle<3>;
} // namespace ThunderEgg::GMG
#endif