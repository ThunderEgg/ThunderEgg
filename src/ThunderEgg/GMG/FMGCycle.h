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

#ifndef THUNDEREGG_GMG_FMGCYCLE_H
#define THUNDEREGG_GMG_FMGCYCLE_H
/**
 * @file
 *
 * @brief FMGCycle class
 */
#include <ThunderEgg/GMG/Cycle.h>
#include <ThunderEgg/GMG/CycleOpts.h>

namespace ThunderEgg::GMG {
/**
 * @brief Implementation of a full multigrid cycle
 */
template<int D>
  requires is_supported_dimension<D>
class FMGCycle : public Cycle<D>
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

  void
  visit(const Level<D>& level, const Vector<D>& f, Vector<D>& u) const override;

public:
  /**
   * @brief Create new FMGCycle
   *
   * @param finest_level the finest level
   * @param opts the options for the cycle pre, post, coarse, and mid sweeps are used
   */
  FMGCycle(const Level<D>& finest_level, const CycleOpts& opts);

  /**
   * @brief Get a clone of this FMGCycle
   *
   * @return FMGCycle<D>* a newly allocated copy
   */
  FMGCycle<D>*
  clone() const override;
};
extern template class FMGCycle<2>;
extern template class FMGCycle<3>;
} // namespace ThunderEgg::GMG
#endif