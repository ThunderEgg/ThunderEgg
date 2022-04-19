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

#ifndef THUNDEREGG_GMG_VCYCLE_H
#define THUNDEREGG_GMG_VCYCLE_H
/**
 * @file
 *
 * @brief VCycle class
 */
#include <ThunderEgg/GMG/Cycle.h>
#include <ThunderEgg/GMG/CycleOpts.h>

namespace ThunderEgg::GMG {
/**
 * @brief Implementation of a V-cycle
 */
template<int D>
  requires is_supported_dimension<D>
class VCycle : public Cycle<D>
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
  /*
   * @brief Implements V-cycle. Pre-smooth, visit coarser level and then post-smooth.
   */
  void
  visit(const Level<D>& level, const Vector<D>& f, Vector<D>& u) const override;

public:
  /**
   * @brief Create new V-cycle
   *
   * @param finest_level a pointer to the finest level
   */
  VCycle(const Level<D>& finest_level, const CycleOpts& opts);

  /**
   * @brief Get a clone of this VCycle
   *
   * @return VCycle<D>* a newly allocated copy
   */
  VCycle<D>*
  clone() const override;
};
extern template class VCycle<2>;
extern template class VCycle<3>;
} // namespace ThunderEgg::GMG
#endif