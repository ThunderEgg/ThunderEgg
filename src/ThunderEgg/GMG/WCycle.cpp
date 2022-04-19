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

#include <ThunderEgg/GMG/WCycle.h>

namespace ThunderEgg::GMG {
/**
 * @brief Implementation of a W-cycle
 */
template<int D>
  requires is_supported_dimension<D>
class WCycle<D>::Implimentation
{
private:
  int num_pre_sweeps = 1;
  int num_post_sweeps = 1;
  int num_coarse_sweeps = 1;
  int num_mid_sweeps = 1;

public:
  /**
   * @brief Create new W-cycle
   *
   * @param finest_level a pointer to the finest level
   */
  Implimentation(const Level<D>& finest_level, const CycleOpts& opts)
    : num_pre_sweeps(opts.pre_sweeps)
    , num_post_sweeps(opts.post_sweeps)
    , num_coarse_sweeps(opts.coarse_sweeps)
    , num_mid_sweeps(opts.mid_sweeps)
  {
  }

  /**
   * @brief Implements W-cycle. Pre-smooth, visit coarser level, smooth, visit coarse level, and
   * then post-smooth.
   */
  void
  visit(const WCycle<D>& cycle, const Level<D>& level, const Vector<D>& f, Vector<D>& u) const
  {
    if (level.coarsest()) {
      for (int i = 0; i < num_coarse_sweeps; i++) {
        level.getSmoother().smooth(f, u);
      }
    } else {
      for (int i = 0; i < num_pre_sweeps; i++) {
        level.getSmoother().smooth(f, u);
      }

      Vector<D> coarser_f = cycle.restrict(level, f, u);

      const Level<D>& coarser_level = level.getCoarser();
      Vector<D> coarser_u = coarser_f.getZeroClone();

      this->visit(cycle, coarser_level, coarser_f, coarser_u);

      coarser_level.getInterpolator().interpolate(coarser_u, u);

      for (int i = 0; i < num_mid_sweeps; i++) {
        level.getSmoother().smooth(f, u);
      }

      coarser_f = cycle.restrict(level, f, u);

      this->visit(cycle, coarser_level, coarser_f, coarser_u);

      coarser_level.getInterpolator().interpolate(coarser_u, u);

      for (int i = 0; i < num_post_sweeps; i++) {
        level.getSmoother().smooth(f, u);
      }
    }
  }
};

template<int D>
  requires is_supported_dimension<D>
void
WCycle<D>::visit(const Level<D>& level, const Vector<D>& f, Vector<D>& u) const
{
  implimentation->visit(*this, level, f, u);
}

template<int D>
  requires is_supported_dimension<D>
WCycle<D>::WCycle(const Level<D>& finest_level, const CycleOpts& opts)
  : Cycle<D>(finest_level)
  , implimentation(new Implimentation(finest_level, opts))
{
}

template<int D>
  requires is_supported_dimension<D>
WCycle<D>*
WCycle<D>::clone() const
{
  return new WCycle<D>(*this);
}

template class WCycle<2>;
template class WCycle<3>;
} // namespace ThunderEgg::GMG