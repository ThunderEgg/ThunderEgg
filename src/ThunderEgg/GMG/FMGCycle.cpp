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

#include <ThunderEgg/GMG/FMGCycle.h>

namespace ThunderEgg::GMG {
/**
 * @brief Implementation of a full multigrid cycle
 */
template<int D>
  requires is_supported_dimension<D>
class FMGCycle<D>::Implimentation
{
private:
  /**
   * @brief Number of pre sweeps
   */
  int num_pre_sweeps = 1;

  /**
   * @brief Number of post sweeps
   */
  int num_post_sweeps = 1;

  /**
   * @brief Number of coarse sweeps
   */
  int num_coarse_sweeps = 1;

  /**
   * @brief Number of middle sweeps
   */
  int num_mid_sweeps = 1;

  /**
   * @brief vcycle visit
   */
  void
  v_visit(const FMGCycle<D>& cycle, const Level<D>& level, const Vector<D>& f, Vector<D>& u) const
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

      for (int i = 0; i < num_post_sweeps; i++) {
        level.getSmoother().smooth(f, u);
      }
    }
  }

public:
  /**
   * @brief Create new FMGCycle
   *
   * @param finest_level the finest level
   * @param opts the options for the cycle pre, post, coarse, and mid sweeps are used
   */
  Implimentation(const Level<D>& finest_level, const CycleOpts& opts)
    : num_pre_sweeps(opts.pre_sweeps)
    , num_post_sweeps(opts.post_sweeps)
    , num_coarse_sweeps(opts.coarse_sweeps)
    , num_mid_sweeps(opts.mid_sweeps)
  {
  }

  void
  visit(const FMGCycle<D>& cycle, const Level<D>& level, const Vector<D>& f, Vector<D>& u) const
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

      v_visit(cycle, coarser_level, coarser_f, coarser_u);

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
FMGCycle<D>::visit(const Level<D>& level, const Vector<D>& f, Vector<D>& u) const
{
  implimentation->visit(*this, level, f, u);
}

template<int D>
  requires is_supported_dimension<D>
FMGCycle<D>::FMGCycle(const Level<D>& finest_level, const CycleOpts& opts)
  : Cycle<D>(finest_level)
  , implimentation(new Implimentation(finest_level, opts))
{
}

template<int D>
  requires is_supported_dimension<D>
FMGCycle<D>*
FMGCycle<D>::clone() const
{
  return new FMGCycle<D>(*this);
}

template class FMGCycle<2>;
template class FMGCycle<3>;
} // namespace ThunderEgg::GMG