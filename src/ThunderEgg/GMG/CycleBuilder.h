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

#ifndef THUNDEREGG_GMG_CYCLEBUILDER_H
#define THUNDEREGG_GMG_CYCLEBUILDER_H
/**
 * @file
 *
 * @brief CycleBuilder class
 */

#include <ThunderEgg/GMG/FMGCycle.h>
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/GMG/VCycle.h>
#include <ThunderEgg/GMG/WCycle.h>
#include <ThunderEgg/RuntimeError.h>
namespace ThunderEgg::GMG {
/**
 * @brief Builder for GMG cycles.
 *
 * The user provides an Operator, Smoother, Restrictor, and Interpolator object for each level.
 *
 * The user is expected to make these calls in the following order:
 *
 * 		builder.addFinestLevel()
 * 		for each intermediate level (if any)
 * 		{
 * 			builder.addIntermediateLevel()
 * 		}
 * 		builder.addCoarsestLevel()
 * 		M = builder.getCycle();
 *
 * If these are called in the wrong order, an exception will be thrown.
 *
 * @tparam D the number of cartesian dimensions
 */
template<int D>
class CycleBuilder
{
private:
  /**
   * @brief has addFinestLevel been called
   */
  bool has_finest = false;
  /**
   * @brief has addCoarsestLevel been called
   */
  bool has_coarsest = false;
  /**
   * @brief The cycle options
   */
  CycleOpts opts;

  /**
   * @brief the finest level
   */
  std::shared_ptr<Level<D>> finest_level;
  /**
   * @brief the last level that was added
   */
  std::shared_ptr<Level<D>> prev_level;

public:
  /**
   * @brief Construct a new CycleBuilder object
   *
   * @param opts_in the options to use for the GMG cycle
   */
  explicit CycleBuilder(const CycleOpts& opts_in)
    : opts(opts_in)
  {}
  /**
   * @brief Add the finest level to the Cycle
   *
   * @param op the Operator for the level
   * @param smoother the Smoother for the level
   * @param restrictor the Restrictor that restricts from this level to the coarser level
   * @param vg the VectorGenerator for the level
   */
  void addFinestLevel(const Operator<D>& op,
                      const Smoother<D>& smoother,
                      const Restrictor<D>& restrictor)
  {
    if (has_finest) {
      throw RuntimeError("addFinestLevel was already called");
    }
    has_finest = true;

    finest_level = std::make_shared<Level<D>>();
    finest_level->setOperator(op);
    finest_level->setSmoother(smoother);
    finest_level->setRestrictor(restrictor);

    prev_level = finest_level;
  }
  /**
   * @brief Add the next intermediate level to the Cycle
   *
   * @param op the Operator for the level
   * @param smoother the Smoother for the level
   * @param restrictor the Restrictor that restricts from this level to the coarser level
   * @param interpolator the Interpolator that restricts from this level to the finer level
   * @param vg the VectorGenerator for the level
   */
  void addIntermediateLevel(const Operator<D>& op,
                            const Smoother<D>& smoother,
                            const Restrictor<D>& restrictor,
                            const Interpolator<D>& interpolator)
  {
    if (!has_finest) {
      throw RuntimeError("addFinestLevel has not been called yet");
    }
    if (has_coarsest) {
      throw RuntimeError("addCoarsestLevel has been called");
    }

    auto new_level = std::make_shared<Level<D>>();
    new_level->setOperator(op);
    new_level->setSmoother(smoother);
    new_level->setInterpolator(interpolator);
    new_level->setRestrictor(restrictor);

    prev_level->setCoarser(new_level);

    prev_level = new_level;
  }
  /**
   * @brief Add the next intermediate level to the Cycle
   *
   * @param op the Operator for the level
   * @param smoother the Smoother for the level
   * @param interpolator the Interpolator that restricts from this level to the finer level
   */
  void addCoarsestLevel(const Operator<D>& op,
                        Smoother<D>& smoother,
                        const Interpolator<D>& interpolator)
  {
    if (!has_finest) {
      throw RuntimeError("addFinestLevel has not been called yet");
    }
    if (has_coarsest) {
      throw RuntimeError("addCoarsestLevel has already been called");
    }
    has_coarsest = true;

    auto new_level = std::make_shared<Level<D>>();
    new_level->setOperator(op);
    new_level->setSmoother(smoother);
    new_level->setInterpolator(interpolator);

    prev_level->setCoarser(new_level);
  }
  /**
   * @brief Get the completed Cycle object
   *
   * @return std::shared_ptr<Cycle<D>> the completed Cycle, will throw an exception if it is
   * incomplete
   */
  std::shared_ptr<Cycle<D>> getCycle() const
  {
    if (!has_coarsest) {
      throw RuntimeError("addCoarsestLevel has not been called");
    }
    std::shared_ptr<Cycle<D>> cycle;
    if (opts.cycle_type == "V") {
      cycle.reset(new VCycle<D>(*finest_level, opts));
    } else if (opts.cycle_type == "W") {
      cycle.reset(new WCycle<D>(*finest_level, opts));
    } else if (opts.cycle_type == "F") {
      cycle.reset(new FMGCycle<D>(*finest_level, opts));
    } else {
      throw RuntimeError("Unsupported Cycle type: " + opts.cycle_type);
    }
    return cycle;
  }
};
extern template class CycleBuilder<2>;
extern template class CycleBuilder<3>;
} // namespace ThunderEgg::GMG

#endif