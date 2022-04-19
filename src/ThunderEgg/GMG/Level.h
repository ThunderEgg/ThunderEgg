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

#ifndef THUNDEREGG_GMG_LEVEL_H
#define THUNDEREGG_GMG_LEVEL_H
/**
 * @file
 *
 * @brief Level class
 */
#include <ThunderEgg/GMG/Interpolator.h>
#include <ThunderEgg/GMG/Restrictor.h>
#include <ThunderEgg/GMG/Smoother.h>
#include <ThunderEgg/Operator.h>
#include <memory>

namespace ThunderEgg::GMG {
/**
 * @brief Represents a level in geometric multi-grid.
 */
template<int D>
  requires is_supported_dimension<D>
class Level
{
private:
  /**
   * @brief Implimentation class
   */
  class Implimentation;

  /**
   * @brief pointer to the implimentation
   */
  std::unique_ptr<Implimentation> implimentation;

public:
  /**
   * @brief Create a Level object.
   */
  Level();

  /**
   * @brief Copy constructor
   *
   * @param other level
   */
  Level(const Level<D>& other);

  /**
   * @brief Destroy the Level object
   */
  ~Level();

  /**
   * @brief Move constructor
   *
   * @param other level
   */
  Level(Level<D>&& other);

  /**
   * @brief Copy assignment
   *
   * @param other level
   */
  Level<D>&
  operator=(const Level<D>& other);

  /**
   * @brief Move assignment
   *
   * @param other level
   */
  Level<D>&
  operator=(Level<D>&& other);

  /**
   * @brief Set the restriction operator for restricting from this level to the coarser level.
   *
   * @param restrictor the restriction operator.
   */
  void
  setRestrictor(const Restrictor<D>& restrictor);

  /**
   * @brief Get the restriction operator for this level.
   *
   * @return Reference to the restrictor
   */
  const Restrictor<D>&
  getRestrictor() const;

  /**
   * @brief Set the interpolation operator for interpolating from this level to the finer level.
   *
   * @param interpolator the interpolation operator.
   */
  void
  setInterpolator(const Interpolator<D>& interpolator);

  /**
   * @brief Get the interpolation operator for this level.
   *
   * @return Reference to the interpolator.
   */
  const Interpolator<D>&
  getInterpolator() const;

  /**
   * @brief Set the operator (matrix) for this level.
   *
   * @param op the operator
   */
  void
  setOperator(const Operator<D>& op);

  /**
   * @brief Get the operator for this level.
   *
   * @return Pointer to the operator.
   */
  const Operator<D>&
  getOperator() const;

  /**
   * @brief Set the smoother for this level.
   *
   * @param smoother the smoother
   */
  void
  setSmoother(const Smoother<D>& smoother);

  /**
   * @brief Get smoother operator for this level.
   *
   * @return Pointer to the smoother operator.
   */
  const Smoother<D>&
  getSmoother() const;

  /**
   * @brief Set pointer to the coarser level.
   *
   * @param coarser the pointer to the coarser level.
   */
  void
  setCoarser(std::shared_ptr<const Level> coarser);

  /**
   * @brief get reference to the coarser level.
   *
   * @return reference to the coarser level.
   */
  const Level&
  getCoarser() const;

  /**
   * @brief Check if this level is the finest level.
   *
   * @return whether or not this level is the finest level.
   */
  bool
  finest() const;

  /**
   * @brief Check if this level is the coarsest level.
   *
   * @return whether or not this level is the coarsest level.
   */
  bool
  coarsest() const;
};
extern template class Level<2>;
extern template class Level<3>;
} // namespace ThunderEgg::GMG
#endif