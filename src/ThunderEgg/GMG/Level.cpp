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

#include <ThunderEgg/GMG/Level.h>

namespace ThunderEgg::GMG {

/**
 * @brief Represents a level in geometric multi-grid.
 */
template<int D>
  requires is_supported_dimension<D>
class Level<D>::Implimentation
{
private:
  /**
   * @brief The operator (matrix) for this level.
   */
  std::shared_ptr<const Operator<D>> op;
  /**
   * @brief The restrictor from this level to the coarser level.
   */
  std::shared_ptr<const Restrictor<D>> restrictor;
  /**
   * @brief The interpolator from this level to the finer level.
   */
  std::shared_ptr<const Interpolator<D>> interpolator;
  /**
   * @brief The smoother for this level.
   */
  std::shared_ptr<const Smoother<D>> smoother;
  /**
   * @brief Pointer to coarser level
   */
  std::shared_ptr<const Level> coarser;

public:
  /**
   * @brief Set the restriction operator for restricting from this level to the coarser level.
   *
   * @param restrictor the restriction operator.
   */
  void
  setRestrictor(const Restrictor<D>& restrictor)
  {
    this->restrictor.reset(restrictor.clone());
  }

  /**
   * @brief Get the restriction operator for this level.
   *
   * @return Reference to the restrictor
   */
  const Restrictor<D>&
  getRestrictor() const
  {
    if (restrictor == nullptr) {
      throw RuntimeError("This level does not have a restrictor");
    }
    return *restrictor;
  }

  /**
   * @brief Set the interpolation operator for interpolating from this level to the finer level.
   *
   * @param interpolator the interpolation operator.
   */
  void
  setInterpolator(const Interpolator<D>& interpolator)
  {
    this->interpolator.reset(interpolator.clone());
  }

  /**
   * @brief Get the interpolation operator for this level.
   *
   * @return Reference to the interpolator.
   */
  const Interpolator<D>&
  getInterpolator() const
  {
    if (interpolator == nullptr) {
      throw RuntimeError("This level does not have an interpolator");
    }
    return *interpolator;
  }

  /**
   * @brief Set the operator (matrix) for this level.
   *
   * @param op the operator
   */
  void
  setOperator(const Operator<D>& op)
  {
    this->op.reset(op.clone());
  }

  /**
   * @brief Get the operator for this level.
   *
   * @return Pointer to the operator.
   */
  const Operator<D>&
  getOperator() const
  {
    if (op == nullptr) {
      throw RuntimeError("This level does not have an Operator");
    }
    return *op;
  }

  /**
   * @brief Set the smoother for this level.
   *
   * @param smoother the smoother
   */
  void
  setSmoother(const Smoother<D>& smoother)
  {
    this->smoother.reset(smoother.clone());
  }

  /**
   * @brief Get smoother operator for this level.
   *
   * @return Pointer to the smoother operator.
   */
  const Smoother<D>&
  getSmoother() const
  {
    if (smoother == nullptr) {
      throw RuntimeError("This level does not have a smoother");
    }
    return *smoother;
  }

  /**
   * @brief Set pointer to the coarser level.
   *
   * @param coarser the pointer to the coarser level.
   */
  void
  setCoarser(std::shared_ptr<const Level> coarser)
  {
    this->coarser = coarser;
  }

  /**
   * @brief get reference to the coarser level.
   *
   * @return reference to the coarser level.
   */
  const Level&
  getCoarser() const
  {
    if (coarser == nullptr) {
      throw RuntimeError("This level does not have a coarser level.");
    }
    return *coarser;
  }

  /**
   * @brief Check if this level is the finest level.
   *
   * @return whether or not this level is the finest level.
   */
  bool
  finest() const
  {
    return interpolator == nullptr;
  }

  /**
   * @brief Check if this level is the coarsest level.
   *
   * @return whether or not this level is the coarsest level.
   */
  bool
  coarsest() const
  {
    return coarser == nullptr;
  }
};

template<int D>
  requires is_supported_dimension<D>
Level<D>::Level()
  : implimentation(new Implimentation())
{
}

template<int D>
  requires is_supported_dimension<D>
Level<D>::~Level() = default;

template<int D>
  requires is_supported_dimension<D>
Level<D>::Level(const Level<D>& other)
  : implimentation(new Implimentation(*other.implimentation))
{
}

template<int D>
  requires is_supported_dimension<D>
Level<D>::Level(Level<D>&& other) = default;

template<int D>
  requires is_supported_dimension<D>
Level<D>&
Level<D>::operator=(const Level<D>& other)
{
  implimentation.reset(new Implimentation(*other.implimentation));
  return *this;
}

template<int D>
  requires is_supported_dimension<D>
Level<D>&
Level<D>::operator=(Level<D>&& other) = default;

template<int D>
  requires is_supported_dimension<D>
void
Level<D>::setRestrictor(const Restrictor<D>& restrictor)
{
  implimentation->setRestrictor(restrictor);
}

template<int D>
  requires is_supported_dimension<D>
const Restrictor<D>&
Level<D>::getRestrictor() const
{
  return implimentation->getRestrictor();
}

template<int D>
  requires is_supported_dimension<D>
void
Level<D>::setInterpolator(const Interpolator<D>& interpolator)
{
  implimentation->setInterpolator(interpolator);
}

template<int D>
  requires is_supported_dimension<D>
const Interpolator<D>&
Level<D>::getInterpolator() const
{
  return implimentation->getInterpolator();
}

template<int D>
  requires is_supported_dimension<D>
void
Level<D>::setOperator(const Operator<D>& op)
{
  implimentation->setOperator(op);
}

template<int D>
  requires is_supported_dimension<D>
const Operator<D>&
Level<D>::getOperator() const
{
  return implimentation->getOperator();
}

template<int D>
  requires is_supported_dimension<D>
void
Level<D>::setSmoother(const Smoother<D>& smoother)
{
  implimentation->setSmoother(smoother);
}

template<int D>
  requires is_supported_dimension<D>
const Smoother<D>&
Level<D>::getSmoother() const
{
  return implimentation->getSmoother();
}

template<int D>
  requires is_supported_dimension<D>
void
Level<D>::setCoarser(std::shared_ptr<const Level> coarser)
{
  implimentation->setCoarser(coarser);
}

template<int D>
  requires is_supported_dimension<D>
const Level<D>&
Level<D>::getCoarser() const
{
  return implimentation->getCoarser();
}

template<int D>
  requires is_supported_dimension<D> bool
Level<D>::finest() const
{
  return implimentation->finest();
}

template<int D>
  requires is_supported_dimension<D> bool
Level<D>::coarsest() const
{
  return implimentation->coarsest();
}

template class ThunderEgg::GMG::Level<2>;
template class ThunderEgg::GMG::Level<3>;
} // namespace ThunderEgg::GMG