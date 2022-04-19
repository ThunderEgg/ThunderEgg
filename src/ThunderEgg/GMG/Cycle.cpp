/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2019-2021 Scott Aiton
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

#include <ThunderEgg/GMG/Cycle.h>

namespace ThunderEgg::GMG {

/**
 * @brief Cycle Implimentation
 */
template<int D>
  requires is_supported_dimension<D>
class Cycle<D>::Implimentation
{
private:
  /**
   * @brief pointer to the finest level
   */
  std::shared_ptr<const Level<D>> finest_level;

public:
  /**
   * @brief Create new cycle object.
   *
   * @param finest_level the finest level object.
   */
  Implimentation(const Level<D>& finest_level)
    : finest_level(new Level<D>(finest_level))
  {
  }

  /**
   * @brief Run one iteration of the cycle.
   *
   * Performs one cycle on on the system `Au=f` where `A` is the operator for the
   * finest level.
   *
   * @param f the RHS vector.
   * @param u the solution vector.
   */
  void
  apply(const Cycle<D>& cycle, const Vector<D>& f, Vector<D>& u) const
  {
    u.setWithGhost(0);
    cycle.visit(*finest_level, f, u);
  }

  /**
   * @brief Get the finest Level
   *
   * @return const Level<D>& the Level
   */
  const Level<D>&
  getFinestLevel() const
  {
    return *finest_level;
  }

  /**
   * @brief Prepare vectors for coarser level.
   *
   * @param level the current level
   * @param f the rhs vector cooresponding to the level
   * @param u the solution vector cooresponding to the level
   * @return Vector<D> the restricted residual vector
   */
  Vector<D> restrict(const Level<D>& level, const Vector<D>& f, const Vector<D>& u) const
  {
    // calculate residual
    Vector<D> r = u.getZeroClone();
    level.getOperator().apply(u, r);
    r.scaleThenAdd(-1, f);
    // create vectors for coarser levels
    return level.getRestrictor().restrict(r);
  }
};

template<int D>
  requires is_supported_dimension<D>
Vector<D> Cycle<D>::restrict(const Level<D>& level, const Vector<D>& f, const Vector<D>& u) const
{
  return implimentation->restrict(level, f, u);
}

template<int D>
  requires is_supported_dimension<D>
Cycle<D>::Cycle(const Level<D>& finest_level)
  : implimentation(new Implimentation(finest_level))
{
}

template<int D>
  requires is_supported_dimension<D>
void
Cycle<D>::apply(const Vector<D>& f, Vector<D>& u) const
{
  implimentation->apply(*this, f, u);
}

template<int D>
  requires is_supported_dimension<D>
const Level<D>&
Cycle<D>::getFinestLevel() const
{
  return implimentation->getFinestLevel();
}

template class Cycle<2>;
template class Cycle<3>;
} // namespace ThunderEgg::GMG