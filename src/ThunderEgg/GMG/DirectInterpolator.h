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

#ifndef THUNDEREGG_GMG_DRCTINTP_H
#define THUNDEREGG_GMG_DRCTINTP_H
/**
 * @file
 *
 * @brief DirectInterpolator class
 */

#include <ThunderEgg/GMG/MPIInterpolator.h>

namespace ThunderEgg::GMG {
/**
 * @brief Directly places values from coarse cell into the corresponding fine cells.
 *
 * This is a piecewise constant interpolation scheme.
 */
template<int D>
  requires is_supported_dimension<D>
class DirectInterpolator : public MPIInterpolator<D>
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
   * @brief Create new DirectInterpolator object.
   *
   * @param coarse_domain the coarser Domain
   * @param fine_domain the finer Domain
   * @param num_components the number of components in each cell
   */
  DirectInterpolator(const Domain<D>& coarse_domain, const Domain<D>& fine_domain);

  /**
   * @brief Destroy the DirectInterpolator object
   */
  ~DirectInterpolator();

  /**
   * @brief Copy constructor
   *
   * @param other other DirectInterpolator
   */
  DirectInterpolator(const DirectInterpolator& other);

  /**
   * @brief Move constructor
   *
   * @param other DirectInterpolator
   */
  DirectInterpolator(DirectInterpolator&& other);

  /**
   * @brief copy assignment
   *
   * @param other DirectInterpolator
   * @return DirectInterpolator& this
   */
  DirectInterpolator&
  operator=(const DirectInterpolator& other);

  /**
   * @brief move assignment
   *
   * @param other DirectInterpolator
   * @return DirectInterpolator& this
   */
  DirectInterpolator&
  operator=(DirectInterpolator&& other);

  /**
   * @brief Clone this interpolator
   *
   * @return DirectInterpolator<D>* a newly allocated copy of this interpolator
   */
  DirectInterpolator<D>*
  clone() const override;

  void
  interpolatePatches(
    const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>>& patches,
    const Vector<D>& coarser_vector,
    Vector<D>& finer_vector) const override;
};
extern template class DirectInterpolator<2>;
extern template class DirectInterpolator<3>;
} // namespace ThunderEgg::GMG
#endif