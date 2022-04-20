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

#ifndef THUNDEREGG_GMG_LINEARRESTRICTOR_H
#define THUNDEREGG_GMG_LINEARRESTRICTOR_H
/**
 * @file
 *
 * @brief LinearRestrictor class
 */
#include <ThunderEgg/GMG/MPIRestrictor.h>

namespace ThunderEgg::GMG {
/**
 * @brief Restrictor that averages the corresponding fine cells into each coarse cell.
 */
template<int D>
  requires is_supported_dimension<D>
class LinearRestrictor : public MPIRestrictor<D>
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
   * @brief Create new LinearRestrictor object.
   *
   * @param fine_domain the finer Domain
   * @param coarse_domain the coarser Domain
   * @param extrapolate_boundary_ghosts set to true if ghost values at the boundaries should be
   * extrapolated
   */
  LinearRestrictor(const Domain<D>& fine_domain,
                   const Domain<D>& coarse_domain,
                   bool extrapolate_boundary_ghosts = false);

  /**
   * @brief Destroy the LinearRestrictor object
   */
  ~LinearRestrictor();

  /**
   * @brief Copy constructor
   *
   * @param other other LinearRestrictor
   */
  LinearRestrictor(const LinearRestrictor& other);

  /**
   * @brief Move constructor
   *
   * @param other LinearRestrictor
   */
  LinearRestrictor(LinearRestrictor&& other);

  /**
   * @brief copy assignment
   *
   * @param other LinearRestrictor
   * @return LinearRestrictor& this
   */
  LinearRestrictor&
  operator=(const LinearRestrictor& other);

  /**
   * @brief move assignment
   *
   * @param other LinearRestrictor
   * @return LinearRestrictor& this
   */
  LinearRestrictor&
  operator=(LinearRestrictor&& other);

  /**
   * @brief Clone this restrictor
   *
   * @return LinearRestrictor<D>* a newly allocated copy of this restrictor
   */
  LinearRestrictor<D>*
  clone() const override;

  void
  restrictPatches(
    const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>>& patches,
    const Vector<D>& finer_vector,
    Vector<D>& coarser_vector) const override;
};
} // namespace ThunderEgg::GMG
// explicit instantiation
extern template class ThunderEgg::GMG::LinearRestrictor<2>;
extern template class ThunderEgg::GMG::LinearRestrictor<3>;
#endif