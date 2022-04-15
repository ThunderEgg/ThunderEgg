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

#ifndef THUNDEREGG_PATCHOPERATOR_H
#define THUNDEREGG_PATCHOPERATOR_H
/**
 * @file
 *
 * @brief PatchOperator class
 */
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/Operator.h>
#include <ThunderEgg/Vector.h>

namespace ThunderEgg {
/**
 * @brief This is an Operator where derived classes only have to implement the two virtual functions
 * that operate on single patch.
 *
 * @tparam D the number of Cartesian dimensions.
 */
template<int D>
  requires is_supported_dimension<D>
class PatchOperator : public Operator<D>
{
private:
  /**
   * @brief the domain that is being solved over
   */
  Domain<D> domain;
  /**
   * @brief The ghost filler, needed for smoothing
   */
  std::shared_ptr<const GhostFiller<D>> ghost_filler;

public:
  /**
   * @brief Construct a new Patch Operator object
   *
   *  This sets the Domain and GhostFiller
   *
   * @param domain  the Domain
   * @param ghost_filler the GhostFiller
   */
  PatchOperator(const Domain<D>& domain, const GhostFiller<D>& ghost_filler);

  /**
   * @brief Clone this patch operator
   *
   * @return PatchOperator<D>* a newly allocated copy of this patch operator
   */
  virtual PatchOperator<D>*
  clone() const override = 0;

  /**
   * @brief Destroy the PatchOperator object
   */
  virtual ~PatchOperator();

  /**
   * @brief Apply the operator to a single patch
   *
   * The ghost values in u will be updated to the latest values, and should not need to be modified
   *
   * @param pinfo  the patch
   * @param u_view the solution
   * @param f_view the left hand side
   */
  virtual void
  applySinglePatch(const PatchInfo<D>& pinfo,
                   const PatchView<const double, D>& u_view,
                   const PatchView<double, D>& f_view) const = 0;

  /**
   * @brief Treat the internal patch boundaries as domain boundaires and modify
   * RHS accordingly.
   *
   * This will be u_viewed in patch solvers to formulate a RHS for the individual patch to solve
   * for.
   *
   * @param pinfo the patch
   * @param u_view the left hand side
   * @param f_view the right hand side
   */
  virtual void
  modifyRHSForInternalBoundaryConditions(const PatchInfo<D>& pinfo,
                                         const PatchView<const double, D>& u_view,
                                         const PatchView<double, D>& f_view) const = 0;

  /**
   * @brief Apply the operator to a single patch
   *
   * The ghost values in u will be updated to the latest values, and should not need to be modified
   *
   * @param pinfo  the patch
   * @param u_view the solution
   * @param f_view the left hand side
   */
  virtual void
  applySinglePatchWithInternalBoundaryConditions(const PatchInfo<D>& pinfo,
                                                 const PatchView<const double, D>& u_view,
                                                 const PatchView<double, D>& f_view) const = 0;

  /**
   * @brief Apply the operator
   *
   * This will update the ghost values in u, and then will call applySinglePatch for each patch
   *
   * @param u the left hand side
   * @param f the right hand side
   */
  void
  apply(const Vector<D>& u, Vector<D>& f) const override;

  /**
   * @brief Get the Domain object associated with this PatchOperator
   */
  const Domain<D>&
  getDomain() const;

  /**
   * @brief Get the GhostFiller object associated with this PatchOperator
   */
  const GhostFiller<D>&
  getGhostFiller() const;
};
} // namespace ThunderEgg

// EXPLICIT INSTANTIATIONS

extern template class ThunderEgg::PatchOperator<2>;
extern template class ThunderEgg::PatchOperator<3>;
#endif
