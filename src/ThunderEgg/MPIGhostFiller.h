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

#ifndef THUNDEREGG_MPIGHOSTFILLER_H
#define THUNDEREGG_MPIGHOSTFILLER_H
/**
 * @file
 *
 * @brief MPIGhostFiller class
 */

#include <ThunderEgg/DimensionalArray.h>
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/GhostFillingType.h>

namespace ThunderEgg {

/**
 * @brief MPIGhostFiller implimentation
 *
 * @tparam D the number of Cartesian dimensions
 */
template<int D>
class MPIGhostFillerImpl;

/**
 * @brief Parallel ghostfiller implimented with MPI
 *
 * There are three functions that have to be overridden in derived classes.
 * fillGhostCellsForNbrPatch, fillGhostCellsForEdgeNbrPatch, fillGhostCellsForCornerNbrPatch, and
 * fillGhostCellsForLocalPatch
 *
 * @tparam D the number of Cartesian dimensions
 */
template<int D>
class MPIGhostFiller : public GhostFiller<D>
{
private:
  /**
   * @brief the fill data
   */
  std::shared_ptr<const MPIGhostFillerImpl<D>> implimentation;

public:
  /**
   * @brief Construct a new MPIGhostFiller object
   *
   * @param domain  the domain being used
   * @param fill_type  the number of side cases to address
   */
  MPIGhostFiller(const Domain<D>& domain, GhostFillingType fill_type);

  /**
   * @brief Fill the ghost cells for the neighboring patch
   *
   * @param pinfo the patch that ghost cells are being filled from
   * @param local_vew the view for patch that ghost cells are being filled from
   * @param nbr_view the view for the neighboring patch, where ghost cells are being
   * filled.
   * @param side the side that the neighboring patch is on
   * @param nbr_type the type of neighbor
   * @param orthant_on_coarse the orthant that the neighbors ghost cells lie on if the neighbor is
   * coarser
   */
  virtual void
  fillGhostCellsForNbrPatch(const PatchInfo<D>& pinfo,
                            const PatchView<const double, D>& local_view,
                            const PatchView<const double, D>& nbr_view,
                            Side<D> side,
                            NbrType nbr_type,
                            Orthant<D - 1> orthant_on_coarse) const = 0;
  /**
   * @brief Fill the edge ghost cells for the neighboring patch
   *
   * @param pinfo the patch that ghost cells are being filled from
   * @param local_view the view for patch that ghost cells are being filled from
   * @param nbr_view  the view for the neighboring patch, where ghost cells are being
   * filled.
   * @param edge the edge that the neighboring patch is on
   * @param nbr_type the type of neighbor
   * @param orthant_on_coarse the orthant that the neighbors ghost cells lie on if the neighbor is
   * coarser
   */
  virtual void
  fillGhostCellsForEdgeNbrPatch(const PatchInfo<D>& pinfo,
                                const PatchView<const double, D>& local_view,
                                const PatchView<const double, D>& nbr_view,
                                Edge edge,
                                NbrType nbr_type,
                                Orthant<1> orthant_on_coarse) const = 0;
  /**
   * @brief Fill the corner ghost cells for the neighboring patch
   *
   * @param pinfo the patch that ghost cells are being filled from
   * @param local_view the view for patch that ghost cells are being filled from
   * @param nbr_view  the view for the neighboring patch, where ghost cells are being
   * filled.
   * @param corner the edge that the neighboring patch is on
   * @param nbr_type the type of neighbor
   */
  virtual void
  fillGhostCellsForCornerNbrPatch(const PatchInfo<D>& pinfo,
                                  const PatchView<const double, D>& local_view,
                                  const PatchView<const double, D>& nbr_view,
                                  Corner<D> corner,
                                  NbrType nbr_type) const = 0;

  /**
   * @brief Perform any on this patches ghost cells.
   *
   * This may be necessary on some schemes because it needs data from the patch itself, not just
   * the neighboring patch
   *
   * @param pinfo the patch
   * @param view the view for the patch
   */
  virtual void
  fillGhostCellsForLocalPatch(const PatchInfo<D>& pinfo,
                              const PatchView<const double, D>& view) const = 0;

  /**
   * @brief Fill ghost cells on a vector
   *
   * @param u  the vector
   */
  void
  fillGhost(const Vector<D>& u) const override;

  /**
   * @brief Get the ghost filling type
   *
   * @return GhostFillingType the type
   */
  GhostFillingType
  getFillType() const;

  /**
   * @brief Get the domain that is being filled for
   *
   * @return std::shared_ptr<const Domain<D>>  the domain
   */
  const Domain<D>&
  getDomain() const;
};
extern template class MPIGhostFiller<2>;
extern template class MPIGhostFiller<3>;
} // namespace ThunderEgg
#endif
