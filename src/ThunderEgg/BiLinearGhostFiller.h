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

#ifndef THUNDEREGG_BILINEARGHOSTFILLER_H
#define THUNDEREGG_BILINEARGHOSTFILLER_H
/**
 * @file
 *
 * @brief BiLinearGhostFiller class
 */

#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/MPIGhostFiller.h>
namespace ThunderEgg {
/**
 * @brief Exchanges ghost cells on patches, uses a BiLinear interpolation scheme for refinement
 * boundaries
 */
class BiLinearGhostFiller : public MPIGhostFiller<2>
{
public:
  /**
   * @brief Construct a new BiLinearGhostFiller object
   *
   * @param domain the domain to fill ghosts for
   * @param fill_type the GhostFillingType
   */
  BiLinearGhostFiller(const Domain<2>& domain, GhostFillingType fill_type);
  /**
   * @brief Clone this BiLinearGhostFiller
   *
   * @return BiLinearGhostFiller* a newly allocated copy of this BiLinearGhostFiller
   */
  BiLinearGhostFiller* clone() const override;

  void fillGhostCellsForNbrPatch(const PatchInfo<2>& pinfo,
                                 const PatchView<const double, 2>& local_view,
                                 const PatchView<const double, 2>& nbr_view,
                                 Side<2> sides,
                                 NbrType nbr_type,
                                 Orthant<1> orthant_on_coarse) const override;

  void fillGhostCellsForEdgeNbrPatch(const PatchInfo<2>& pinfo,
                                     const PatchView<const double, 2>& local_view,
                                     const PatchView<const double, 2>& nbr_view,
                                     Edge edge,
                                     NbrType nbr_type,
                                     Orthant<1> orthant_on_coarse) const override;

  void fillGhostCellsForCornerNbrPatch(const PatchInfo<2>& pinfo,
                                       const PatchView<const double, 2>& local_view,
                                       const PatchView<const double, 2>& nbr_view,
                                       Corner<2> corner,
                                       NbrType nbr_type) const override;

  void fillGhostCellsForLocalPatch(const PatchInfo<2>& pinfo,
                                   const PatchView<const double, 2>& view) const override;
};
} // namespace ThunderEgg
#endif
