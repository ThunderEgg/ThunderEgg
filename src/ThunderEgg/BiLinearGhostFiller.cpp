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

#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/RuntimeError.h>
namespace ThunderEgg {
BiLinearGhostFiller::BiLinearGhostFiller(const Domain<2>& domain, GhostFillingType fill_type)
  : MPIGhostFiller<2>(domain, fill_type)
{}
BiLinearGhostFiller*
BiLinearGhostFiller::clone() const
{
  return new BiLinearGhostFiller(*this);
}
namespace {
/**
 * @brief This is just a simple copy of values
 *
 * @param local_view the local pach
 * @param nbr_view the neighbor patch
 * @param side the side that the neighbor patch is on
 */
void
FillGhostForNormalNbr(const PatchView<const double, 2>& local_view,
                      const PatchView<const double, 2>& nbr_view,
                      const Side<2> side)
{
  View<const double, 2> local_slice = local_view.getSliceOn(side, { 0 });
  View<double, 2> nbr_ghosts = nbr_view.getGhostSliceOn(side.opposite(), { 0 });
  Loop::OverInteriorIndexes<2>(
    nbr_ghosts, [&](const std::array<int, 2>& coord) { nbr_ghosts[coord] = local_slice[coord]; });
}
/**
 * @brief Fill the ghost values for a coarse neighbor
 *
 * This is only part of the value. Values from the interior of the coarse neighbor will also be
 * needed.
 *
 * Those values are added in FillLocalGhostsForFineNbr
 *
 * @param local_view the local patch
 * @param nbr_view the neighbor patch
 * @param side the side that the neighbor patch is on
 * @param orthant the orthant of the coarser neighbor face that this patch lies on
 */
void
FillGhostForCoarseNbr(const PatchView<const double, 2>& local_view,
                      const PatchView<const double, 2>& nbr_view,
                      const Side<2> side,
                      const Orthant<1> orthant)
{
  int offset = 0;
  if (orthant == Orthant<1>::upper()) {
    offset = local_view.getEnd()[!side.getAxisIndex()] + 1;
  }
  View<const double, 2> local_slice = local_view.getSliceOn(side, { 0 });
  View<double, 2> nbr_ghosts = nbr_view.getGhostSliceOn(side.opposite(), { 0 });
  Loop::OverInteriorIndexes<2>(nbr_ghosts, [&](const std::array<int, 2>& coord) {
    nbr_ghosts[{ (coord[0] + offset) / 2, coord[1] }] += 2.0 / 3.0 * local_slice[coord];
  });
}
/**
 * @brief Fill the ghost values for a fine neighbor
 *
 * This is only part of the value. Values from the interior of the fine neighbor will also be
 * needed.
 *
 * Those values are added in FillLocalGhostsForCoarseNbr
 *
 * @param local_view the local patch
 * @param nbr_view the neighbor patch
 * @param side the side that the neighbor patch is on
 * @param orthant the orthant of the this patches face that the finer neighbor patch lies on
 */
void
FillGhostForFineNbr(const PatchView<const double, 2>& local_view,
                    const PatchView<const double, 2>& nbr_view,
                    const Side<2> side,
                    const Orthant<1> orthant)
{
  int offset = 0;
  if (orthant == Orthant<1>::upper()) {
    offset = local_view.getEnd()[!side.getAxisIndex()] + 1;
  }
  View<const double, 2> local_slice = local_view.getSliceOn(side, { 0 });
  View<double, 2> nbr_ghosts = nbr_view.getGhostSliceOn(side.opposite(), { 0 });
  Loop::OverInteriorIndexes<2>(nbr_ghosts, [&](const std::array<int, 2>& coord) {
    nbr_ghosts[coord] += 2.0 / 3.0 * local_slice[{ (coord[0] + offset) / 2, coord[1] }];
  });
}
/**
 * @brief Fill in the extra information needed for this patches ghost cells when there is a coarser
 * neighbor
 *
 * This is only part of the value.
 *
 * Those values compliment FillGhostForFineNbr
 *
 * @param pinfo the patchinfo object
 * @param view the patch data
 * @param side the side that the neighbor patch is on
 */
void
FillLocalGhostsForCoarseNbr(const PatchInfo<2>& pinfo,
                            const PatchView<const double, 2>& view,
                            const Side<2> side)
{
  View<const double, 2> local_slice = view.getSliceOn(side, { 0 });
  View<double, 2> local_ghosts = view.getGhostSliceOn(side, { 0 });
  int offset = 0;
  if (pinfo.getCoarseNbrInfo(side).orth_on_coarse == Orthant<1>::upper()) {
    offset = view.getEnd()[!side.getAxisIndex()] + 1;
  }
  Loop::OverInteriorIndexes<2>(local_ghosts, [&](const std::array<int, 2>& coord) {
    local_ghosts[coord] += 2.0 / 3.0 * local_slice[coord];
    if ((coord[0] + offset) % 2 == 0) {
      local_ghosts[{ coord[0] + 1, coord[1] }] += -1.0 / 3.0 * local_slice[coord];
    } else {
      local_ghosts[{ coord[0] - 1, coord[1] }] += -1.0 / 3.0 * local_slice[coord];
    }
  });
}
/**
 * @brief Fill in the extra information needed for this patches ghost cells when there is a finer
 * neighbor
 *
 * This is only part of the value.
 *
 * Those values compliment FillGhostForCoarseNbr
 *
 * @param view the patch data
 * @param side the side that the neighbor patch is on
 */
void
FillLocalGhostsForFineNbr(const PatchView<const double, 2>& view, const Side<2> side)
{
  View<const double, 2> local_slice = view.getSliceOn(side, { 0 });
  View<double, 2> local_ghosts = view.getGhostSliceOn(side, { 0 });
  Loop::OverInteriorIndexes<2>(local_ghosts, [&](const std::array<int, 2>& coord) {
    local_ghosts[coord] += -1.0 / 3.0 * local_slice[coord];
  });
}
/**
 * @brief This is just a simple copy of values
 *
 * @param local_view the local pach
 * @param nbr_view the neighbor patch
 * @param corner the corner that the neighbor patch is on
 */
void
FillGhostForCornerNormalNbr(const PatchView<const double, 2>& local_view,
                            const PatchView<const double, 2>& nbr_view,
                            Corner<2> corner)
{
  View<const double, 1> local_slice = local_view.getSliceOn(corner, { 0, 0 });
  View<double, 1> nbr_ghost = nbr_view.getGhostSliceOn(corner.opposite(), { 0, 0 });
  for (int c = local_slice.getStart()[0]; c <= local_slice.getEnd()[0]; c++) {
    nbr_ghost[{ c }] = local_slice[{ c }];
  }
}
/**
 * @brief Fill the ghost values for a coarse neighbor
 *
 * This is only part of the value. Values from the interior of the coarse neighbor will also be
 * needed.
 *
 * Those values are added in FillLocalGhostsForCornerFineNbr
 *
 * @param local_view the local patch
 * @param nbr_view the neighbor patch
 * @param corner the corner that the neighbor patch is on
 */
void
FillGhostForCornerCoarseNbr(const PatchView<const double, 2>& local_view,
                            const PatchView<const double, 2>& nbr_view,
                            Corner<2> corner)
{
  View<const double, 1> local_slice = local_view.getSliceOn(corner, { 0, 0 });
  View<double, 1> nbr_ghost = nbr_view.getGhostSliceOn(corner.opposite(), { 0, 0 });
  for (int c = local_slice.getStart()[0]; c <= local_slice.getEnd()[0]; c++) {
    nbr_ghost[{ c }] += 4.0 * local_slice[{ c }] / 3.0;
  }
}
/**
 * @brief Fill the ghost values for a fine neighbor
 *
 * This is only part of the value. Values from the interior of the fine neighbor will also be
 * needed.
 *
 * Those values are added in FillLocalGhostsForCornerCoarseNbr
 *
 * @param local_view the local patch
 * @param nbr_view the neighbor patch
 * @param corner the corner that the neighbor patch is on
 */
void
FillGhostForCornerFineNbr(const PatchView<const double, 2>& local_view,
                          const PatchView<const double, 2>& nbr_view,
                          Corner<2> corner)
{
  View<const double, 1> local_slice = local_view.getSliceOn(corner, { 0, 0 });
  View<double, 1> nbr_ghost = nbr_view.getGhostSliceOn(corner.opposite(), { 0, 0 });
  for (int c = local_slice.getStart()[0]; c <= local_slice.getEnd()[0]; c++) {
    nbr_ghost[{ c }] += 2.0 * local_slice[{ c }] / 3.0;
  }
}
/**
 * @brief Fill in the extra information needed for this patches ghost cells when there is a coarser
 * neighbor
 *
 * This is only part of the value. Values from the interior of the coarse neighbor will also be
 * needed.
 *
 * Those values compliment FillGhostForCornerFineNbr
 *
 * @param view the patch data
 * @param corner the corner that the neighbor patch is on
 */
void
FillLocalGhostsForCornerCoarseNbr(const PatchView<const double, 2>& view, Corner<2> corner)
{
  View<const double, 1> local_slice = view.getSliceOn(corner, { 0, 0 });
  View<double, 1> local_ghosts = view.getGhostSliceOn(corner, { 0, 0 });
  for (int c = local_slice.getStart()[0]; c <= local_slice.getEnd()[0]; c++) {
    local_ghosts[{ c }] += local_slice[{ c }] / 3.0;
  }
}
/**
 * @brief Fill in the extra information needed for this patches ghost cells when there is a finer
 * neighbor
 *
 * This is only part of the value.
 *
 * Those values compliment FillGhostForCornerCoarseNbr
 *
 * @param view the patch data
 * @param corner the corner that the neighbor patch is on
 */
void
FillLocalGhostsForCornerFineNbr(const PatchView<const double, 2>& view, Corner<2> corner)
{
  View<const double, 1> local_slice = view.getSliceOn(corner, { 0, 0 });
  View<double, 1> local_ghosts = view.getGhostSliceOn(corner, { 0, 0 });
  for (int c = local_slice.getStart()[0]; c <= local_slice.getEnd()[0]; c++) {
    local_ghosts[{ c }] += -local_slice[{ c }] / 3.0;
  }
}
/**
 * @brief Add in extra information needed from the local patch on the sides
 *
 * @param pinfo the pinfo object
 * @param view the patch data
 */
void
FillLocalGhostCellsOnSides(const PatchInfo<2>& pinfo, const PatchView<const double, 2>& view)
{
  for (Side<2> side : Side<2>::getValues()) {
    if (pinfo.hasNbr(side)) {
      switch (pinfo.getNbrType(side)) {
        case NbrType::Normal:
          // nothing needs to be done
          break;
        case NbrType::Coarse:
          FillLocalGhostsForCoarseNbr(pinfo, view, side);
          break;
        case NbrType::Fine:
          FillLocalGhostsForFineNbr(view, side);
          break;
        default:
          throw RuntimeError("Unsupported Nbr Type");
      }
    }
  }
}
/**
 * @brief Add in extra information needed from the local patch on the corner
 *
 * @param pinfo the pinfo object
 * @param view the patch view
 */
void
FillLocalGhostCellsOnCorners(const PatchInfo<2>& pinfo, const PatchView<const double, 2>& view)
{
  for (Corner<2> corner : Corner<2>::getValues()) {
    if (pinfo.hasNbr(corner)) {
      switch (pinfo.getNbrType(corner)) {
        case NbrType::Normal:
          // nothing needs to be done
          break;
        case NbrType::Coarse:
          FillLocalGhostsForCornerCoarseNbr(view, corner);
          break;
        case NbrType::Fine:
          FillLocalGhostsForCornerFineNbr(view, corner);
          break;
        default:
          throw RuntimeError("Unsupported Nbr Type");
      }
    }
  }
}
} // namespace
void
BiLinearGhostFiller::fillGhostCellsForNbrPatch(const PatchInfo<2>& pinfo,
                                               const PatchView<const double, 2>& local_view,
                                               const PatchView<const double, 2>& nbr_view,
                                               Side<2> side,
                                               NbrType nbr_type,
                                               Orthant<1> orthant_on_coarse) const
{
  switch (nbr_type) {
    case NbrType::Normal:
      FillGhostForNormalNbr(local_view, nbr_view, side);
      break;
    case NbrType::Coarse:
      FillGhostForCoarseNbr(local_view, nbr_view, side, orthant_on_coarse);
      break;
    case NbrType::Fine:
      FillGhostForFineNbr(local_view, nbr_view, side, orthant_on_coarse);
      break;
    default:
      throw RuntimeError("Unsupported Nbr Type");
  }
}
void
BiLinearGhostFiller::fillGhostCellsForEdgeNbrPatch(const PatchInfo<2>& pinfo,
                                                   const PatchView<const double, 2>& local_view,
                                                   const PatchView<const double, 2>& nbr_view,
                                                   Edge edge,
                                                   NbrType nbr_type,
                                                   Orthant<1> orthant_on_coarse) const
{
  // 2D, edges not needed
}

void
BiLinearGhostFiller::fillGhostCellsForCornerNbrPatch(const PatchInfo<2>& pinfo,
                                                     const PatchView<const double, 2>& local_view,
                                                     const PatchView<const double, 2>& nbr_view,
                                                     Corner<2> corner,
                                                     NbrType nbr_type) const
{
  switch (nbr_type) {
    case NbrType::Normal:
      FillGhostForCornerNormalNbr(local_view, nbr_view, corner);
      break;
    case NbrType::Coarse:
      FillGhostForCornerCoarseNbr(local_view, nbr_view, corner);
      break;
    case NbrType::Fine:
      FillGhostForCornerFineNbr(local_view, nbr_view, corner);
      break;
    default:
      throw RuntimeError("Unsupported Nbr Type");
  }
}
void
BiLinearGhostFiller::fillGhostCellsForLocalPatch(const PatchInfo<2>& pinfo,
                                                 const PatchView<const double, 2>& view) const
{
  switch (this->getFillType()) {
    case GhostFillingType::Corners: // Fill corners and faces
      FillLocalGhostCellsOnCorners(pinfo, view);
      [[fallthrough]];
    case GhostFillingType::Faces:
      FillLocalGhostCellsOnSides(pinfo, view);
      break;
    default:
      throw RuntimeError("Unsupported GhostFillingType");
  }
}
} // namespace ThunderEgg