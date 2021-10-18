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
class LinearRestrictor : public MPIRestrictor<D>
{
private:
  /**
   * @brief true if ghost values at boundaries should be extrapolated
   */
  bool extrapolate_boundary_ghosts;

  /**
   * @brief Extrapolate to the ghosts on the parent patch
   *
   * @param pinfo the patch
   * @param fine_view the finer patch
   * @param coarse_view the coarser patch
   */
  void extrapolateBoundaries(const PatchInfo<D>& pinfo,
                             const PatchView<const double, D>& fine_view,
                             const PatchView<double, D>& coarse_view) const
  {
    Orthant<D> orth = pinfo.orth_on_parent;
    std::array<int, D> starts;
    for (size_t i = 0; i < D; i++) {
      starts[i] = orth.isLowerOnAxis(i) ? 0 : coarse_view.getEnd()[i] + 1;
    }
    // extrapolate ghost values
    for (Side<D> s : pinfo.orth_on_parent.getExteriorSides()) {
      // if (!pinfo.hasNbr(s)) {
      View<const double, D> fine_ghost = fine_view.getSliceOn(s, { -1 });
      View<const double, D> fine_interior = fine_view.getSliceOn(s, { 0 });
      View<double, D> coarse_ghost = coarse_view.getSliceOn(s, { -1 });
      Loop::OverInteriorIndexes<D>(fine_ghost, [&](const std::array<int, D>& coord) {
        std::array<int, D> coarse_coord;
        for (size_t x = 0; x < s.getAxisIndex(); x++) {
          coarse_coord[x] = (coord[x] + starts[x]) / 2;
        }
        for (size_t x = s.getAxisIndex() + 1; x < D; x++) {
          coarse_coord[x - 1] = (coord[x - 1] + starts[x]) / 2;
        }
        coarse_coord[D - 1] = coord[D - 1];
        coarse_ghost[coarse_coord] += (3 * fine_ghost[coord] - fine_interior[coord]) / (1 << D);
      });
      //}
    }
    if constexpr (D >= 2) {
      Corner<D> c(pinfo.orth_on_parent.getIndex());
      if (!pinfo.hasNbr(c)) {
        std::array<int, D> neg_one;
        neg_one.fill(-1);
        std::array<int, D> zero;
        zero.fill(0);
        View<const double, 1> fine_ghost = fine_view.getSliceOn(c, neg_one);
        View<const double, 1> fine_interior = fine_view.getSliceOn(c, zero);
        View<double, 1> coarse_ghost = coarse_view.getSliceOn(c, neg_one);
        Loop::OverInteriorIndexes<1>(fine_ghost, [&](const std::array<int, 1>& coord) {
          coarse_ghost[coord] += (3 * fine_ghost[coord] - fine_interior[coord]);
        });
      }
    }
  }
  /**
   * @brief Restrict to a coarser patch
   *
   * @param pinfo the patch
   * @param parent_index the index of the coarser patch
   * @param finer_vector the finer vector
   * @param coarser_vector the coarser vector
   */
  void restrictToCoarserParent(const PatchInfo<D>& pinfo,
                               int parent_index,
                               const Vector<D>& finer_vector,
                               Vector<D>& coarser_vector) const
  {
    PatchView<double, D> coarse_view = coarser_vector.getPatchView(parent_index);
    PatchView<const double, D> fine_view = finer_vector.getPatchView(pinfo.local_index);
    // get starting index in coarser patch
    Orthant<D> orth = pinfo.orth_on_parent;
    std::array<int, D> starts;
    for (size_t i = 0; i < D; i++) {
      starts[i] = orth.isLowerOnAxis(i) ? 0 : (coarse_view.getEnd()[i] + 1);
    }

    // interpolate interior values
    Loop::OverInteriorIndexes<D + 1>(fine_view, [&](const std::array<int, D + 1>& coord) {
      std::array<int, D + 1> coarse_coord;
      for (size_t x = 0; x < D; x++) {
        coarse_coord[x] = (coord[x] + starts[x]) / 2;
      }
      coarse_coord[D] = coord[D];
      coarse_view[coarse_coord] += fine_view[coord] / (1 << D);
    });

    if (extrapolate_boundary_ghosts) {
      extrapolateBoundaries(pinfo, fine_view, coarse_view);
    }
  }
  /**
   * @brief Copy to a parent patch
   *
   * @param pinfo the patch
   * @param parent_index the index of the coarser patch
   * @param finer_vector the finer vector
   * @param coarser_vector the coarser vector
   */

  void copyToParent(const PatchInfo<D>& pinfo,
                    int parent_index,
                    const Vector<D>& finer_vector,
                    Vector<D>& coarser_vector) const
  {
    PatchView<double, D> coarse_view = coarser_vector.getPatchView(parent_index);
    PatchView<const double, D> fine_view = finer_vector.getPatchView(pinfo.local_index);
    // just copy the values
    if (extrapolate_boundary_ghosts) {
      Loop::OverAllIndexes<D + 1>(fine_view, [&](const std::array<int, D + 1>& coord) {
        coarse_view[coord] += fine_view[coord];
      });
    } else {
      Loop::OverInteriorIndexes<D + 1>(fine_view, [&](const std::array<int, D + 1>& coord) {
        coarse_view[coord] += fine_view[coord];
      });
    }
  }

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
                   bool extrapolate_boundary_ghosts = false)
    : MPIRestrictor<D>(coarse_domain, fine_domain)
    , extrapolate_boundary_ghosts(extrapolate_boundary_ghosts)
  {}

  /**
   * @brief Clone this restrictor
   *
   * @return LinearRestrictor<D>* a newly allocated copy of this restrictor
   */
  LinearRestrictor<D>* clone() const override { return new LinearRestrictor<D>(*this); }
  void restrictPatches(
    const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>>& patches,
    const Vector<D>& finer_vector,
    Vector<D>& coarser_vector) const override
  {
    for (const auto& pair : patches) {
      if (pair.second.get().hasCoarseParent()) {
        restrictToCoarserParent(pair.second.get(), pair.first, finer_vector, coarser_vector);
      } else {
        copyToParent(pair.second.get(), pair.first, finer_vector, coarser_vector);
      }
    }
  }
};
} // namespace ThunderEgg::GMG
// explicit instantiation
extern template class ThunderEgg::GMG::LinearRestrictor<2>;
extern template class ThunderEgg::GMG::LinearRestrictor<3>;
#endif