/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020      Scott Aiton
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

#include <ThunderEgg/GMG/DirectInterpolator.h>
#include <ThunderEgg/GMG/InterLevelComm.h>

namespace ThunderEgg::GMG {
template<int D>
  requires is_supported_dimension<D>
class DirectInterpolator<D>::Implimentation
{
public:
  void
  interpolatePatches(
    const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>>& patches,
    const Vector<D>& coarser_vector,
    Vector<D>& finer_vector) const
  {
    for (auto pair : patches) {
      const PatchInfo<D>& pinfo = pair.second.get();
      PatchView<const double, D> coarse_view = coarser_vector.getPatchView(pair.first);
      PatchView<double, D> fine_view = finer_vector.getPatchView(pinfo.local_index);

      if (pinfo.hasCoarseParent()) {
        Orthant<D> orth = pinfo.orth_on_parent;
        std::array<int, D> starts;
        for (size_t i = 0; i < D; i++) {
          starts[i] = orth.isLowerOnAxis(i) ? 0 : (coarse_view.getEnd()[i] + 1);
        }

        Loop::OverInteriorIndexes<D + 1>(fine_view, [&](const std::array<int, D + 1>& coord) {
          std::array<int, D + 1> coarse_coord;
          for (size_t x = 0; x < D; x++) {
            coarse_coord[x] = (coord[x] + starts[x]) / 2;
          }
          coarse_coord[D] = coord[D];
          fine_view[coord] += coarse_view[coarse_coord];
        });
      } else {
        Loop::OverInteriorIndexes<D + 1>(fine_view, [&](const std::array<int, D + 1>& coord) {
          fine_view[coord] += coarse_view[coord];
        });
      }
    }
  }
};

template<int D>
  requires is_supported_dimension<D>
DirectInterpolator<D>::DirectInterpolator(const Domain<D>& coarse_domain,
                                          const Domain<D>& fine_domain)
  : MPIInterpolator<D>(coarse_domain, fine_domain)
  , implimentation(new Implimentation())
{
}

template<int D>
  requires is_supported_dimension<D>
DirectInterpolator<D>::~DirectInterpolator() = default;

template<int D>
  requires is_supported_dimension<D>
DirectInterpolator<D>::DirectInterpolator(const DirectInterpolator& other)
  : MPIInterpolator<D>(other)
  , implimentation(new Implimentation(*other.implimentation))
{
}

template<int D>
  requires is_supported_dimension<D>
DirectInterpolator<D>::DirectInterpolator(DirectInterpolator&& other) = default;

template<int D>
  requires is_supported_dimension<D>
DirectInterpolator<D>&
DirectInterpolator<D>::operator=(const DirectInterpolator<D>& other)
{
  implimentation.reset(new Implimentation(*other.implimentation));
  return *this;
}

template<int D>
  requires is_supported_dimension<D>
DirectInterpolator<D>&
DirectInterpolator<D>::operator=(DirectInterpolator<D>&& other) = default;

template<int D>
  requires is_supported_dimension<D>
DirectInterpolator<D>*
DirectInterpolator<D>::clone() const
{
  return new DirectInterpolator<D>(*this);
}

template<int D>
  requires is_supported_dimension<D>
void
DirectInterpolator<D>::interpolatePatches(
  const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>>& patches,
  const Vector<D>& coarser_vector,
  Vector<D>& finer_vector) const
{
  implimentation->interpolatePatches(patches, coarser_vector, finer_vector);
}
template class DirectInterpolator<2>;
template class DirectInterpolator<3>;
} // namespace ThunderEgg::GMG