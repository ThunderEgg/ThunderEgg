/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2021      Scott Aiton
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

#ifndef THUNDEREGG_COMPONENTARRAY_H
#define THUNDEREGG_COMPONENTARRAY_H
#include <ThunderEgg/ComponentView.h>
/**
 * @file
 *
 * @brief ComponentArray class
 */
namespace ThunderEgg {
/**
 * @brief Array for acessing data of a patch. It supports variable striding
 *
 * @tparam D number of cartesian dimensions
 */
template<int D>
class ComponentArray
{
private:
  std::vector<double> vector;
  ComponentView<double, D> view;

public:
  /**
   * @brief Construct a new View object
   *
   * @param lengths the lengths in each direction
   * @param num_ghost_cells the number of ghost cells on each side of the patch
   */
  ComponentArray(const std::array<int, D>& lengths, int num_ghost_cells)
  {
    std::array<int, D> strides;
    strides[0] = 1;
    for (int i = 1; i < D; i++) {
      strides[i] = strides[i - 1] * (lengths[i - 1] + 2 * num_ghost_cells);
    }
    int size = strides[D - 1] * (lengths[D - 1] + 2 * num_ghost_cells);
    vector.resize(size);
    view = ComponentView<double, D>(vector.data(), strides, lengths, num_ghost_cells);
  }

  /**
   * @brief Copy constructor
   *
   * @param other the array to copy
   */
  ComponentArray(const ComponentArray<D>& other)
    : vector(other.vector)
    , view(vector.data(),
           other.getStrides(),
           other.getGhostStart(),
           other.getStart(),
           other.getEnd(),
           other.getGhostEnd())

  {}

  /**
   * @brief Copy assignment
   *
   * @param other the array to copy
   * @return PatchArray<D>&  this
   */
  ComponentArray<D>& operator=(const ComponentArray<D>& other)
  {
    vector = other.vector;
    view = ComponentView<double, D>(vector.data(),
                                    other.getStrides(),
                                    other.getGhostStart(),
                                    other.getStart(),
                                    other.getEnd(),
                                    other.getGhostEnd());
    return *this;
  }
  /**
   * @brief Get the slice on a given face
   *
   * @tparam M the dimension of the face
   * @param f the face
   * @param offset offset the offset of the value {0,..,0} is the non ghost cell touching the face.
   * {-1,..,-1} is the first ghost cell touching that face
   * @return View<M> a view to the slice on the face
   */
  template<int M>
  inline View<double, M> getSliceOn(Face<D, M> f, const std::array<int, D - M>& offset)
  {
    return view.template getSliceOn<M>(f, offset);
  }

  /**
   * @brief Get the slice on a given face
   *
   * @tparam M the dimension of the face
   * @param f the face
   * @param offset offset the offset of the value {0,..,0} is the non ghost cell touching the face.
   * {-1,..,-1} is the first ghost cell touching that face
   * @return ConstView<M> a view to the slice on the face
   */
  template<int M>
  inline View<const double, M> getSliceOn(Face<D, M> f, const std::array<int, D - M>& offset) const
  {
    return View<const double, M>(view.template getSliceOn<M>(f, offset));
  }

  /**
   * @brief Get the gosts slice on a given face
   *
   * @tparam M the dimension of the face
   * @param f the face
   * @param offset offset the offset of the value {0,..,0} first ghost cell slice touching that face
   * face
   * @return View<M> a view to the slice on the face
   */
  template<int M>
  inline View<double, M> getGhostSliceOn(Face<D, M> f,
                                         const std::array<size_t, D - M>& offset) const
  {
    return view.template getGhostSliceOn<M>(f, offset);
  }

  inline const double& operator[](const std::array<int, D>& coord) const { return view[coord]; }
  template<class... Types>
  inline const double& operator()(Types... args) const
  {
    return view(args...);
  }
  inline void set(const std::array<int, D>& coord, double value) const { view.set(coord, value); }
  inline double& operator[](const std::array<int, D>& coord) { return view[coord]; }
  template<class... Types>
  inline double& operator()(Types... args)
  {
    return view(args...);
  }
  inline void set(const std::array<int, D>& coord, double value) { view.set(coord, value); }

  /**
   * @brief Get the strides of the patch in each direction
   */
  inline const std::array<int, D>& getStrides() const { return view.getStrides(); }
  /**
   * @brief Get the coordinate of the first element
   */
  inline const std::array<int, D>& getStart() const { return view.getStart(); }
  /**
   * @brief Get the coordinate of the last element
   */
  inline const std::array<int, D>& getEnd() const { return view.getEnd(); }
  /**
   * @brief Get the coordinate of the first ghost cell element
   */
  inline const std::array<int, D>& getGhostStart() const { return view.getGhostStart(); }
  /**
   * @brief Get the coordinate of the last ghost cell element
   */
  inline const std::array<int, D>& getGhostEnd() const { return view.getGhostEnd(); }
};
extern template class ComponentArray<1>;
extern template class ComponentArray<2>;
extern template class ComponentArray<3>;
} // namespace ThunderEgg
#endif