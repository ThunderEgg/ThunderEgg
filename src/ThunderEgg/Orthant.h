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

#ifndef THUNDEREGG_ORTHANT_H
#define THUNDEREGG_ORTHANT_H
/**
 * @file
 *
 * @brief Orthant class
 */
#include <ThunderEgg/Face.h>
#include <array>
#include <iostream>
#include <numeric>
#include <vector>

namespace ThunderEgg {

template<int D>
constexpr bool is_supported_side_dimension = D >= 1 && D <= 3;

template<int D>
constexpr bool is_supported_orthant_dimension = D >= 0 && D <= 3;

/**
 * @brief An enum-style class that represents the octants of a cube.
 *
 * The octants are named in the following way:
 *
 * bse means Bottom-South-West, which is in the corner of where the bottom, south, and west sides
 * meet.
 */
template<int D>
  requires is_supported_orthant_dimension<D>
class Orthant
{
private:
  /**
   * @brief the value of the enum.
   */
  unsigned char val = num_orthants;

public:
  static constexpr size_t num_orthants = 1 << D;

  /**
   * @brief Create new Orthant<D> with given value.
   *
   * @param val_in the value
   */
  explicit Orthant(const unsigned char val_in);

  /**
   * @brief Default constructor that initializes the value to null().
   */
  Orthant();

  /**
   * @brief null value
   */
  static Orthant<D>
  null();

  /**
   * @brief Lower half of line.
   */
  static Orthant<1>
  lower();

  /**
   * @brief Upper half of line.
   */
  static Orthant<1>
  upper();

  /**
   * @brief South-West quadrant of square.
   */
  static Orthant<2>
  sw();

  /**
   * @brief South-East quadrant of square.
   */
  static Orthant<2>
  se();

  /**
   * @brief North-West quadrant of square.
   */
  static Orthant<2>
  nw();

  /**
   * @brief North-East quadrant of square.
   */
  static Orthant<2>
  ne();

  /**
   * @brief Bottom-South-West octant of cube.
   */
  static Orthant<3>
  bsw();

  /**
   * @brief Bottom-South-East octant of cube.
   */
  static Orthant<3>
  bse();

  /**
   * @brief Bottom-North-West octant of cube.
   */
  static Orthant<3>
  bnw();

  /**
   * @brief Bottom-North-East octant of cube.
   */
  static Orthant<3>
  bne();

  /**
   * @brief Top-South-West octant of cube.
   */
  static Orthant<3>
  tsw();

  /**
   * @brief Top-South-East octant of cube.
   */
  static Orthant<3>
  tse();

  /**
   * @brief Top-North-West octant of cube.
   */
  static Orthant<3>
  tnw();

  /**
   * @brief Top-North-East octant of cube.
   */
  static Orthant<3>
  tne();

  /**
   * @brief Range class for Orthant
   *
   * It provides the begin and end fuctions for iterator loops
   */
  class Range
  {
  public:
    /**
     * @brief Input iterator for Orthant values
     */
    class Iterator : public std::input_iterator_tag
    {
    private:
      /**
       * @brief The current side
       */
      Orthant<D> o;

    public:
      /**
       * @brief Construct a new Iterator object with the given Orthant value
       *
       * @param o_in the orthant
       */
      explicit Iterator(Orthant<D> o_in);

      /**
       * @brief Increment the Orthant value
       *
       * @return const Orthant<D>& the resulting value
       */
      const Orthant<D>&
      operator++();

      /**
       * @brief Get a reference to the Orthant object
       *
       * @return const Orthant<D>&  the reference
       */
      const Orthant<D>&
      operator*() const;

      /**
       * @brief Get a pointer to the Orthant object
       *
       * @return const Orthant<D>* the pointer
       */
      const Orthant<D>*
      operator->() const;

      /**
       * @brief Check the iterators reference the same value
       *
       * @param b the other iterator
       * @return true if the same
       * @return false if different
       */
      bool
      operator==(const Iterator& b) const;

      /**
       * @brief Check the iterators don't reference the same value
       *
       * @param b the other iterator
       * @return true if different
       * @return false if the same
       */
      bool
      operator!=(const Iterator& b) const;
    };

    /**
     * @brief Returns an iterator with the lowest Orthant value
     *
     * @return Iterator the iterator
     */
    Iterator
    begin();

    /**
     * @brief Returns an iterator with Orthant<D>::null()
     *
     * @return Iterator the iterator
     */
    Iterator
    end();
  };

  /**
   * @brief Get a range of values that can be iterated over
   *
   * @return Range the range of values
   */
  static Range
  getValues();

  /**
   * @brief Get the integer value of the octant.
   *
   * @return The integer value.
   */
  size_t
  getIndex() const;

  /**
   * @brief Return the octant that neighbors this octant on a particular side.
   *
   * @param s the side of the octant that you want the neighbor of.
   *
   * @return  The octant that neighbors on that side.
   */
  Orthant<D>
  getNbrOnSide(Side<D> s) const;

  /**
   * @brief Get the sides of the octant that are on the interior of the cube.
   *
   * @return The sides of the octant that are on the interior of the cube.
   */
  std::array<Side<D>, D>
  getInteriorSides() const;

  /**
   * @brief Get the sides of the octant that are on the exterior of the cube.
   *
   * @return The sides of the octant that are on the exterior of the cube.
   */
  std::array<Side<D>, D>
  getExteriorSides() const;

  /**
   * @brief Return whether or not the octant lies on a particular side of a cube.
   *
   * @param s the side of the cube.j
   *
   * @return Whether or not it lies on that side.
   */
  bool
  isOnSide(Side<D> s) const;

  bool
  isHigherOnAxis(size_t axis) const;

  bool
  isLowerOnAxis(size_t axis) const;

  /**
   * @brief From the point of view of an axis, get orthant that this orthant lies on in the D-1
   * dimension
   *
   * @param axis the axis
   * @return Orthant<D - 1> the resulting orthant
   */
  template<int Dm1 = D - 1>
    requires(Dm1 == D - 1) && (D - 1 >= 0)
  Orthant<Dm1> collapseOnAxis(size_t axis) const;

  /**
   * @brief Get an array of all Orthant<D> values that lie on a particular side of the cube.
   *
   * When the two axis that the side lies on are arranged in the following way, the octants are
   * returned in the following order:
   *
   *   ^
   *   |
   *   |  2  |  3
   *   |-----+-----
   *   |  0  |  1
   *   +----------->
   *
   * @return The array.
   */

  static std::array<Orthant, num_orthants / 2>
  getValuesOnSide(Side<D> s);

  /**
   * @brief Equals operator.
   *
   * @param other The other octant.
   *
   * @return Whether or not the value of this octant equals the value other octant.
   */
  bool
  operator==(const Orthant<D>& other) const;

  /**
   * @brief Not Equals operator.
   *
   * @param other The other octant.
   *
   * @return Whether or not the value of this octant is not equal the value other octant.
   */
  bool
  operator!=(const Orthant<D>& other) const;

  /**
   * @brief Less Tan operator.
   *
   * @param other The other octant.
   *
   * @return Whether or not the value of this octant is less than the value other octant.
   */
  bool
  operator<(const Orthant<D>& other) const;
};

/**
 * @brief ostream operator that prints a string representation of Orthant<0> enum.
 *
 * For example, Orthant<0>::null() will print out "Orthant<0>::null()".
 *
 * @param os the ostream
 * @param o the orthant
 *
 * @return  the ostream
 */
std::ostream&
operator<<(std::ostream& os, const Orthant<0>& o);

/**
 * @brief ostream operator that prints a string representation of Orthant<1> enum.
 *
 * For example, Orthant<1>::lower() will print out "Orthant<1>::lower()".
 *
 * @param os the ostream
 * @param o the orthant
 *
 * @return  the ostream
 */
std::ostream&
operator<<(std::ostream& os, const Orthant<1>& o);

/**
 * @brief ostream operator that prints a string representation of quadrant enum.
 *
 * For example, Orthant<2>::sw() will print out "Orthant<2>::sw()".
 *
 * @param os the ostream
 * @param o the orthant
 *
 * @return  the ostream
 */
std::ostream&
operator<<(std::ostream& os, const Orthant<2>& o);

/**
 * @brief ostream operator that prints a string representation of orthant enum.
 *
 * For example, Orthant<3>::bsw() will print out "Orthant<3>::bsw()".
 *
 * @param os the ostream
 * @param o the orthant
 *
 * @return  the ostream
 */
std::ostream&
operator<<(std::ostream& os, const Orthant<3>& o);

void
to_json(tpl::nlohmann::json& j, const Orthant<0>& o);
void
to_json(tpl::nlohmann::json& j, const Orthant<1>& o);
void
to_json(tpl::nlohmann::json& j, const Orthant<2>& o);
void
to_json(tpl::nlohmann::json& j, const Orthant<3>& o);
void
from_json(const tpl::nlohmann::json& j, Orthant<0>& o);
void
from_json(const tpl::nlohmann::json& j, Orthant<1>& o);
void
from_json(const tpl::nlohmann::json& j, Orthant<2>& o);
void
from_json(const tpl::nlohmann::json& j, Orthant<3>& o);

// EXPLICIT INSTANTIATIONS

extern template class Orthant<0>;

extern template class Orthant<1>;

extern template class Orthant<2>;

extern template class Orthant<3>;

} // namespace ThunderEgg
#endif
