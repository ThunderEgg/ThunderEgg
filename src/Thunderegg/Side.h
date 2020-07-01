/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
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

#ifndef THUNDEREGG_SIDE_H
#define THUNDEREGG_SIDE_H
#include <array>
#include <iostream>
#include <numeric>
#include <vector>
namespace Thunderegg
{
/**
 * @brief An enum-style class that represents the sides of a cube.
 *
 * The sides of the cube are named in the following way:
 *
 * Orthogonal to axis | Lower on axis | Higher on axis
 * ------------------ | ------------- | --------------
 *  x-axis            | west          | east
 *  y-axis            | south         | north
 *  z-axis            | bottom        | top
 *
 */
template <size_t D> class Side
{
	private:
	/**
	 * @brief the value of the enum
	 */
	int val = -1;

	public:
	// enum definitions
	static constexpr int west   = 0b000;
	static constexpr int east   = 0b001;
	static constexpr int south  = 0b010;
	static constexpr int north  = 0b011;
	static constexpr int bottom = 0b100;
	static constexpr int top    = 0b101;

	static constexpr int num_sides = 2 * D;

	/**
	 * @brief Default constructor that initializes the value to -1
	 */
	Side() = default;
	/**
	 * @brief Initialize new Side with given value
	 *
	 * @param val the value
	 */
	Side(const int val)
	{
		this->val = val;
	}
	/**
	 * @brief Get the integer value of the side.
	 *
	 * @return The value of the side.
	 */
	int toInt() const
	{
		return val;
	}
	operator int() const
	{
		return val;
	}
	operator size_t() const
	{
		return val;
	}
	/**
	 * @brief Get an array of all side values, in increasing order.
	 *
	 * @return The array.
	 */
	static std::array<Side<D>, num_sides> getValues();
	/**
	 * @brief Return whether or not the side of the cube is lower on the axis that is orthogonal to
	 * it.
	 *
	 * For example: The x-axis is orthogonal to both the west and east sides. Since west is lower on
	 * the axis, west will return true and east will return false.
	 *
	 * @return Whether or not it is lower on the axis.
	 */
	inline bool isLowerOnAxis() const
	{
		// is least-significant bit set?
		return !(val & 0x1);
	}
	/**
	 * @brief Return whether or not the side of the cube is higher on the axis that is orthogonal to
	 * it.
	 *
	 * For example: The x-axis is orthogonal to both the west and east sides. Since west is lower on
	 * the axis, west will return false and east will return true.
	 *
	 * @return Whether or not it is higher on the axis.
	 */
	inline bool isHigherOnAxis() const
	{
		// is least-significant bit set?
		return (val & 0x1);
	}
	/**
	 * @brief Return the axis that the side lies on.
	 */
	inline int axis() const
	{
		return val / 2;
	}
	/**
	 * @brief Return the opposite side of the cube.
	 *
	 * For example: the opposite of east is west.
	 *
	 * @return The opposite side.
	 */
	Side opposite() const;
	/**
	 * @brief Compare the enum values.
	 *
	 * @param other The other side.
	 *
	 * @return Whether or not the value of this side is lower than the other side.
	 */
	bool operator<(const Side &other) const
	{
		return val < other.val;
	}
	/**
	 * @brief Equals operator.
	 *
	 * @param other The other side.
	 *
	 * @return Whether or not the value of this side equals the value other side.
	 */
	bool operator==(const Side &other) const
	{
		return val == other.val;
	}
	/**
	 * @brief Equals operator.
	 *
	 * @param other The other value.
	 *
	 * @return Whether or not the value of this side equals the value of the integer.
	 */
	bool operator==(const int &other) const
	{
		return val == other;
	}
	/**
	 * @brief Not Equals operator.
	 *
	 * @param other The other side.
	 *
	 * @return Whether or not the value of this side equals the value other side.
	 */
	bool operator!=(const Side &other) const
	{
		return val != other.val;
	}
	/**
	 * @brief Not Equals operator.
	 *
	 * @param other The other value.
	 *
	 * @return Whether or not the value of this side equals the value of the integer.
	 */
	bool operator!=(const int &other) const
	{
		return val != other;
	}
};
template <size_t D> inline std::array<Side<D>, Side<D>::num_sides> Side<D>::getValues()
{
	std::array<Side<D>, Side<D>::num_sides> retval;
	std::iota(retval.begin(), retval.end(), 0);
	return retval;
}
template <size_t D> inline Side<D> Side<D>::opposite() const
{
	Side<D> retval = *this;
	retval.val ^= 0x1;
	return retval;
}
template <size_t D> constexpr int Side<D>::west;
template <size_t D> constexpr int Side<D>::east;
template <size_t D> constexpr int Side<D>::south;
template <size_t D> constexpr int Side<D>::north;
template <size_t D> constexpr int Side<D>::bottom;
template <size_t D> constexpr int Side<D>::top;
template <size_t D> constexpr int Side<D>::num_sides;
/**
 * @brief ostream operator that prints a string representation of side enum.
 *
 * For example, Side::west will print out "Side::west".
 *
 * @param os the ostream
 * @param s the side to print out.
 *
 * @return  the ostream
 */
inline std::ostream &operator<<(std::ostream &os, const Side<2> &s)
{
	switch (s.toInt()) {
		case Side<2>::west:
			os << "Side::west";
			break;
		case Side<2>::east:
			os << "Side::east";
			break;
		case Side<2>::south:
			os << "Side::south";
			break;
		case Side<2>::north:
			os << "Side::north";
			break;
	}
	return os;
}
/**
 * @brief ostream operator that prints a string representation of side enum.
 *
 * For example, Side::west will print out "Side::west".
 *
 * @param os the ostream
 * @param s the side to print out.
 *
 * @return  the ostream
 */
inline std::ostream &operator<<(std::ostream &os, const Side<3> &s)
{
	switch (s.toInt()) {
		case Side<3>::west:
			os << "Side::west";
			break;
		case Side<3>::east:
			os << "Side::east";
			break;
		case Side<3>::south:
			os << "Side::south";
			break;
		case Side<3>::north:
			os << "Side::north";
			break;
		case Side<3>::bottom:
			os << "Side::bottom";
			break;
		case Side<3>::top:
			os << "Side::top";
			break;
	}
	return os;
}
} // namespace Thunderegg
#endif
