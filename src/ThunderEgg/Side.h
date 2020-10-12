/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
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
#include <ThunderEgg/tpl/json.hpp>
#include <array>
#include <iostream>
#include <numeric>
#include <vector>
namespace ThunderEgg
{
/**
 * @brief An enum-style class that represents the sides of a patch
 *
 * The sides of the patch are named in the following way:
 *
 * Orthogonal to axis | Lower on axis | Higher on axis
 * ------------------ | ------------- | --------------
 *  x-axis            | west          | east
 *  y-axis            | south         | north
 *  z-axis            | bottom        | top
 *
 * @tparam D the number of cartesian dimensions
 */
template <int D> class Side
{
	private:
	/**
	 * @brief the value of the enum
	 */
	unsigned char val = num_sides;

	public:
	/**
	 * @brief The number of possible side values
	 *
	 */
	static constexpr size_t num_sides = 2 * D;
	/**
	 * @brief Default constructor that sets the value to null()
	 */
	Side() : val(num_sides) {}

	/**
	 * @brief Construct a new Side object
	 *
	 * @param val_in the value of the side object
	 */
	explicit Side(unsigned char val_in) : val(val_in) {}
	/**
	 * @brief Get the side that is lower on the specified axis
	 *
	 * @param axis the axis
	 * @return Side<D> the resulting side
	 */
	static Side<D> LowerSideOnAxis(size_t axis)
	{
		return Side<D>(2 * axis);
	}
	/**
	 * @brief Get the side that is higher on the specified axis
	 *
	 * @param axis the axis
	 * @return Side<D> the resulting side
	 */
	static Side<D> HigherSideOnAxis(size_t axis)
	{
		return Side<D>(2 * axis + 1);
	}
	/**
	 * @brief Construct a new side object for the western side
	 */
	static Side<D> west()
	{
		return Side(0);
	}
	/**
	 * @brief Construct a new side object for the eastern side
	 */
	static Side<D> east()
	{
		return Side(1);
	}
	/**
	 * @brief Construct a new side object for the southern side
	 */
	static Side<D> south()
	{
		return Side(2);
	}
	/**
	 * @brief Construct a new side object for the northern side
	 */
	static Side<D> north()
	{
		return Side(3);
	}
	/**
	 * @brief Construct a new side object for the bottom side
	 */
	static Side<D> bottom()
	{
		return Side(4);
	}
	/**
	 * @brief Construct a new side object for the top side
	 */
	static Side<D> top()
	{
		return Side(5);
	}
	/**
	 * @brief Construct a new side object with a null value.
	 *
	 * The value of the side will be num_sides
	 */
	static Side<D> null()
	{
		return Side(num_sides);
	}
	/**
	 * @brief Range class for Side
	 *
	 * It provides the begin and end fuctions for iterator loops
	 */
	class Range
	{
		public:
		/**
		 * @brief Input iterator for Side values
		 */
		class Iterator : public std::input_iterator_tag
		{
			private:
			/**
			 * @brief The current side
			 */
			Side<D> s;

			public:
			/**
			 * @brief Construct a new Iterator object with the given Side value
			 *
			 * @param s_in the side
			 */
			explicit Iterator(Side<D> s_in) : s(s_in) {}
			/**
			 * @brief Increment the side value
			 *
			 * @return const Side<D>& the resulting value
			 */
			const Side<D> &operator++()
			{
				++s.val;
				return s;
			}
			/**
			 * @brief Get a reference to the side object
			 *
			 * @return const Side<D>&  the reference
			 */
			const Side<D> &operator*() const
			{
				return s;
			}
			/**
			 * @brief Get a pointer to the side object
			 *
			 * @return const Side<D>* the pointer
			 */
			const Side<D> *operator->() const
			{
				return &s;
			}
			/**
			 * @brief Check the iterators reference the same value
			 *
			 * @param b the other iterator
			 * @return true if the same
			 * @return false if different
			 */
			bool operator==(const Iterator &b) const
			{
				return s.val == b.s.val;
			}
			/**
			 * @brief Check the iterators don't reference the same value
			 *
			 * @param b the other iterator
			 * @return true if different
			 * @return false if the same
			 */
			bool operator!=(const Iterator &b) const
			{
				return s.val != b.s.val;
			}
		};
		/**
		 * @brief Returns an iterator with the lowest Side value
		 *
		 * @return Iterator the iterator
		 */
		Iterator begin()
		{
			return Iterator(Side<D>(0));
		}
		/**
		 * @brief Returns an iterator with Side<D>::null()
		 *
		 * @return Iterator the iterator
		 */
		Iterator end()
		{
			return Iterator(null());
		}
	};
	/**
	 * @brief Get a range of values that can be iterated over
	 *
	 * @return Range the range of values
	 */
	static Range getValues()
	{
		return Range();
	}

	/**
	 * @brief Get the integer index of the side
	 *
	 * @return The index of the side.
	 */
	inline size_t getIndex() const
	{
		return val;
	}
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
	 * @brief Return the index of the axis that the side lies on.
	 *
	 * 0 will be x axis, 1 will be y axis, etc.
	 */
	inline size_t getAxisIndex() const
	{
		return val >> 1;
	}
	/**
	 * @brief Return the opposite side of the cube.
	 *
	 * For example: the opposite of east is west.
	 *
	 * @return The opposite side.
	 */
	Side opposite() const
	{
		Side<D> retval = *this;
		// flip least significant bit
		retval.val ^= 0x1;
		return retval;
	}
	/**
	 * @brief Compare the enum indexes.
	 *
	 * @param other The other side.
	 *
	 * @return Whether or not the index of this side is lower than the other side.
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
};
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
inline std::ostream &operator<<(std::ostream &os, const Side<1> &s)
{
	if (s == Side<1>::east()) {
		os << "Side<1>::east()";
	} else if (s == Side<1>::west()) {
		os << "Side<1>::west()";
	} else if (s == Side<1>::null()) {
		os << "Side<1>::null()";
	} else {
		os << "Side<1> undefined value: " << s.getIndex();
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
inline std::ostream &operator<<(std::ostream &os, const Side<2> &s)
{
	if (s == Side<2>::east()) {
		os << "Side<2>::east()";
	} else if (s == Side<2>::west()) {
		os << "Side<2>::west()";
	} else if (s == Side<2>::south()) {
		os << "Side<2>::south()";
	} else if (s == Side<2>::north()) {
		os << "Side<2>::north()";
	} else if (s == Side<2>::null()) {
		os << "Side<2>::null()";
	} else {
		os << "Side<2> undefined value: " << s.getIndex();
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
	if (s == Side<3>::east()) {
		os << "Side<3>::east()";
	} else if (s == Side<3>::west()) {
		os << "Side<3>::west()";
	} else if (s == Side<3>::south()) {
		os << "Side<3>::south()";
	} else if (s == Side<3>::north()) {
		os << "Side<3>::north()";
	} else if (s == Side<3>::bottom()) {
		os << "Side<3>::bottom()";
	} else if (s == Side<3>::top()) {
		os << "Side<3>::top()";
	} else if (s == Side<3>::null()) {
		os << "Side<3>::null()";
	} else {
		os << "Side<3> undefined value: " << s.getIndex();
	}
	return os;
}
void to_json(nlohmann::json &j, const Side<1> &s);
void to_json(nlohmann::json &j, const Side<2> &s);
void to_json(nlohmann::json &j, const Side<3> &s);
void from_json(const nlohmann::json &j, Side<1> &s);
void from_json(const nlohmann::json &j, Side<2> &s);
void from_json(const nlohmann::json &j, Side<3> &s);
} // namespace ThunderEgg
#endif
