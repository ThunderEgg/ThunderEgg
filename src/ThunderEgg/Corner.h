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

#ifndef THUNDEREGG_CORNER_H
#define THUNDEREGG_CORNER_H
#include <ThunderEgg/Face.h>
#include <array>
#include <iostream>
#include <numeric>
#include <vector>
namespace ThunderEgg
{
/**
 * @brief An enum-style class that represents the corners of a patch
 *
 * The corners are named in the following way:
 *
 * bse means Bottom-South-West, which is in the corner of where the bottom, south, and west sides
 * meet.
 */
template <int D> class Corner
{
	private:
	/**
	 * @brief the value of the enum.
	 */
	unsigned char val = num_corners;

	public:
	static constexpr size_t num_corners    = D > 1 ? 1 << D : 0;
	static constexpr size_t dimensionality = 0;
	/**
	 * @brief Create new Corner<D> with given value.
	 *
	 * @param val_in the value
	 */
	explicit Corner(const unsigned char val_in) : val(val_in) {}
	/**
	 * @brief Default constructor that initializes the value to null().
	 */
	Corner() {}
	/**
	 * @brief null value
	 */
	static Corner<D> null()
	{
		return Corner<D>(num_corners);
	}
	/**
	 * @brief South-West quadrant of square.
	 */
	static Corner<2> sw()
	{
		return Corner<2>(0b00);
	}
	/**
	 * @brief South-East quadrant of square.
	 */
	static Corner<2> se()
	{
		return Corner<2>(0b01);
	}
	/**
	 * @brief North-West quadrant of square.
	 */
	static Corner<2> nw()
	{
		return Corner<2>(0b10);
	}
	/**
	 * @brief North-East quadrant of square.
	 */
	static Corner<2> ne()
	{
		return Corner<2>(0b11);
	}
	/**
	 * @brief Bottom-South-West corner of cube.
	 */
	static Corner<3> bsw()
	{
		return Corner<3>(0b000);
	}
	/**
	 * @brief Bottom-South-East corner of cube.
	 */
	static Corner<3> bse()
	{
		return Corner<3>(0b001);
	}
	/**
	 * @brief Bottom-North-West corner of cube.
	 */
	static Corner<3> bnw()
	{
		return Corner<3>(0b010);
	}
	/**
	 * @brief Bottom-North-East corner of cube.
	 */
	static Corner<3> bne()
	{
		return Corner<3>(0b011);
	}
	/**
	 * @brief Top-South-West corner of cube.
	 */
	static Corner<3> tsw()
	{
		return Corner<3>(0b100);
	}
	/**
	 * @brief Top-South-East corner of cube.
	 */
	static Corner<3> tse()
	{
		return Corner<3>(0b101);
	}
	/**
	 * @brief Top-North-West corner of cube.
	 */
	static Corner<3> tnw()
	{
		return Corner<3>(0b110);
	}
	/**
	 * @brief Top-North-East corner of cube.
	 */
	static Corner<3> tne()
	{
		return Corner<3>(0b111);
	}
	/**
	 * @brief Range class for Corner
	 *
	 * It provides the begin and end fuctions for iterator loops
	 */
	class Range
	{
		public:
		/**
		 * @brief Input iterator for Corner values
		 */
		class Iterator : public std::input_iterator_tag
		{
			private:
			/**
			 * @brief The current side
			 */
			Corner<D> o;

			public:
			/**
			 * @brief Construct a new Iterator object with the given Corner value
			 *
			 * @param o_in the corner
			 */
			explicit Iterator(Corner<D> o_in) : o(o_in) {}
			/**
			 * @brief Increment the Corner value
			 *
			 * @return const Corner<D>& the resulting value
			 */
			const Corner<D> &operator++()
			{
				++o.val;
				return o;
			}
			/**
			 * @brief Get a reference to the Corner object
			 *
			 * @return const Corner<D>&  the reference
			 */
			const Corner<D> &operator*() const
			{
				return o;
			}
			/**
			 * @brief Get a pointer to the Corner object
			 *
			 * @return const Corner<D>* the pointer
			 */
			const Corner<D> *operator->() const
			{
				return &o;
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
				return o.val == b.o.val;
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
				return o.val != b.o.val;
			}
		};
		/**
		 * @brief Returns an iterator with the lowest Corner value
		 *
		 * @return Iterator the iterator
		 */
		Iterator begin()
		{
			return Iterator(Corner<D>(0));
		}
		/**
		 * @brief Returns an iterator with Corner<D>::null()
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
	 * @brief Get the integer value of the corner.
	 *
	 * @return The integer value.
	 */
	size_t getIndex() const
	{
		return val;
	}
	/**
	 * @brief Return the corner that neighbors this corner on a particular side.
	 *
	 * @param s the side of the corner that you want the neighbor of.
	 *
	 * @return  The corner that neighbors on that side.
	 */
	Corner<D> getNbrOnSide(Side<D> s) const
	{
		Corner<D> retval = *this;
		// flip the bit for that side
		retval.val ^= (0x1 << s.getAxisIndex());
		return retval;
	}
	/**
	 * @brief Get the sides of the corner that are on the interior of the cube.
	 *
	 * @return The sides of the corner that are on the interior of the cube.
	 */
	std::array<Side<D>, D> getInteriorSides() const
	{
		std::array<Side<D>, D> retval;
		for (size_t i = 0; i < D; i++) {
			size_t side = 2 * i;
			if (!((1 << i) & val)) {
				side |= 1;
			}
			retval[i] = Side<D>(side);
		}
		return retval;
	}
	/**
	 * @brief Get the sides of the corner that are on the exterior of the cube.
	 *
	 * @return The sides of the corner that are on the exterior of the cube.
	 */
	std::array<Side<D>, D> getExteriorSides() const
	{
		std::array<Side<D>, D> retval;
		for (size_t i = 0; i < D; i++) {
			size_t side = 2 * i;
			if ((1 << i) & val) {
				side |= 1;
			}
			retval[i] = Side<D>(side);
		}
		return retval;
	}
	/**
	 * @brief Return whether or not the corner lies on a particular side of a cube.
	 *
	 * @param s the side of the cube.j
	 *
	 * @return Whether or not it lies on that side.
	 */
	bool isOnSide(Side<D> s) const
	{
		int  idx        = s.getIndex() / 2;
		int  remainder  = s.getIndex() % 2;
		bool is_bit_set = val & (0x1 << idx);
		return is_bit_set == remainder;
	}
	/**
	 * @brief Return if the corner is lower on a given axis
	 *
	 * @param axis the axis
	 * @return if corner is lower on a given axis
	 */
	bool isLowerOnAxis(size_t axis) const
	{
		return !(val & (0b1 << axis));
	}
	/**
	 * @brief Return if the corner is higher on a given axis
	 *
	 * @param axis the axis
	 * @return if corner is higher on a given axis
	 */
	bool isHigherOnAxis(size_t axis) const
	{
		return val & (0b1 << axis);
	}
	/**
	 * @brief From the point of view of an axis, get corner that this corner lies on in the D-1
	 * dimension
	 *
	 * @param axis the axis
	 * @return Corner<D - 1> the resulting corner
	 */
	Corner<D - 1> collapseOnAxis(size_t axis) const
	{
		size_t upper_mask = (~0x0U) << axis;
		return Corner<D - 1>(((val >> 1) & upper_mask) | (val & ~upper_mask));
	}
	/**
	 * @brief get the corner on the opposite side
	 *
	 * @return Corner<D> the corner on the opposite side
	 */
	Corner<D> opposite() const
	{
		unsigned char mask = ~((~0x0U) << D);
		return Corner<D>(val ^ mask);
	}
	/**
	 * @brief Get an array of all Corner<D> values that lie on a particular side of the cube.
	 *
	 * When the two axis that the side lies on are arranged in the following way, the corners are
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
	static std::array<Corner, num_corners / 2> getValuesOnSide(Side<D> s)
	{
		unsigned int bit_to_insert = s.getAxisIndex();
		unsigned int set_bit       = s.isLowerOnAxis() ? 0 : 1;
		unsigned int lower_mask    = ~((~0x0U) << bit_to_insert);
		unsigned int upper_mask    = (~0x0U) << (bit_to_insert + 1);

		std::array<Corner<D>, Corner<D>::num_corners / 2> retval;
		for (size_t i = 0; i < Corner<D>::num_corners / 2; i++) {
			size_t value = (i << 1) & upper_mask;
			value |= i & lower_mask;
			value |= set_bit << bit_to_insert;
			retval[i] = Corner<D>(value);
		}
		return retval;
	}
	/**
	 * @brief Equals operator.
	 *
	 * @param other The other corner.
	 *
	 * @return Whether or not the value of this corner equals the value other corner.
	 */
	bool operator==(const Corner<D> &other) const
	{
		return val == other.val;
	}
	/**
	 * @brief Not Equals operator.
	 *
	 * @param other The other corner.
	 *
	 * @return Whether or not the value of this corner is not equal the value other corner.
	 */
	bool operator!=(const Corner<D> &other) const
	{
		return val != other.val;
	}
	/**
	 * @brief Less Tan operator.
	 *
	 * @param other The other corner.
	 *
	 * @return Whether or not the value of this corner is less than the value other corner.
	 */
	bool operator<(const Corner<D> &other) const
	{
		return val < other.val;
	}
};

/**
 * @brief ostream operator that prints a string representation of Corner<1> enum.
 *
 * For example, Corner<1>::lower() will print out "Corner<1>::lower()".
 *
 * @param os the ostream
 * @param o the corner
 *
 * @return  the ostream
 */
inline std::ostream &operator<<(std::ostream &os, const Corner<1> &o)
{
	if (o == Corner<1>::null()) {
		os << "Corner<1>::null()";
	} else {
		os << "Corner<1> invalid value: " << o.getIndex();
	}
	return os;
}
/**
 * @brief ostream operator that prints a string representation of quadrant enum.
 *
 * For example, Corner<2>::sw() will print out "Corner<2>::sw()".
 *
 * @param os the ostream
 * @param o the corner
 *
 * @return  the ostream
 */
inline std::ostream &operator<<(std::ostream &os, const Corner<2> &o)
{
	if (o == Corner<2>::sw()) {
		os << "Corner<2>::sw()";
	} else if (o == Corner<2>::se()) {
		os << "Corner<2>::se()";
	} else if (o == Corner<2>::nw()) {
		os << "Corner<2>::nw()";
	} else if (o == Corner<2>::ne()) {
		os << "Corner<2>::ne()";
	} else if (o == Corner<2>::null()) {
		os << "Corner<2>::null()";
	} else {
		os << "Corner<2> invalid value: " << o.getIndex();
	}
	return os;
}
/**
 * @brief ostream operator that prints a string representation of corner enum.
 *
 * For example, Corner<3>::bsw() will print out "Corner<3>::bsw()".
 *
 * @param os the ostream
 * @param o the corner
 *
 * @return  the ostream
 */
inline std::ostream &operator<<(std::ostream &os, const Corner<3> &o)
{
	if (o == Corner<3>::bsw()) {
		os << "Corner<3>::bsw()";
	} else if (o == Corner<3>::bse()) {
		os << "Corner<3>::bse()";
	} else if (o == Corner<3>::bnw()) {
		os << "Corner<3>::bnw()";
	} else if (o == Corner<3>::bne()) {
		os << "Corner<3>::bne()";
	} else if (o == Corner<3>::tsw()) {
		os << "Corner<3>::tsw()";
	} else if (o == Corner<3>::tse()) {
		os << "Corner<3>::tse()";
	} else if (o == Corner<3>::tnw()) {
		os << "Corner<3>::tnw()";
	} else if (o == Corner<3>::tne()) {
		os << "Corner<3>::tne()";
	} else if (o == Corner<3>::null()) {
		os << "Corner<3>::null()";
	} else {
		os << "Corner<3> invalid value: " << o.getIndex();
	}
	return os;
}
void to_json(nlohmann::json &j, const Corner<1> &o);
void to_json(nlohmann::json &j, const Corner<2> &o);
void to_json(nlohmann::json &j, const Corner<3> &o);
void from_json(const nlohmann::json &j, Corner<1> &o);
void from_json(const nlohmann::json &j, Corner<2> &o);
void from_json(const nlohmann::json &j, Corner<3> &o);
} // namespace ThunderEgg
#endif
