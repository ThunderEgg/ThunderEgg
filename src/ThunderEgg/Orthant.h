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

#ifndef THUNDEREGG_ORTHANT_H
#define THUNDEREGG_ORTHANT_H
#include <ThunderEgg/Side.h>
#include <array>
#include <iostream>
#include <numeric>
#include <vector>
namespace ThunderEgg
{
/**
 * @brief An enum-style class that represents the octants of a cube.
 *
 * The octants are named in the following way:
 *
 * bse means Bottom-South-West, which is in the corner of where the bottom, south, and west sides
 * meet.
 */
template <size_t D> class Orthant
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
	explicit Orthant(const unsigned char val_in) : val(val_in) {}
	/**
	 * @brief Default constructor that initializes the value to null().
	 */
	Orthant() {}
	/**
	 * @brief null value
	 */
	static Orthant<D> null()
	{
		return Orthant<D>(num_orthants);
	}
	/**
	 * @brief Lower half of line.
	 */
	static Orthant<1> lower()
	{
		return Orthant<1>(0b0);
	}
	/**
	 * @brief Upper half of line.
	 */
	static Orthant<1> upper()
	{
		return Orthant<1>(0b1);
	}
	/**
	 * @brief South-West quadrant of square.
	 */
	static Orthant<2> sw()
	{
		return Orthant<2>(0b00);
	}
	/**
	 * @brief South-East quadrant of square.
	 */
	static Orthant<2> se()
	{
		return Orthant<2>(0b01);
	}
	/**
	 * @brief North-West quadrant of square.
	 */
	static Orthant<2> nw()
	{
		return Orthant<2>(0b10);
	}
	/**
	 * @brief North-East quadrant of square.
	 */
	static Orthant<2> ne()
	{
		return Orthant<2>(0b11);
	}
	/**
	 * @brief Bottom-South-West octant of cube.
	 */
	static Orthant<3> bsw()
	{
		return Orthant<3>(0b000);
	}
	/**
	 * @brief Bottom-South-East octant of cube.
	 */
	static Orthant<3> bse()
	{
		return Orthant<3>(0b001);
	}
	/**
	 * @brief Bottom-North-West octant of cube.
	 */
	static Orthant<3> bnw()
	{
		return Orthant<3>(0b010);
	}
	/**
	 * @brief Bottom-North-East octant of cube.
	 */
	static Orthant<3> bne()
	{
		return Orthant<3>(0b011);
	}
	/**
	 * @brief Top-South-West octant of cube.
	 */
	static Orthant<3> tsw()
	{
		return Orthant<3>(0b100);
	}
	/**
	 * @brief Top-South-East octant of cube.
	 */
	static Orthant<3> tse()
	{
		return Orthant<3>(0b101);
	}
	/**
	 * @brief Top-North-West octant of cube.
	 */
	static Orthant<3> tnw()
	{
		return Orthant<3>(0b110);
	}
	/**
	 * @brief Top-North-East octant of cube.
	 */
	static Orthant<3> tne()
	{
		return Orthant<3>(0b111);
	}
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
			explicit Iterator(Orthant<D> o_in) : o(o_in) {}
			/**
			 * @brief Increment the Orthant value
			 *
			 * @return const Orthant<D>& the resulting value
			 */
			const Orthant<D> &operator++()
			{
				++o.val;
				return o;
			}
			/**
			 * @brief Get a reference to the Orthant object
			 *
			 * @return const Orthant<D>&  the reference
			 */
			const Orthant<D> &operator*() const
			{
				return o;
			}
			/**
			 * @brief Get a pointer to the Orthant object
			 *
			 * @return const Orthant<D>* the pointer
			 */
			const Orthant<D> *operator->() const
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
		 * @brief Returns an iterator with the lowest Orthant value
		 *
		 * @return Iterator the iterator
		 */
		Iterator begin()
		{
			return Iterator(Orthant<D>(0));
		}
		/**
		 * @brief Returns an iterator with Orthant<D>::null()
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
	 * @brief Get the integer value of the octant.
	 *
	 * @return The integer value.
	 */
	size_t getIndex() const
	{
		return val;
	}
	/**
	 * @brief Return the octant that neighbors this octant on a particular side.
	 *
	 * @param s the side of the octant that you want the neighbor of.
	 *
	 * @return  The octant that neighbors on that side.
	 */
	Orthant<D> getNbrOnSide(Side<D> s) const
	{
		Orthant<D> retval = *this;
		// flip the bit for that side
		retval.val ^= (0x1 << s.getAxisIndex());
		return retval;
	}
	/**
	 * @brief Get the sides of the octant that are on the interior of the cube.
	 *
	 * @return The sides of the octant that are on the interior of the cube.
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
	 * @brief Get the sides of the octant that are on the exterior of the cube.
	 *
	 * @return The sides of the octant that are on the exterior of the cube.
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
	 * @brief Return whether or not the octant lies on a particular side of a cube.
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
	 * @brief From the point of view of an axis, get orthant that this orthant lies on in the D-1
	 * dimension
	 *
	 * @param axis the axis
	 * @return Orthant<D - 1> the resulting orthant
	 */
	Orthant<D - 1> collapseOnAxis(size_t axis) const
	{
		size_t upper_mask = (~0x0U) << axis;
		return Orthant<D - 1>(((val >> 1) & upper_mask) | (val & ~upper_mask));
	}
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
	static std::array<Orthant, num_orthants / 2> getValuesOnSide(Side<D> s)
	{
		unsigned int bit_to_insert = s.getAxisIndex();
		unsigned int set_bit       = s.isLowerOnAxis() ? 0 : 1;
		unsigned int lower_mask    = ~((~0x0U) << bit_to_insert);
		unsigned int upper_mask    = (~0x0U) << (bit_to_insert + 1);

		std::array<Orthant<D>, Orthant<D>::num_orthants / 2> retval;
		for (size_t i = 0; i < Orthant<D>::num_orthants / 2; i++) {
			size_t value = (i << 1) & upper_mask;
			value |= i & lower_mask;
			value |= set_bit << bit_to_insert;
			retval[i] = Orthant<D>(value);
		}
		return retval;
	}
	/**
	 * @brief Equals operator.
	 *
	 * @param other The other octant.
	 *
	 * @return Whether or not the value of this octant equals the value other octant.
	 */
	bool operator==(const Orthant<D> &other) const
	{
		return val == other.val;
	}
	/**
	 * @brief Not Equals operator.
	 *
	 * @param other The other octant.
	 *
	 * @return Whether or not the value of this octant is not equal the value other octant.
	 */
	bool operator!=(const Orthant<D> &other) const
	{
		return val != other.val;
	}
	/**
	 * @brief Less Tan operator.
	 *
	 * @param other The other octant.
	 *
	 * @return Whether or not the value of this octant is less than the value other octant.
	 */
	bool operator<(const Orthant<D> &other) const
	{
		return val < other.val;
	}
};

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
inline std::ostream &operator<<(std::ostream &os, const Orthant<1> &o)
{
	if (o == Orthant<1>::lower()) {
		os << "Orthant<1>::lower()";
	} else if (o == Orthant<1>::upper()) {
		os << "Orthant<1>::upper()";
	} else if (o == Orthant<1>::null()) {
		os << "Orthant<1>::null()";
	} else {
		os << "Orthant<1> invalid value: " << o.getIndex();
	}
	return os;
}
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
inline std::ostream &operator<<(std::ostream &os, const Orthant<2> &o)
{
	if (o == Orthant<2>::sw()) {
		os << "Orthant<2>::sw()";
	} else if (o == Orthant<2>::se()) {
		os << "Orthant<2>::se()";
	} else if (o == Orthant<2>::nw()) {
		os << "Orthant<2>::nw()";
	} else if (o == Orthant<2>::ne()) {
		os << "Orthant<2>::ne()";
	} else if (o == Orthant<2>::null()) {
		os << "Orthant<2>::null()";
	} else {
		os << "Orthant<2> invalid value: " << o.getIndex();
	}
	return os;
}
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
inline std::ostream &operator<<(std::ostream &os, const Orthant<3> &o)
{
	if (o == Orthant<3>::bsw()) {
		os << "Orthant<3>::bsw()";
	} else if (o == Orthant<3>::bse()) {
		os << "Orthant<3>::bse()";
	} else if (o == Orthant<3>::bnw()) {
		os << "Orthant<3>::bnw()";
	} else if (o == Orthant<3>::bne()) {
		os << "Orthant<3>::bne()";
	} else if (o == Orthant<3>::tsw()) {
		os << "Orthant<3>::tsw()";
	} else if (o == Orthant<3>::tse()) {
		os << "Orthant<3>::tse()";
	} else if (o == Orthant<3>::tnw()) {
		os << "Orthant<3>::tnw()";
	} else if (o == Orthant<3>::tne()) {
		os << "Orthant<3>::tne()";
	} else if (o == Orthant<3>::null()) {
		os << "Orthant<3>::null()";
	} else {
		os << "Orthant<3> invalid value: " << o.getIndex();
	}
	return os;
}

} // namespace ThunderEgg
#endif
