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

#ifndef THUNDEREGG_FACE_H
#define THUNDEREGG_FACE_H
/**
 * @file
 *
 * @brief Face class
 */

#include <ThunderEgg/tpl/json.hpp>
#include <array>
namespace ThunderEgg
{
/**
 * @brief Enum-style class for the faces of an n-dimensional cube
 *
 * @tparam D the dimension
 * @tparam M the dimensionality of the face
 *
 * For 3D, this class represents corners (`M = 0`), edges (`M = 1`), and faces (`M = 2`).
 */
template <int D, int M> class Face
{
	private:
	/**
	 * @brief The value of the face
	 */
	unsigned char value = number_of;

	/**
	 * @brief 2^N
	 */
	static constexpr size_t twopow(size_t N)
	{
		return 0b1 << N;
	}

	/**
	 * @brief N!
	 */
	static constexpr size_t factorial(size_t N)
	{
		return N == 0 ? 1 : N * factorial(N - 1);
	}

	/**
	 * @brief N choose K
	 */
	static constexpr size_t choose(size_t N, size_t K)
	{
		return factorial(N) / (factorial(K) * factorial(N - K));
	}

	static constexpr size_t sum(int N)
	{
		return N == -1 ? 0 : (twopow(D - N) * choose(D, N) + sum(N - 1));
	}

	public:
	static constexpr int dimensionality = M;
	/**
	 * @brief the number of faces
	 */
	static constexpr size_t number_of = twopow(D - M) * choose(D, M);

	static constexpr size_t sum_of_faces = M <= 0 ? 0 : sum(M - 1);
	/**
	 * @brief Construct a new Face object
	 *
	 * @param value the value of the face
	 */
	explicit Face(unsigned char value) : value(value)
	{
		static_assert(D <= 3 && D >= 0, "Only up to 3 dimensions supported");
		static_assert(M >= 0 && M < D, "Invalid M value");
	}

	/**
	 * @brief Construct a new Face object with a null value
	 */
	Face()
	{
		static_assert(D <= 3 && D >= 0, "Only up to 3 dimensions supported");
		static_assert(M >= 0 && M < D, "Invalid M value");
	}

	/**
	 * @brief null value
	 *
	 * @return Face<D, M>  the null value
	 */
	static Face<D, M> null()
	{
		return Face<D, M>(number_of);
	}

	/*
	 * these can be cleaned up when version is bumped to c++20
	 */

	/**
	 * @brief west side
	 */
	template <int N = 0> static auto west() -> typename std::enable_if<D <= 3 && M == D - 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b000);
	}
	/**
	 * @brief east side
	 */
	template <int N = 0> static auto east() -> typename std::enable_if<D <= 3 && M == D - 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b001);
	}
	/**
	 * @brief south side
	 */
	template <int N = 0> static auto south() -> typename std::enable_if<D <= 3 && D >= 2 && M == D - 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b010);
	}
	/**
	 * @brief north side
	 */
	template <int N = 0> static auto north() -> typename std::enable_if<D <= 3 && D >= 2 && M == D - 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b011);
	}
	/**
	 * @brief bottom side
	 */
	template <int N = 0> static auto bottom() -> typename std::enable_if<D <= 3 && D >= 3 && M == D - 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b100);
	}
	/**
	 * @brief top side
	 */
	template <int N = 0> static auto top() -> typename std::enable_if<D <= 3 && D >= 3 && M == D - 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b101);
	}

	/**
	 * @brief southwest corner
	 */
	template <int N = 0> static auto sw() -> typename std::enable_if<D == 2 && M == 0 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b00);
	}
	/**
	 * @brief southeast corner
	 */
	template <int N = 0> static auto se() -> typename std::enable_if<D == 2 && M == 0 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b01);
	}
	/**
	 * @brief northwest corner
	 */
	template <int N = 0> static auto nw() -> typename std::enable_if<D == 2 && M == 0 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b10);
	}
	/**
	 * @brief northeast corner
	 */
	template <int N = 0> static auto ne() -> typename std::enable_if<D == 2 && M == 0 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b11);
	}

	/**
	 * @brief bottom-south-west corner
	 */
	template <int N = 0> static auto bsw() -> typename std::enable_if<D == 3 && M == 0 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b000);
	}
	/**
	 * @brief bottom-south-east corner
	 */
	template <int N = 0> static auto bse() -> typename std::enable_if<D == 3 && M == 0 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b001);
	}
	/**
	 * @brief bottom-north-west corner
	 */
	template <int N = 0> static auto bnw() -> typename std::enable_if<D == 3 && M == 0 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b010);
	}
	/**
	 * @brief bottom-north-east corner
	 */
	template <int N = 0> static auto bne() -> typename std::enable_if<D == 3 && M == 0 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b011);
	}
	/**
	 * @brief top-south-west corner
	 */
	template <int N = 0> static auto tsw() -> typename std::enable_if<D == 3 && M == 0 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b100);
	}
	/**
	 * @brief top-south-east corner
	 */
	template <int N = 0> static auto tse() -> typename std::enable_if<D == 3 && M == 0 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b101);
	}
	/**
	 * @brief top-north-west corner
	 */
	template <int N = 0> static auto tnw() -> typename std::enable_if<D == 3 && M == 0 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b110);
	}
	/**
	 * @brief top-north-east corner
	 */
	template <int N = 0> static auto tne() -> typename std::enable_if<D == 3 && M == 0 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b111);
	}

	/**
	 * @brief bottom-south edge
	 */
	template <int N = 0> static auto bs() -> typename std::enable_if<D == 3 && M == 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b0000);
	}
	/**
	 * @brief bottom-north edge
	 */
	template <int N = 0> static auto bn() -> typename std::enable_if<D == 3 && M == 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b0001);
	}
	/**
	 * @brief top-south edge
	 */
	template <int N = 0> static auto ts() -> typename std::enable_if<D == 3 && M == 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b0010);
	}
	/**
	 * @brief top-north edge
	 */
	template <int N = 0> static auto tn() -> typename std::enable_if<D == 3 && M == 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b0011);
	}
	/**
	 * @brief bottom-west edge
	 */
	template <int N = 0> static auto bw() -> typename std::enable_if<D == 3 && M == 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b0100);
	}
	/**
	 * @brief bottom-east edge
	 */
	template <int N = 0> static auto be() -> typename std::enable_if<D == 3 && M == 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b0101);
	}
	/**
	 * @brief top-west edge
	 */
	template <int N = 0> static auto tw() -> typename std::enable_if<D == 3 && M == 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b0110);
	}
	/**
	 * @brief top-east edge
	 */
	template <int N = 0> static auto te() -> typename std::enable_if<D == 3 && M == 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b0111);
	}
	/**
	 * @brief south-west edge
	 */
	template <int N = 0> static auto sw() -> typename std::enable_if<D == 3 && M == 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b1000);
	}
	/**
	 * @brief south-east edge
	 */
	template <int N = 0> static auto se() -> typename std::enable_if<D == 3 && M == 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b1001);
	}
	/**
	 * @brief north-west edge
	 */
	template <int N = 0> static auto nw() -> typename std::enable_if<D == 3 && M == 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b1010);
	}
	/**
	 * @brief north-east edge
	 */
	template <int N = 0> static auto ne() -> typename std::enable_if<D == 3 && M == 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b1011);
	}

	/**
	 * @brief Range class for Face
	 *
	 * This provides the begin and end fuctions for iterator loops
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
			Face<D, M> s;

			public:
			/**
			 * @brief Construct a new Iterator object with the given Side value
			 *
			 * @param s_in the side
			 */
			explicit Iterator(Face<D, M> s) : s(s) {}
			/**
			 * @brief Increment the side value
			 *
			 * @return const Face<D,M>& the resulting value
			 */
			const Face<D, M> &operator++()
			{
				++s.value;
				return s;
			}
			/**
			 * @brief Get a reference to the side object
			 *
			 * @return const Face<D,M>&  the reference
			 */
			const Face<D, M> &operator*() const
			{
				return s;
			}
			/**
			 * @brief Get a pointer to the side object
			 *
			 * @return const Face<D,M>* the pointer
			 */
			const Face<D, M> *operator->() const
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
				return s.value == b.s.value;
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
				return s.value != b.s.value;
			}
		};
		/**
		 * @brief Returns an iterator with the lowest Side value
		 *
		 * @return Iterator the iterator
		 */
		Iterator begin()
		{
			return Iterator(Face<D, M>(0));
		}
		/**
		 * @brief Returns an iterator with Edge<D,M>::null()
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
	 * @brief Get the index for this Face
	 *
	 * @return size_t the index
	 */
	size_t getIndex() const
	{
		return value;
	}

	/**
	 * @brief Get the axis index of this side
	 */
	template <int N = 0> auto getAxisIndex() const -> typename std::enable_if<D <= 3 && D >= 1 && M == D - 1 && N == N, size_t>::type
	{
		return value >> 1;
	}

	/**
	 * @brief Return if this side is lower on it's axis
	 */
	template <int N = 0> auto isLowerOnAxis() const -> typename std::enable_if<D <= 3 && D >= 1 && M == D - 1 && N == N, bool>::type
	{
		return !(value & 0b1);
	}

	/**
	 * @brief Return if this side is higher on it's axis
	 */
	template <int N = 0> auto isHigherOnAxis() const -> typename std::enable_if<D <= 3 && D >= 1 && M == D - 1 && N == N, bool>::type
	{
		return value & 0b1;
	}

	/**
	 * @brief Get the face on the opposite side of the hypercube
	 *
	 * @return Face<D, M> the face on the opposite side
	 */
	Face<D, M> opposite() const
	{
		return Face<D, M>(value ^ ~((~0u) << (D - M)));
	}

	/**
	 * @brief Get the Sides that this Face lies on
	 *
	 * @return std::array<Face<D, D - 1>, D - M> the sides
	 */
	std::array<Face<D, D - 1>, D - M> getSides() const
	{
		std::array<Face<D, D - 1>, D - M> sides;
		if constexpr (D > 1 && M == 0) {
			for (size_t i = 0; i < D; i++) {
				size_t side = 2 * i;
				if ((1 << i) & value) {
					side |= 1;
				}
				sides[i] = Face<D, D - 1>(side);
			}
		} else if constexpr (D == 3 && M == 1) {
			size_t not_axis   = value >> 2;
			size_t curr_index = 0;
			for (unsigned char i = 0; i < D; i++) {
				if (i != not_axis) {
					unsigned char bit = (value & (0b1 << curr_index)) >> curr_index;
					sides[curr_index] = Face<D, D - 1>((i << 1) ^ bit);
					curr_index++;
				}
			}
		} else {
			sides[0] = Face<D, D - 1>(value);
		}
		return sides;
	}

	/**
	 * @brief Compare the enum indexes.
	 *
	 * @param other The other side.
	 *
	 * @return Whether or not the index of this face is lower than the other face.
	 */
	bool operator<(const Face<D, M> &other) const
	{
		return value < other.value;
	}

	/**
	 * @brief Equals operator.
	 *
	 * @param other The other side.
	 *
	 * @return Whether or not the value of this face equals the value other face.
	 */
	bool operator==(const Face<D, M> &other) const
	{
		return value == other.value;
	}

	/**
	 * @brief Not Equals operator.
	 *
	 * @param other The other side.
	 *
	 * @return Whether or not the value of this face equals the value other face.
	 */
	bool operator!=(const Face<D, M> &other) const
	{
		return value != other.value;
	}
};

/**
 * @brief Side class
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> using Side = Face<D, D - 1>;
/**
 * @brief Get the higher side on a given axis
 *
 * @tparam D the number of Cartesian dimensions
 * @param axis the axis
 * @return Side<D> the side
 */
template <int D> Side<D> HigherSideOnAxis(size_t axis)
{
	return Side<D>((axis << 1) + 1);
}
/**
 * @brief Get the lower side on a given axis
 *
 * @tparam D the number of Cartesian dimensions
 * @param axis the axis
 * @return Side<D> the side
 */
template <int D> Side<D> LowerSideOnAxis(size_t axis)
{
	return Side<D>(axis << 1);
}
void to_json(tpl::nlohmann::json &j, const Side<1> &s);
void to_json(tpl::nlohmann::json &j, const Side<2> &s);
void to_json(tpl::nlohmann::json &j, const Side<3> &s);
void from_json(const tpl::nlohmann::json &j, Side<1> &s);
void from_json(const tpl::nlohmann::json &j, Side<2> &s);
void from_json(const tpl::nlohmann::json &j, Side<3> &s);
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
std::ostream &operator<<(std::ostream &os, const Side<1> &s);
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
std::ostream &operator<<(std::ostream &os, const Side<2> &s);
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
std::ostream &operator<<(std::ostream &os, const Side<3> &s);

/**
 * @brief Edge class
 */
using Edge = Face<3, 1>;
/**
 * @brief ostream operator that prints a string representation of edge enum.
 *
 * For example, Edge::sw() will print out "Edge::sw()".
 *
 * @param os the ostream
 * @param o the edge
 *
 * @return  the ostream
 */
std::ostream &operator<<(std::ostream &os, const Edge &o);
void          to_json(tpl::nlohmann::json &j, const Edge &o);
void          from_json(const tpl::nlohmann::json &j, Edge &o);

/**
 * @brief Corner class
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> using Corner = Face<D, 0>;
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
std::ostream &operator<<(std::ostream &os, const Corner<2> &o);
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
std::ostream &operator<<(std::ostream &os, const Corner<3> &o);
void          to_json(tpl::nlohmann::json &j, const Corner<2> &o);
void          to_json(tpl::nlohmann::json &j, const Corner<3> &o);
void          from_json(const tpl::nlohmann::json &j, Corner<2> &o);
void          from_json(const tpl::nlohmann::json &j, Corner<3> &o);
} // namespace ThunderEgg
#endif
