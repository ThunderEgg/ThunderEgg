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

#ifndef THUNDEREGG_FACE_H
#define THUNDEREGG_FACE_H
#include <ThunderEgg/tpl/json.hpp>
#include <array>
namespace ThunderEgg
{
template <int D, int M> class Face
{
	private:
	/**
	 * @brief The value of the face
	 */
	unsigned char value;

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

	public:
	static constexpr int dimensionality = M;
	/**
	 * @brief the number of faces
	 */
	static constexpr size_t number_of = twopow(D - M) * choose(D, M);

	/**
	 * @brief Construct a new Face object
	 *
	 * @param value the value of the face to give
	 */
	Face(unsigned char value) : value(value)
	{
		static_assert(D <= 3 && D > 0, "Only up to 3 dimensions supported");
		static_assert(M >= 0 && M < D, "Invalid M value");
	}

	Face() : value(number_of)
	{
		static_assert(D <= 3 && D > 0, "Only up to 3 dimensions supported");
		static_assert(M >= 0 && M < D, "Invalid M value");
	}

	static Face<D, M> null()
	{
		return Face<D, M>(number_of);
	}

	/*
	 * these can be cleaned up when version is bumped to c++20
	 */

	template <int N = 0> static auto west() -> typename std::enable_if<D <= 3 && M == D - 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b000);
	}
	template <int N = 0> static auto east() -> typename std::enable_if<D <= 3 && M == D - 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b001);
	}
	template <int N = 0> static auto south() -> typename std::enable_if<D <= 3 && D >= 2 && M == D - 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b010);
	}
	template <int N = 0> static auto north() -> typename std::enable_if<D <= 3 && D >= 2 && M == D - 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b011);
	}
	template <int N = 0> static auto bottom() -> typename std::enable_if<D <= 3 && D >= 3 && M == D - 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b100);
	}
	template <int N = 0> static auto top() -> typename std::enable_if<D <= 3 && D >= 3 && M == D - 1 && N == N, Face<D, M>>::type
	{
		return Face<D, M>(0b101);
	}

	/**
	 * @brief Range class for Edge
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
		return {Face<D, D - 1>(value)};
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

template <int D> using Side = Face<D, D - 1>;
template <int D> Side<D> HigherSideOnAxis(int i)
{
	return Side<D>((i << 1) + 1);
}
template <int D> Side<D> LowerSideOnAxis(int i)
{
	return Side<D>(i << 1);
}
void to_json(nlohmann::json &j, const Side<1> &s);
void to_json(nlohmann::json &j, const Side<2> &s);
void to_json(nlohmann::json &j, const Side<3> &s);
void from_json(const nlohmann::json &j, Side<1> &s);
void from_json(const nlohmann::json &j, Side<2> &s);
void from_json(const nlohmann::json &j, Side<3> &s);
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
} // namespace ThunderEgg
#endif
