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

#ifndef THUNDEREGG_EDGE_H
#define THUNDEREGG_EDGE_H
#include <ThunderEgg/Side.h>
#include <array>
#include <iostream>
#include <numeric>
#include <vector>
namespace ThunderEgg
{
/**
 * @brief An enum-style class that represents the edges of a 3D patch
 *
 * The edges are named in the following way:
 *
 * bs means Bottom-South-West, which is in the edge of where the bottom and south sides
 * meet.
 */
template <int D> class Edge
{
	private:
	/**
	 * @brief the value of the enum.
	 */
	unsigned char val = num_edges;

	public:
	/**
	 * @brief The number of edges
	 */
	static constexpr size_t num_edges = D > 2 ? 12 : 0;

	/**
	 * @brief Dimensionality
	 */
	static constexpr size_t dimensionality = D > 2 ? 1 : 0;

	/**
	 * @brief Create new Edge<D> with given value.
	 *
	 * @param val_in the value
	 */
	explicit Edge(const unsigned char val_in) : val(val_in) {}
	/**
	 * @brief Default constructor that initializes the value to null().
	 */
	Edge() {}
	/**
	 * @brief null value
	 */
	static Edge<D> null()
	{
		return Edge<D>(num_edges);
	}
	/**
	 * @brief Bottom-South edge of cube.
	 */
	static Edge<3> bs()
	{
		return Edge<3>(0b0000);
	}
	/**
	 * @brief Bottom-South edge of cube.
	 */
	static Edge<3> bn()
	{
		return Edge<3>(0b0001);
	}
	/**
	 * @brief Top-North edge of cube.
	 */
	static Edge<3> ts()
	{
		return Edge<3>(0b0010);
	}
	/**
	 * @brief Top-North edge of cube.
	 */
	static Edge<3> tn()
	{
		return Edge<3>(0b0011);
	}
	/**
	 * @brief Bottom-West edge of cube.
	 */
	static Edge<3> bw()
	{
		return Edge<3>(0b0100);
	}
	/**
	 * @brief Bottom-East edge of cube.
	 */
	static Edge<3> be()
	{
		return Edge<3>(0b0101);
	}
	/**
	 * @brief Top-West edge of cube.
	 */
	static Edge<3> tw()
	{
		return Edge<3>(0b0110);
	}
	/**
	 * @brief Top-East edge of cube.
	 */
	static Edge<3> te()
	{
		return Edge<3>(0b0111);
	}
	/**
	 * @brief South-West edge of cube.
	 */
	static Edge<3> sw()
	{
		return Edge<3>(0b1000);
	}
	/**
	 * @brief South-West edge of cube.
	 */
	static Edge<3> se()
	{
		return Edge<3>(0b1001);
	}
	/**
	 * @brief North-East edge of cube.
	 */
	static Edge<3> nw()
	{
		return Edge<3>(0b1010);
	}
	/**
	 * @brief North-East edge of cube.
	 */
	static Edge<3> ne()
	{
		return Edge<3>(0b1011);
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
		 * @brief Input iterator for Edge values
		 */
		class Iterator : public std::input_iterator_tag
		{
			private:
			/**
			 * @brief The current side
			 */
			Edge<D> o;

			public:
			/**
			 * @brief Construct a new Iterator object with the given Edge value
			 *
			 * @param o_in the edge
			 */
			explicit Iterator(Edge<D> o_in) : o(o_in) {}
			/**
			 * @brief Increment the Edge value
			 *
			 * @return const Edge<D>& the resulting value
			 */
			const Edge<D> &operator++()
			{
				++o.val;
				return o;
			}
			/**
			 * @brief Get a reference to the Edge object
			 *
			 * @return const Edge<D>&  the reference
			 */
			const Edge<D> &operator*() const
			{
				return o;
			}
			/**
			 * @brief Get a pointer to the Edge object
			 *
			 * @return const Edge<D>* the pointer
			 */
			const Edge<D> *operator->() const
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
		 * @brief Returns an iterator with the lowest Edge value
		 *
		 * @return Iterator the iterator
		 */
		Iterator begin()
		{
			return Iterator(Edge<D>(0));
		}
		/**
		 * @brief Returns an iterator with Edge<D>::null()
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
	 * @brief Get the edge on the opposite side
	 *
	 * @return Edge<D> the edge on the opposite side
	 */
	Edge<D> opposite() const
	{
		return Edge<D>(val ^ 0b11);
	}
	/**
	 * @brief Get the integer value of the edge.
	 *
	 * @return The integer value.
	 */
	size_t getIndex() const
	{
		return val;
	}
	/**
	 * @brief Get the axis that the edge lies along
	 *
	 * @return The axis
	 */
	size_t getAxisIndex() const
	{
		return val / 4;
	}
	/**
	 * @brief Get the Sides that the edge lies along
	 *
	 * @return std::array<Side<D>, 2>  the sides that the edge lies along
	 */
	std::array<Side<D>, 2> getSides() const
	{
		if (val == bs().getIndex()) {
			return {Side<D>::south(), Side<D>::bottom()};
		} else if (val == tn().getIndex()) {
			return {Side<D>::north(), Side<D>::top()};
		} else if (val == bn().getIndex()) {
			return {Side<D>::north(), Side<D>::bottom()};
		} else if (val == ts().getIndex()) {
			return {Side<D>::south(), Side<D>::top()};
		} else if (val == bw().getIndex()) {
			return {Side<D>::west(), Side<D>::bottom()};
		} else if (val == te().getIndex()) {
			return {Side<D>::east(), Side<D>::top()};
		} else if (val == be().getIndex()) {
			return {Side<D>::east(), Side<D>::bottom()};
		} else if (val == tw().getIndex()) {
			return {Side<D>::west(), Side<D>::top()};
		} else if (val == sw().getIndex()) {
			return {Side<D>::west(), Side<D>::south()};
		} else if (val == ne().getIndex()) {
			return {Side<D>::east(), Side<D>::north()};
		} else if (val == se().getIndex()) {
			return {Side<D>::east(), Side<D>::south()};
		} else if (val == nw().getIndex()) {
			return {Side<D>::west(), Side<D>::north()};
		} else {
			return {Side<D>::null(), Side<D>::null()};
		}
	}
	/**
	 * @brief Equals operator.
	 *
	 * @param other The other edge.
	 *
	 * @return Whether or not the value of this edge equals the value other edge.
	 */
	bool operator==(const Edge<D> &other) const
	{
		return val == other.val;
	}
	/**
	 * @brief Not Equals operator.
	 *
	 * @param other The other edge.
	 *
	 * @return Whether or not the value of this edge is not equal the value other edge.
	 */
	bool operator!=(const Edge<D> &other) const
	{
		return val != other.val;
	}
	/**
	 * @brief Less Than operator.
	 *
	 * @param other The other edge.
	 *
	 * @return Whether or not the value of this edge is less than the value other edge.
	 */
	bool operator<(const Edge<D> &other) const
	{
		return val < other.val;
	}
};

/**
 * @brief ostream operator that prints a string representation of Edge<1> enum.
 *
 * For example, Edge<1>::lower() will print out "Edge<1>::lower()".
 *
 * @param os the ostream
 * @param o the edge
 *
 * @return  the ostream
 */
inline std::ostream &operator<<(std::ostream &os, const Edge<1> &o)
{
	if (o == Edge<1>::null()) {
		os << "Edge<1>::null()";
	} else {
		os << "Edge<1> invalid value: " << o.getIndex();
	}
	return os;
}
/**
 * @brief ostream operator that prints a string representation of quadrant enum.
 *
 * For example, Edge<2>::sw() will print out "Edge<2>::sw()".
 *
 * @param os the ostream
 * @param o the edge
 *
 * @return  the ostream
 */
inline std::ostream &operator<<(std::ostream &os, const Edge<2> &o)
{
	if (o == Edge<2>::null()) {
		os << "Edge<2>::null()";
	} else {
		os << "Edge<2> invalid value: " << o.getIndex();
	}
	return os;
}
/**
 * @brief ostream operator that prints a string representation of edge enum.
 *
 * For example, Edge<3>::bsw() will print out "Edge<3>::bsw()".
 *
 * @param os the ostream
 * @param o the edge
 *
 * @return  the ostream
 */
inline std::ostream &operator<<(std::ostream &os, const Edge<3> &o)
{
	if (o == Edge<3>::bs()) {
		os << "Edge<3>::bs()";
	} else if (o == Edge<3>::bn()) {
		os << "Edge<3>::bn()";
	} else if (o == Edge<3>::ts()) {
		os << "Edge<3>::ts()";
	} else if (o == Edge<3>::tw()) {
		os << "Edge<3>::tw()";
	} else if (o == Edge<3>::bw()) {
		os << "Edge<3>::bw()";
	} else if (o == Edge<3>::be()) {
		os << "Edge<3>::be()";
	} else if (o == Edge<3>::tw()) {
		os << "Edge<3>::tw()";
	} else if (o == Edge<3>::te()) {
		os << "Edge<3>::te()";
	} else if (o == Edge<3>::sw()) {
		os << "Edge<3>::sw()";
	} else if (o == Edge<3>::se()) {
		os << "Edge<3>::se()";
	} else if (o == Edge<3>::nw()) {
		os << "Edge<3>::nw()";
	} else if (o == Edge<3>::ne()) {
		os << "Edge<3>::ne()";
	} else if (o == Edge<3>::null()) {
		os << "Edge<3>::null()";
	} else {
		os << "Edge<3> invalid value: " << o.getIndex();
	}
	return os;
}
void to_json(nlohmann::json &j, const Edge<1> &o);
void to_json(nlohmann::json &j, const Edge<2> &o);
void to_json(nlohmann::json &j, const Edge<3> &o);
void from_json(const nlohmann::json &j, Edge<1> &o);
void from_json(const nlohmann::json &j, Edge<2> &o);
void from_json(const nlohmann::json &j, Edge<3> &o);
} // namespace ThunderEgg
#endif
