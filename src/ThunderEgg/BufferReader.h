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

#ifndef THUNDEREGG_BUFFERREADER_H
#define THUNDEREGG_BUFFERREADER_H
/**
 * @file
 *
 * @brief BufferReader class
 */

#include <ThunderEgg/Serializable.h>
#include <cstddef>
#include <iostream>
#include <type_traits>

namespace ThunderEgg
{
/**
 * @brief Class that is used to help read serialized objects from a buffer.
 */
class BufferReader
{
	private:
	/**
	 * @brief The pointer to the buffer
	 */
	char *buffer = nullptr;
	/**
	 * @brief The current position in the buffer
	 */
	int pos = 0;
	/**
	 * @brief return if T is derived from Serializable
	 *
	 * @tparam T the type
	 * @return true if derived from Serializable
	 * @return false if not derived from Serializable
	 */
	template <typename T> static constexpr bool isSerializable()
	{
		return std::is_base_of<Serializable, T>::value;
	}

	public:
	/**
	 * @brief Create a new BufferReader with given buffer.
	 *
	 * @param buffer the pointer to the beginning of the buffer.
	 */
	BufferReader(char *buffer)
	{
		this->buffer = buffer;
	}
	/**
	 * @brief get the current position in the buffer
	 *
	 * @return the current position
	 */
	int getPos()
	{
		return pos;
	}
	/**
	 * @brief Get an object of the buffer.
	 *
	 * @param obj the Serializable object.
	 *
	 * @return  this BufferReader
	 */
	BufferReader &operator>>(Serializable &obj)
	{
		pos += obj.deserialize(buffer + pos);
		return *this;
	}
	/**
	 * @brief Get an object from the buffer.
	 *
	 * @tparam T the type of the object.
	 * @param obj the object. This object must be in serialized form.
	 *
	 * @return  this BufferReader
	 */
	template <typename T> typename std::enable_if<!isSerializable<T>(), BufferReader>::type &operator>>(T &obj)
	{
		for (size_t i = 0; i < sizeof(T); i++) {
			reinterpret_cast<char *>(&obj)[i] = buffer[pos];
			pos++;
		}
		return *this;
	}
};
} // namespace ThunderEgg
#endif
