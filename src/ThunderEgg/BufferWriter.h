/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2018-2021 Scott Aiton
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

#ifndef THUNDEREGG_BUFFERWRITER_H
#define THUNDEREGG_BUFFERWRITER_H
/**
 * @file
 *
 * @brief BufferWriter class
 */

#include <ThunderEgg/Serializable.h>
#include <cstddef>
#include <iostream>
#include <type_traits>

namespace ThunderEgg {
/**
 * @brief Class that is used to help serialize objects into a buffer.
 */
class BufferWriter
{
private:
  /**
   * @brief The pointer to the buffer
   */
  char* buffer = nullptr;
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
  template<typename T>
  static constexpr bool isSerializable()
  {
    return std::is_base_of<Serializable, T>::value;
  }

public:
  /**
   * @brief Create a new BufferWriter with the buffer set to nullptr. This is helpful for
   * determining the size needed for the buffer.
   */
  BufferWriter() = default;
  /**
   * @brief Create a new BufferWriter with given buffer.
   *
   * @param buffer the pointer to the beginning of the buffer.
   */
  BufferWriter(char* buffer) { this->buffer = buffer; }
  /**
   * @brief get the current position in the buffer
   *
   * @return  the current position
   */
  int getPos() { return pos; }
  /**
   * @brief  Add object to the buffer.
   *
   * @param obj the Serializable object.
   *
   * @return  this BufferWriter
   */
  BufferWriter& operator<<(const Serializable& obj)
  {
    pos += obj.serialize(buffer == nullptr ? nullptr : (buffer + pos));
    return *this;
  }
  /**
   * @brief Add an object to the buffer.
   *
   * @tparam T the type of the object.
   * @param obj the object. This object must be in serialized form.
   *
   * @return  this BufferWriter
   */
  template<typename T>
  typename std::enable_if<!isSerializable<T>(), BufferWriter>::type& operator<<(const T& obj)
  {
    if (buffer != nullptr) {
      for (size_t i = 0; i < sizeof(T); i++) {
        buffer[pos] = reinterpret_cast<const char*>(&obj)[i];
        pos++;
      }
    } else {
      pos += sizeof(T);
    }
    return *this;
  }
};
} // namespace ThunderEgg
#endif
