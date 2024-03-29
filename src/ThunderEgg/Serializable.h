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

#ifndef SERIALIZABLE_H
#define SERIALIZABLE_H
/**
 * @file
 *
 * @brief Serializable class
 */

#include <memory>
namespace ThunderEgg {
/**
 * @brief Interface for serializing objects
 */
class Serializable
{
public:
  /**
   * @brief Destroy the Serializable object
   */
  virtual ~Serializable() = default;
  /**
   * @brief Serialize object into buffer
   *
   * @param buffer the buffer. Can be set to nullptr if you just want the size
   * @return int the size of the serialized object
   */
  virtual int serialize(char* buffer) const = 0;
  /**
   * @brief Deserialize an object
   *
   * @param buffer the buffer
   * @return int the size of object that was deserialized
   */
  virtual int deserialize(char* buffer) = 0;
};
} // namespace ThunderEgg
#endif
