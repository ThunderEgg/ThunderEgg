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

#ifndef THUNDEREGG_COMMUNICATOR_H
#define THUNDEREGG_COMMUNICATOR_H
/**
 * @file
 *
 * @brief Communicator class
 */
#include <ThunderEgg/RuntimeError.h>
#include <mpi.h>
#include <memory>

namespace ThunderEgg {
/**
 * @brief wrapper arount MPI_Comm, provides proper copy operators. Classes that have a communicator
 * are meant to store a Communicator object instead of a raw MPI_Comm
 *
 * This class will only call MPI_Comm_dup with the MPI_Comm constructor.
 * When sharing the commmunicator whithin the ThunderEgg library,
 * MPI_Comm_dup will not be called, and they will share the same MPI_Comm.
 *
 */
class Communicator
{
private:
  /**
   * @brief The communicator
   */
  std::shared_ptr<MPI_Comm> comm;

public:
  /**
   * @brief Construct a new Communicator with a null communicator
   */
  Communicator() = default;
  /**
   * @brief Destroy the Communicator object
   */
  ~Communicator() = default;
  /**
   * @brief Construct a new Communicator from a specified MPI_Comm
   *
   * This will call MPI_Comm_dup on the communicator.
   *
   * @param comm the comm
   */
  explicit Communicator(MPI_Comm comm);
  /**
   * @brief Get the raw MPI_Comm object
   *
   * @return MPI_Comm the comm
   */
  MPI_Comm getMPIComm() const;
  /**
   * @brief Get the size of the communicator
   *
   * @return int the size
   */
  int getSize() const;
  /**
   * @brief Get the rank of this processor
   *
   * @return int the rank
   */
  int getRank() const;
};
} // namespace ThunderEgg
#endif
