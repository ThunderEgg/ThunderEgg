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

#include "Communicator.h"
#include <utility>

namespace ThunderEgg {
namespace {
void
CheckErr(int err)
{
  if (err != MPI_SUCCESS) {
    std::string message = "MPI Call failed with error: ";
    char err_string[MPI_MAX_ERROR_STRING];
    int err_string_length;
    MPI_Error_string(err, err_string, &err_string_length);
    message += err_string;
    throw RuntimeError(message);
  }
}
} // namespace

Communicator::Communicator(MPI_Comm comm)
{
  MPI_Comm* comm_dup = new MPI_Comm;
  CheckErr(MPI_Comm_dup(comm, comm_dup));
  this->comm.reset(comm_dup, [](MPI_Comm* comm_ptr) {
    int finalized;
    MPI_Finalized(&finalized);
    if (*comm_ptr != MPI_COMM_NULL && !finalized) {
      MPI_Comm_free(comm_ptr);
    }
  });
}

MPI_Comm
Communicator::getMPIComm() const
{
  if (comm == nullptr) {
    throw RuntimeError("Null communicator");
  }
  return *comm;
}
int
Communicator::getRank() const
{
  if (comm == nullptr) {
    throw RuntimeError("Null communicator");
  }
  int rank;
  CheckErr(MPI_Comm_rank(*comm, &rank));
  return rank;
}
int
Communicator::getSize() const
{
  if (comm == nullptr) {
    throw RuntimeError("Null communicator");
  }
  int size;
  CheckErr(MPI_Comm_size(*comm, &size));
  return size;
}
}; // namespace ThunderEgg