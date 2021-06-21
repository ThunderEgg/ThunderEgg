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

#include "Communicator.h"

namespace ThunderEgg
{
namespace
{
void CheckErr(int err)
{
	if (err != MPI_SUCCESS) {
		std::string message = "MPI Call failed with error: ";
		char        err_string[MPI_MAX_ERROR_STRING];
		int         err_string_length;
		MPI_Error_string(err, err_string, &err_string_length);
		message += err_string;
		throw RuntimeError(message);
	}
}
} // namespace
Communicator::Communicator(const Communicator &other)
{
	if (other.comm != MPI_COMM_NULL) {
		CheckErr(MPI_Comm_dup(other.comm, &this->comm));
	}
}
Communicator &Communicator::operator=(const Communicator &other)
{
	if (other.comm != MPI_COMM_NULL) {
		CheckErr(MPI_Comm_dup(other.comm, &this->comm));
	}
	return *this;
}
Communicator::Communicator(MPI_Comm comm)
{
	CheckErr(MPI_Comm_dup(comm, &this->comm));
}
Communicator::~Communicator()
{
	if (comm != MPI_COMM_NULL) {
		MPI_Comm_free(&comm);
	}
}
MPI_Comm Communicator::getMPIComm() const
{
	if (comm == MPI_COMM_NULL) {
		throw RuntimeError("Null communicator");
	}
	return comm;
}
int Communicator::getRank() const
{
	int rank;
	CheckErr(MPI_Comm_rank(comm, &rank));
	return rank;
}
int Communicator::getSize() const
{
	int size;
	CheckErr(MPI_Comm_size(comm, &size));
	return size;
}
}; // namespace ThunderEgg