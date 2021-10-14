/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2019-2021 Scott Aiton
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
#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>
#include <mpi.h>
#if TEST_P4EST
#include <sc.h>
#endif
#if TEST_PETSC
#include <petscsys.h>
#endif

int
main(int argc, char* argv[])
{
#if TEST_P4EST
  sc_set_log_defaults(NULL, NULL, SC_LP_SILENT);
#endif
  // global setup...
#if TEST_PETSC
  PetscInitialize(nullptr, nullptr, nullptr, nullptr);
#else
  MPI_Init(nullptr, nullptr);
#endif

  int result = Catch::Session().run(argc, argv);

  // abort if failure, some tests can hang otherwise
  if (result > 0) {
    MPI_Abort(MPI_COMM_WORLD, result);
  }

  // global clean-up...
#if TEST_PETSC
  PetscFinalize();
#else
  MPI_Finalize();
#endif
  return result;
}
