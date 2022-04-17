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

#include <mpi.h>
#if TEST_P4EST
#include <sc.h>
#endif
#if TEST_PETSC
#include <petscsys.h>
#endif

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>
int
main(int argc, char* argv[])
{
  bool catch_add_tests = false;
  for (int i = 0; i < argc; i++) {
    catch_add_tests = strcmp(argv[i], "--list-tests") == 0;
    if (catch_add_tests)
      break;
    catch_add_tests = strcmp(argv[i], "--list-reporters") == 0;
    if (catch_add_tests)
      break;
  }
#if TEST_P4EST
  sc_set_log_defaults(NULL, NULL, SC_LP_SILENT);
#endif
  // global setup...
#if TEST_PETSC
  PetscInitialize(nullptr, nullptr, nullptr, nullptr);
#else
  MPI_Init(nullptr, nullptr);
#endif
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int result = 0;
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  if (!catch_add_tests || rank == 0) {
    result = context.run();
  }

  // abort if failure, some tests can hang otherwise
  if (!catch_add_tests && (result > 0 || context.shouldExit())) {
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
