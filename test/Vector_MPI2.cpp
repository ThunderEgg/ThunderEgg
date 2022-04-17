/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
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

#include "Vector_MOCKS.h"
#include "utils/DomainReader.h"

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Vector<3> twoNorm")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int nz : { 1, 4, 5 }) {
            for (int num_local_patches : { 1, 13 }) {
              Communicator comm(MPI_COMM_WORLD);
              array<int, 3> ns = { nx, ny, nz };

              Vector<3> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

              INFO("num_ghost_cells:   " << num_ghost_cells);
              INFO("nx:                " << nx);
              INFO("ny:                " << ny);
              INFO("nz:                " << nz);
              INFO("num_local_patches: " << num_local_patches);
              INFO("num_components:    " << num_components);

              int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells) * num_components * num_local_patches;
              double* vec_data = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                vec_data[i] = 10 - (x - 0.75) * (x - 0.75);
              }
              vec.setWithGhost(1);
              vec.shift(28);

              double expected_norm = 0;
              for (int i = 0; i < vec.getNumLocalPatches(); i++) {
                for (int c = 0; c < vec.getNumComponents(); c++) {
                  auto ld = vec.getComponentView(c, i);
                  Loop::Nested<3>(ld.getStart(), ld.getEnd(), [&](std::array<int, 3>& coord) { expected_norm += ld[coord] * ld[coord]; });
                }
              }
              double global_expected_norm;
              MPI_Allreduce(&expected_norm, &global_expected_norm, 1, MPI_DOUBLE, MPI_SUM, comm.getMPIComm());
              global_expected_norm = sqrt(global_expected_norm);

              CHECK(vec.twoNorm() == doctest::Approx(global_expected_norm));
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> infNorm")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int nz : { 1, 4, 5 }) {
            for (int num_local_patches : { 1, 13 }) {
              Communicator comm(MPI_COMM_WORLD);
              array<int, 3> ns = { nx, ny, nz };

              Vector<3> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

              INFO("num_ghost_cells:   " << num_ghost_cells);
              INFO("nx:                " << nx);
              INFO("ny:                " << ny);
              INFO("nz:                " << nz);
              INFO("num_local_patches: " << num_local_patches);
              INFO("num_components:    " << num_components);

              int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells) * num_components * num_local_patches;
              double* vec_data = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                vec_data[i] = 10 - (x - 0.75) * (x - 0.75);
              }
              vec.setWithGhost(1);
              vec.shift(28);

              double expected_norm = 0;
              for (int i = 0; i < vec.getNumLocalPatches(); i++) {
                for (int c = 0; c < vec.getNumComponents(); c++) {
                  auto ld = vec.getComponentView(c, i);
                  Loop::Nested<3>(ld.getStart(), ld.getEnd(), [&](std::array<int, 3>& coord) { expected_norm = max(abs(ld[coord]), expected_norm); });
                }
              }
              double global_expected_norm;
              MPI_Allreduce(&expected_norm, &global_expected_norm, 1, MPI_DOUBLE, MPI_MAX, comm.getMPIComm());

              CHECK(vec.infNorm() == global_expected_norm);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> dot")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int nz : { 1, 4, 5 }) {
            for (int num_local_patches : { 1, 13 }) {
              Communicator comm(MPI_COMM_WORLD);
              array<int, 3> ns = { nx, ny, nz };

              Vector<3> a(comm, ns, num_components, num_local_patches, num_ghost_cells);
              Vector<3> b(comm, ns, num_components, num_local_patches, num_ghost_cells);

              INFO("num_ghost_cells:   " << num_ghost_cells);
              INFO("nx:                " << nx);
              INFO("ny:                " << ny);
              INFO("nz:                " << nz);
              INFO("num_local_patches: " << num_local_patches);
              INFO("num_components:    " << num_components);

              int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells) * num_components * num_local_patches;
              double* a_data = &a.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                a_data[i] = 10 - (x - 0.75) * (x - 0.75);
              }

              double* b_data = &b.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                b_data[i] = (x - 0.5) * (x - 0.5);
              }

              double expected_value = 0;
              for (int i = 0; i < a.getNumLocalPatches(); i++) {
                for (int c = 0; c < a.getNumComponents(); c++) {
                  auto a_ld = a.getComponentView(c, i);
                  auto b_ld = b.getComponentView(c, i);
                  Loop::Nested<3>(b_ld.getStart(), b_ld.getEnd(), [&](std::array<int, 3>& coord) { expected_value += b_ld[coord] * a_ld[coord]; });
                }
              }
              double global_expected_value;
              MPI_Allreduce(&expected_value, &global_expected_value, 1, MPI_DOUBLE, MPI_SUM, comm.getMPIComm());

              CHECK(a.dot(b) == doctest::Approx(global_expected_value));
            }
          }
        }
      }
    }
  }
}
