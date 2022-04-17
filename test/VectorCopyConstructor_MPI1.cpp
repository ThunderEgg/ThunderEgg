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

#include "utils/DomainReader.h"
#include <ThunderEgg/Vector.h>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

#define MESHES "mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"

TEST_CASE("Vector<2> copy from domain constructor")
{
  for (auto mesh_file : { MESHES }) {
    for (int num_components : { 1, 2, 3 }) {
      for (auto num_ghost_cells : { 0, 1, 5 }) {
        for (int nx : { 1, 4, 5 }) {
          for (int ny : { 1, 4, 5 }) {
            Communicator comm(MPI_COMM_WORLD);

            int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
            int patch_stride = component_stride * num_components;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost_cells);
            Domain<2> domain = domain_reader.getFinerDomain();

            INFO("num_ghost_cells:   " << num_ghost_cells);
            INFO("nx:                " << nx);
            INFO("ny:                " << ny);

            Vector<2> vec_to_copy(domain, num_components);

            double* view_to_copy = &vec_to_copy.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_components * domain.getNumLocalPatches();
            for (size_t i = 0; i < size; i++) {
              double x = (i + 0.5) / size;
              view_to_copy[i] = 10 - (x - 0.75) * (x - 0.75);
            }

            Vector<2> vec(vec_to_copy);

            CHECK_EQ(vec.getNumGhostCells(), num_ghost_cells);

            int result;
            int err = MPI_Comm_compare(vec.getCommunicator().getMPIComm(), comm.getMPIComm(), &result);
            REQUIRE_EQ(err, MPI_SUCCESS);
            CHECK_EQ(result, MPI_CONGRUENT);

            CHECK_EQ(vec.getNumLocalPatches(), vec_to_copy.getNumLocalPatches());
            CHECK_EQ(vec.getNumComponents(), vec_to_copy.getNumComponents());
            CHECK_EQ(vec.getNumLocalCells(), vec_to_copy.getNumLocalCells());

            double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            CHECK_NE(view, view_to_copy);
            for (int i = 0; i < domain.getNumLocalPatches(); i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<double, 2> ld = vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }

            const Vector<2>& const_vec = vec;

            for (int i = 0; i < domain.getNumLocalPatches(); i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<const double, 2> ld = const_vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }
            for (int i = 0; i < domain.getNumLocalPatches(); i++) {
              PatchView<double, 2> view = vec.getPatchView(i);
              PatchView<double, 2> view_to_copy = vec_to_copy.getPatchView(i);
              Loop::OverAllIndexes<3>(view, [&](const std::array<int, 3> coord) { CHECK_EQ(view[coord], view_to_copy[coord]); });
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> copy from managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 4, 5 }) {
            Communicator comm(MPI_COMM_WORLD);

            int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
            int patch_stride = component_stride * num_components;

            INFO("num_ghost_cells:   " << num_ghost_cells);
            INFO("nx:                " << nx);
            INFO("ny:                " << ny);

            Vector<2> vec_to_copy(comm, { nx, ny }, num_components, num_local_patches, num_ghost_cells);

            double* view_to_copy = &vec_to_copy.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_components * num_local_patches;
            for (size_t i = 0; i < size; i++) {
              double x = (i + 0.5) / size;
              view_to_copy[i] = 10 - (x - 0.75) * (x - 0.75);
            }

            Vector<2> vec(vec_to_copy);

            CHECK_EQ(vec.getNumGhostCells(), num_ghost_cells);

            int result;
            int err = MPI_Comm_compare(vec.getCommunicator().getMPIComm(), comm.getMPIComm(), &result);
            REQUIRE_EQ(err, MPI_SUCCESS);
            CHECK_EQ(result, MPI_CONGRUENT);

            CHECK_EQ(vec.getNumLocalPatches(), vec_to_copy.getNumLocalPatches());
            CHECK_EQ(vec.getNumComponents(), vec_to_copy.getNumComponents());
            CHECK_EQ(vec.getNumLocalCells(), vec_to_copy.getNumLocalCells());

            double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            CHECK_NE(view, view_to_copy);
            for (int i = 0; i < num_local_patches; i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<double, 2> ld = vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }

            const Vector<2>& const_vec = vec;

            for (int i = 0; i < num_local_patches; i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<const double, 2> ld = const_vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }
            for (int i = 0; i < num_local_patches; i++) {
              PatchView<double, 2> view = vec.getPatchView(i);
              PatchView<double, 2> view_to_copy = vec_to_copy.getPatchView(i);
              Loop::OverAllIndexes<3>(view, [&](const std::array<int, 3> coord) { CHECK_EQ(view[coord], view_to_copy[coord]); });
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> copy from unmanaged constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 4, 5 }) {
            Communicator comm(MPI_COMM_WORLD);

            int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
            int patch_stride = component_stride * num_components;

            INFO("num_ghost_cells:   " << num_ghost_cells);
            INFO("nx:                " << nx);
            INFO("ny:                " << ny);

            array<int, 3> lengths = { nx, ny, num_components };
            array<int, 3> strides = { 1, nx + 2 * num_ghost_cells, (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) };
            double data[patch_stride * num_local_patches];
            vector<double*> patch_starts(num_local_patches);
            for (int i = 0; i < num_local_patches; i++) {
              patch_starts[i] = data + i * patch_stride;
            }
            Vector<2> vec_to_copy(comm, patch_starts, strides, lengths, num_ghost_cells);

            double* view_to_copy = &vec_to_copy.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_components * num_local_patches;
            for (size_t i = 0; i < size; i++) {
              double x = (i + 0.5) / size;
              view_to_copy[i] = 10 - (x - 0.75) * (x - 0.75);
            }

            Vector<2> vec(vec_to_copy);

            CHECK_EQ(vec.getNumGhostCells(), num_ghost_cells);

            int result;
            int err = MPI_Comm_compare(vec.getCommunicator().getMPIComm(), comm.getMPIComm(), &result);
            REQUIRE_EQ(err, MPI_SUCCESS);
            CHECK_EQ(result, MPI_CONGRUENT);

            CHECK_EQ(vec.getNumLocalPatches(), vec_to_copy.getNumLocalPatches());
            CHECK_EQ(vec.getNumComponents(), vec_to_copy.getNumComponents());
            CHECK_EQ(vec.getNumLocalCells(), vec_to_copy.getNumLocalCells());

            double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            CHECK_NE(view, view_to_copy);
            for (int i = 0; i < num_local_patches; i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<double, 2> ld = vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }

            const Vector<2>& const_vec = vec;

            for (int i = 0; i < num_local_patches; i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<const double, 2> ld = const_vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }
            for (int i = 0; i < num_local_patches; i++) {
              PatchView<double, 2> view = vec.getPatchView(i);
              PatchView<double, 2> view_to_copy = vec_to_copy.getPatchView(i);
              Loop::OverAllIndexes<3>(view, [&](const std::array<int, 3> coord) { CHECK_EQ(view[coord], view_to_copy[coord]); });
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> copy assign from domain constructor")
{
  for (auto mesh_file : { MESHES }) {
    for (int num_components : { 1, 2, 3 }) {
      for (auto num_ghost_cells : { 0, 1, 5 }) {
        for (int nx : { 1, 4, 5 }) {
          for (int ny : { 1, 4, 5 }) {
            Communicator comm(MPI_COMM_WORLD);

            int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
            int patch_stride = component_stride * num_components;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost_cells);
            Domain<2> domain = domain_reader.getFinerDomain();

            INFO("num_ghost_cells:   " << num_ghost_cells);
            INFO("nx:                " << nx);
            INFO("ny:                " << ny);

            Vector<2> vec_to_copy(domain, num_components);

            double* view_to_copy = &vec_to_copy.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_components * domain.getNumLocalPatches();
            for (size_t i = 0; i < size; i++) {
              double x = (i + 0.5) / size;
              view_to_copy[i] = 10 - (x - 0.75) * (x - 0.75);
            }

            Vector<2> vec;
            vec = vec_to_copy;

            CHECK_EQ(vec.getNumGhostCells(), num_ghost_cells);

            int result;
            int err = MPI_Comm_compare(vec.getCommunicator().getMPIComm(), comm.getMPIComm(), &result);
            REQUIRE_EQ(err, MPI_SUCCESS);
            CHECK_EQ(result, MPI_CONGRUENT);

            CHECK_EQ(vec.getNumLocalPatches(), vec_to_copy.getNumLocalPatches());
            CHECK_EQ(vec.getNumComponents(), vec_to_copy.getNumComponents());
            CHECK_EQ(vec.getNumLocalCells(), vec_to_copy.getNumLocalCells());

            double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            CHECK_NE(view, view_to_copy);
            for (int i = 0; i < domain.getNumLocalPatches(); i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<double, 2> ld = vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }

            const Vector<2>& const_vec = vec;

            for (int i = 0; i < domain.getNumLocalPatches(); i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<const double, 2> ld = const_vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }
            for (int i = 0; i < domain.getNumLocalPatches(); i++) {
              PatchView<double, 2> view = vec.getPatchView(i);
              PatchView<double, 2> view_to_copy = vec_to_copy.getPatchView(i);
              Loop::OverAllIndexes<3>(view, [&](const std::array<int, 3> coord) { CHECK_EQ(view[coord], view_to_copy[coord]); });
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> copy assign from managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 4, 5 }) {
            Communicator comm(MPI_COMM_WORLD);

            int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
            int patch_stride = component_stride * num_components;

            INFO("num_ghost_cells:   " << num_ghost_cells);
            INFO("nx:                " << nx);
            INFO("ny:                " << ny);

            Vector<2> vec_to_copy(comm, { nx, ny }, num_components, num_local_patches, num_ghost_cells);

            double* view_to_copy = &vec_to_copy.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_components * num_local_patches;
            for (size_t i = 0; i < size; i++) {
              double x = (i + 0.5) / size;
              view_to_copy[i] = 10 - (x - 0.75) * (x - 0.75);
            }

            Vector<2> vec;
            vec = vec_to_copy;

            CHECK_EQ(vec.getNumGhostCells(), num_ghost_cells);

            int result;
            int err = MPI_Comm_compare(vec.getCommunicator().getMPIComm(), comm.getMPIComm(), &result);
            REQUIRE_EQ(err, MPI_SUCCESS);
            CHECK_EQ(result, MPI_CONGRUENT);

            CHECK_EQ(vec.getNumLocalPatches(), vec_to_copy.getNumLocalPatches());
            CHECK_EQ(vec.getNumComponents(), vec_to_copy.getNumComponents());
            CHECK_EQ(vec.getNumLocalCells(), vec_to_copy.getNumLocalCells());

            double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            CHECK_NE(view, view_to_copy);
            for (int i = 0; i < num_local_patches; i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<double, 2> ld = vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }

            const Vector<2>& const_vec = vec;

            for (int i = 0; i < num_local_patches; i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<const double, 2> ld = const_vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }
            for (int i = 0; i < num_local_patches; i++) {
              PatchView<double, 2> view = vec.getPatchView(i);
              PatchView<double, 2> view_to_copy = vec_to_copy.getPatchView(i);
              Loop::OverAllIndexes<3>(view, [&](const std::array<int, 3> coord) { CHECK_EQ(view[coord], view_to_copy[coord]); });
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> copy assign from unmanaged constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 4, 5 }) {
            Communicator comm(MPI_COMM_WORLD);

            int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
            int patch_stride = component_stride * num_components;

            INFO("num_ghost_cells:   " << num_ghost_cells);
            INFO("nx:                " << nx);
            INFO("ny:                " << ny);

            array<int, 3> lengths = { nx, ny, num_components };
            array<int, 3> strides = { 1, nx + 2 * num_ghost_cells, (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) };
            double data[patch_stride * num_local_patches];
            vector<double*> patch_starts(num_local_patches);
            for (int i = 0; i < num_local_patches; i++) {
              patch_starts[i] = data + i * patch_stride;
            }
            Vector<2> vec_to_copy(comm, patch_starts, strides, lengths, num_ghost_cells);

            double* view_to_copy = &vec_to_copy.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            int size = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * num_components * num_local_patches;
            for (size_t i = 0; i < size; i++) {
              double x = (i + 0.5) / size;
              view_to_copy[i] = 10 - (x - 0.75) * (x - 0.75);
            }

            Vector<2> vec;
            vec = vec_to_copy;

            CHECK_EQ(vec.getNumGhostCells(), num_ghost_cells);

            int result;
            int err = MPI_Comm_compare(vec.getCommunicator().getMPIComm(), comm.getMPIComm(), &result);
            REQUIRE_EQ(err, MPI_SUCCESS);
            CHECK_EQ(result, MPI_CONGRUENT);

            CHECK_EQ(vec.getNumLocalPatches(), vec_to_copy.getNumLocalPatches());
            CHECK_EQ(vec.getNumComponents(), vec_to_copy.getNumComponents());
            CHECK_EQ(vec.getNumLocalCells(), vec_to_copy.getNumLocalCells());

            double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            CHECK_NE(view, view_to_copy);
            for (int i = 0; i < num_local_patches; i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<double, 2> ld = vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }

            const Vector<2>& const_vec = vec;

            for (int i = 0; i < num_local_patches; i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<const double, 2> ld = const_vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }
            for (int i = 0; i < num_local_patches; i++) {
              PatchView<double, 2> view = vec.getPatchView(i);
              PatchView<double, 2> view_to_copy = vec_to_copy.getPatchView(i);
              Loop::OverAllIndexes<3>(view, [&](const std::array<int, 3> coord) { CHECK_EQ(view[coord], view_to_copy[coord]); });
            }
          }
        }
      }
    }
  }
}
