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

TEST_CASE("Vector<2> getNumGhostCells domain constructor")
{
  for (auto mesh_file : { MESHES }) {
    for (int num_components : { 1, 2, 3 }) {
      for (auto num_ghost_cells : { 0, 1, 5 }) {
        for (int nx : { 1, 4, 5 }) {
          for (int ny : { 1, 4, 5 }) {
            Communicator comm(MPI_COMM_WORLD);

            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost_cells);
            Domain<2> domain = domain_reader.getFinerDomain();

            INFO("num_ghost_cells:   " << num_ghost_cells);
            INFO("nx:                " << nx);
            INFO("ny:                " << ny);

            Vector<2> vec(domain, num_components);

            CHECK(vec.getNumGhostCells() == num_ghost_cells);
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getMPIComm domain constructor")
{
  for (auto mesh_file : { MESHES }) {
    for (int num_components : { 1, 2, 3 }) {
      for (auto num_ghost_cells : { 0, 1, 5 }) {
        for (int nx : { 1, 4, 5 }) {
          for (int ny : { 1, 4, 5 }) {
            Communicator comm(MPI_COMM_WORLD);

            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost_cells);
            Domain<2> domain = domain_reader.getFinerDomain();

            INFO("num_ghost_cells:   " << num_ghost_cells);
            INFO("nx:                " << nx);
            INFO("ny:                " << ny);

            Vector<2> vec(domain, num_components);

            int result;
            int err = MPI_Comm_compare(vec.getCommunicator().getMPIComm(), comm.getMPIComm(), &result);
            REQUIRE(err == MPI_SUCCESS);
            CHECK(result == MPI_CONGRUENT);
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getNumLocalPatches domain constructor")
{
  for (auto mesh_file : { MESHES }) {
    for (int num_components : { 1, 2, 3 }) {
      for (auto num_ghost_cells : { 0, 1, 5 }) {
        for (int nx : { 1, 4, 5 }) {
          for (int ny : { 1, 4, 5 }) {
            Communicator comm(MPI_COMM_WORLD);

            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost_cells);
            Domain<2> domain = domain_reader.getFinerDomain();

            INFO("num_ghost_cells:   " << num_ghost_cells);
            INFO("nx:                " << nx);
            INFO("ny:                " << ny);

            Vector<2> vec(domain, num_components);

            CHECK(vec.getNumLocalPatches() == domain.getNumLocalPatches());
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getNumComponents domain constructor")
{
  for (auto mesh_file : { MESHES }) {
    for (int num_components : { 1, 2, 3 }) {
      for (auto num_ghost_cells : { 0, 1, 5 }) {
        for (int nx : { 1, 4, 5 }) {
          for (int ny : { 1, 4, 5 }) {
            Communicator comm(MPI_COMM_WORLD);

            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost_cells);
            Domain<2> domain = domain_reader.getFinerDomain();

            INFO("num_ghost_cells:   " << num_ghost_cells);
            INFO("nx:                " << nx);
            INFO("ny:                " << ny);

            Vector<2> vec(domain, num_components);

            CHECK(vec.getNumComponents() == num_components);
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getNumLocalCells domain constructor")
{
  for (auto mesh_file : { MESHES }) {
    for (int num_components : { 1, 2, 3 }) {
      for (auto num_ghost_cells : { 0, 1, 5 }) {
        for (int nx : { 1, 4, 5 }) {
          for (int ny : { 1, 4, 5 }) {
            Communicator comm(MPI_COMM_WORLD);

            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost_cells);
            Domain<2> domain = domain_reader.getFinerDomain();

            INFO("num_ghost_cells:   " << num_ghost_cells);
            INFO("nx:                " << nx);
            INFO("ny:                " << ny);

            Vector<2> vec(domain, num_components);

            CHECK(vec.getNumLocalCells() == domain.getNumLocalCells());
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getComponentView domain constructor")
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

            Vector<2> vec(domain, num_components);

            double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            for (int i = 0; i < domain.getNumLocalPatches(); i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<double, 2> ld = vec.getComponentView(c, i);
                CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i + c * component_stride);
                CHECK(&ld[ld.getGhostEnd()] == view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }
            if (ENABLE_DEBUG) {
              CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
              CHECK_THROWS_AS(vec.getComponentView(num_components, 0), RuntimeError);
              CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
              CHECK_THROWS_AS(vec.getComponentView(0, domain.getNumLocalPatches()), RuntimeError);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getComponentView const domain constructor")
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

            const Vector<2> vec(domain, num_components);

            const double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            for (int i = 0; i < domain.getNumLocalPatches(); i++) {
              INFO("i:                 " << i);
              for (int c = 0; c < num_components; c++) {
                INFO("c:                 " << c);
                ComponentView<const double, 2> ld = vec.getComponentView(c, i);
                CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i + c * component_stride);
                CHECK(&ld[ld.getGhostEnd()] == view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }
            if (ENABLE_DEBUG) {
              CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
              CHECK_THROWS_AS(vec.getComponentView(num_components, 0), RuntimeError);
              CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
              CHECK_THROWS_AS(vec.getComponentView(0, domain.getNumLocalPatches()), RuntimeError);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getPatchView domain constructor")
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

            Vector<2> vec(domain, num_components);

            double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            for (int i = 0; i < domain.getNumLocalPatches(); i++) {
              INFO("i:                 " << i);
              PatchView<double, 2> ld = vec.getPatchView(i);
              CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i);
              CHECK(&ld[ld.getGhostEnd()] == view + patch_stride * (i + 1) - 1);
            }
            if (ENABLE_DEBUG) {
              CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
              CHECK_THROWS_AS(vec.getPatchView(domain.getNumLocalPatches()), RuntimeError);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getPatchView const domain constructor")
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

            const Vector<2> vec(domain, num_components);

            const double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            for (int i = 0; i < domain.getNumLocalPatches(); i++) {
              INFO("i:                 " << i);
              PatchView<const double, 2> ld = vec.getPatchView(i);
              CHECK(&ld[ld.getGhostStart()] == view + patch_stride * i);
              CHECK(&ld[ld.getGhostEnd()] == view + patch_stride * (i + 1) - 1);
            }
            if (ENABLE_DEBUG) {
              CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
              CHECK_THROWS_AS(vec.getPatchView(domain.getNumLocalPatches()), RuntimeError);
            }
          }
        }
      }
    }
  }
}
