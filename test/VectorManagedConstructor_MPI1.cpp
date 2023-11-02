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

#include <ThunderEgg/Vector.h>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Vector<1> getNumGhostCells managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int num_local_patches : { 1, 13 }) {
          Communicator comm(MPI_COMM_WORLD);

          array<int, 1> ns = { nx };

          Vector<1> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

          CHECK_EQ(vec.getNumGhostCells(), num_ghost_cells);
        }
      }
    }
  }
}
TEST_CASE("Vector<1> getMPIComm managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int num_local_patches : { 1, 13 }) {
          Communicator comm(MPI_COMM_WORLD);

          array<int, 1> ns = { nx };

          Vector<1> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

          int result;
          int err = MPI_Comm_compare(vec.getCommunicator().getMPIComm(), comm.getMPIComm(), &result);
          REQUIRE_EQ(err, MPI_SUCCESS);
          CHECK_EQ(result, MPI_IDENT);
        }
      }
    }
  }
}
TEST_CASE("Vector<1> getNumLocalPatches managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int num_local_patches : { 1, 13 }) {
          Communicator comm(MPI_COMM_WORLD);

          array<int, 1> ns = { nx };

          Vector<1> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

          CHECK_EQ(vec.getNumLocalPatches(), num_local_patches);
        }
      }
    }
  }
}
TEST_CASE("Vector<1> getNumComponents managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int num_local_patches : { 1, 13 }) {
          Communicator comm(MPI_COMM_WORLD);

          array<int, 1> ns = { nx };

          Vector<1> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

          CHECK_EQ(vec.getNumComponents(), num_components);
        }
      }
    }
  }
}
TEST_CASE("Vector<1> getNumLocalCells managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int num_local_patches : { 1, 13 }) {
          Communicator comm(MPI_COMM_WORLD);

          array<int, 1> ns = { nx };

          Vector<1> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

          CHECK_EQ(vec.getNumLocalCells(), nx * num_local_patches);
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getNumGhostCells managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 13 }) {
            Communicator comm(MPI_COMM_WORLD);

            array<int, 2> ns = { nx, ny };

            Vector<2> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

            CHECK_EQ(vec.getNumGhostCells(), num_ghost_cells);
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getMPIComm managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 13 }) {
            Communicator comm(MPI_COMM_WORLD);

            array<int, 2> ns = { nx, ny };

            Vector<2> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

            int result;
            int err = MPI_Comm_compare(vec.getCommunicator().getMPIComm(), comm.getMPIComm(), &result);
            REQUIRE_EQ(err, MPI_SUCCESS);
            CHECK_EQ(result, MPI_IDENT);
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getNumLocalPatches managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 13 }) {
            Communicator comm(MPI_COMM_WORLD);

            array<int, 2> ns = { nx, ny };

            Vector<2> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

            CHECK_EQ(vec.getNumLocalPatches(), num_local_patches);
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getNumComponents managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 13 }) {
            Communicator comm(MPI_COMM_WORLD);

            array<int, 2> ns = { nx, ny };

            Vector<2> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

            CHECK_EQ(vec.getNumComponents(), num_components);
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getNumLocalCells managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 13 }) {
            Communicator comm(MPI_COMM_WORLD);

            array<int, 2> ns = { nx, ny };

            Vector<2> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

            CHECK_EQ(vec.getNumLocalCells(), nx * ny * num_local_patches);
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getNumGhostCells managed constructor")
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

              CHECK_EQ(vec.getNumGhostCells(), num_ghost_cells);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getMPIComm managed constructor")
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

              int result;
              int err = MPI_Comm_compare(vec.getCommunicator().getMPIComm(), comm.getMPIComm(), &result);
              REQUIRE_EQ(err, MPI_SUCCESS);
              CHECK_EQ(result, MPI_IDENT);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getNumLocalPatches managed constructor")
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

              CHECK_EQ(vec.getNumLocalPatches(), num_local_patches);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getNumComponents managed constructor")
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

              CHECK_EQ(vec.getNumComponents(), num_components);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getNumLocalCells managed constructor")
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

              CHECK_EQ(vec.getNumLocalCells(), nx * ny * nz * num_local_patches);
            }
          }
        }
      }
    }
  }
}

TEST_CASE("Vector<1> getComponentView.h managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int num_local_patches : { 1, 13 }) {
          Communicator comm(MPI_COMM_WORLD);

          array<int, 1> ns = { nx };
          int component_stride = (nx + 2 * num_ghost_cells);
          int patch_stride = component_stride * num_components;

          Vector<1> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

          double* view = &vec.getPatchView(0)(-num_ghost_cells, 0);
          for (int i = 0; i < num_local_patches; i++) {
            for (int c = 0; c < num_components; c++) {
              ComponentView<double, 1> ld = vec.getComponentView(c, i);
              CHECK_EQ(&ld[ld.getGhostStart()], view + (patch_stride * i + c * component_stride));
              CHECK_EQ(&ld[ld.getGhostEnd()], view + (patch_stride * i + (c + 1) * component_stride) - 1);
            }
          }
          if (ENABLE_DEBUG) {
            CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
            CHECK_THROWS_AS(vec.getComponentView(num_components, 0), RuntimeError);
            CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
            CHECK_THROWS_AS(vec.getComponentView(0, num_local_patches), RuntimeError);
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<1> getComponentView const managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int num_local_patches : { 1, 13 }) {
          Communicator comm(MPI_COMM_WORLD);

          array<int, 1> ns = { nx };
          int component_stride = (nx + 2 * num_ghost_cells);
          int patch_stride = component_stride * num_components;

          const Vector<1> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

          const double* view = &vec.getPatchView(0)(-num_ghost_cells, 0);
          for (int i = 0; i < num_local_patches; i++) {
            for (int c = 0; c < num_components; c++) {
              ComponentView<const double, 1> ld = vec.getComponentView(c, i);
              CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
              CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
            }
          }
          if (ENABLE_DEBUG) {
            CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
            CHECK_THROWS_AS(vec.getComponentView(num_components, 0), RuntimeError);
            CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
            CHECK_THROWS_AS(vec.getComponentView(0, num_local_patches), RuntimeError);
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getComponentView managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 13 }) {
            Communicator comm(MPI_COMM_WORLD);

            array<int, 2> ns = { nx, ny };
            int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
            int patch_stride = component_stride * num_components;

            Vector<2> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

            double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            for (int i = 0; i < num_local_patches; i++) {
              for (int c = 0; c < num_components; c++) {
                ComponentView<double, 2> ld = vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }
            if (ENABLE_DEBUG) {
              CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
              CHECK_THROWS_AS(vec.getComponentView(num_components, 0), RuntimeError);
              CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
              CHECK_THROWS_AS(vec.getComponentView(0, num_local_patches), RuntimeError);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getComponentView const managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 13 }) {
            Communicator comm(MPI_COMM_WORLD);

            array<int, 2> ns = { nx, ny };
            int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
            int patch_stride = component_stride * num_components;

            const Vector<2> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

            const double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            for (int i = 0; i < num_local_patches; i++) {
              for (int c = 0; c < num_components; c++) {
                ComponentView<const double, 2> ld = vec.getComponentView(c, i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
              }
            }
            if (ENABLE_DEBUG) {
              CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
              CHECK_THROWS_AS(vec.getComponentView(num_components, 0), RuntimeError);
              CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
              CHECK_THROWS_AS(vec.getComponentView(0, num_local_patches), RuntimeError);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getComponentView managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int nz : { 1, 4, 5 }) {
            for (int num_local_patches : { 1, 13 }) {
              Communicator comm(MPI_COMM_WORLD);

              array<int, 3> ns = { nx, ny, nz };
              int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells);
              int patch_stride = component_stride * num_components;

              Vector<3> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

              double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (int i = 0; i < num_local_patches; i++) {
                for (int c = 0; c < num_components; c++) {
                  ComponentView<double, 3> ld = vec.getComponentView(c, i);
                  CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                  CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
                }
              }
              if (ENABLE_DEBUG) {
                CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
                CHECK_THROWS_AS(vec.getComponentView(num_components, 0), RuntimeError);
                CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
                CHECK_THROWS_AS(vec.getComponentView(0, num_local_patches), RuntimeError);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getComponentView const managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int nz : { 1, 4, 5 }) {
            for (int num_local_patches : { 1, 13 }) {
              Communicator comm(MPI_COMM_WORLD);

              array<int, 3> ns = { nx, ny, nz };
              int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells);
              int patch_stride = component_stride * num_components;

              const Vector<3> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

              const double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (int i = 0; i < num_local_patches; i++) {
                for (int c = 0; c < num_components; c++) {
                  ComponentView<const double, 3> ld = vec.getComponentView(c, i);
                  CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i + c * component_stride);
                  CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * i + (c + 1) * component_stride - 1);
                }
              }
              if (ENABLE_DEBUG) {
                CHECK_THROWS_AS(vec.getComponentView(-1, 0), RuntimeError);
                CHECK_THROWS_AS(vec.getComponentView(num_components, 0), RuntimeError);
                CHECK_THROWS_AS(vec.getComponentView(0, -1), RuntimeError);
                CHECK_THROWS_AS(vec.getComponentView(0, num_local_patches), RuntimeError);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<1> getPatchView managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int num_local_patches : { 1, 13 }) {
          Communicator comm(MPI_COMM_WORLD);

          array<int, 1> ns = { nx };
          int component_stride = (nx + 2 * num_ghost_cells);
          int patch_stride = component_stride * num_components;

          Vector<1> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

          double* view = &vec.getPatchView(0)(-num_ghost_cells, 0);
          for (int i = 0; i < num_local_patches; i++) {
            PatchView<double, 1> ld = vec.getPatchView(i);
            CHECK_EQ(&ld[ld.getGhostStart()], view + (patch_stride * i));
            CHECK_EQ(&ld[ld.getGhostEnd()], view + (patch_stride * (i + 1) - 1));
          }
          if (ENABLE_DEBUG) {
            CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
            CHECK_THROWS_AS(vec.getPatchView(num_local_patches), RuntimeError);
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<1> getPatchView const managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int num_local_patches : { 1, 13 }) {
          Communicator comm(MPI_COMM_WORLD);

          array<int, 1> ns = { nx };
          int component_stride = (nx + 2 * num_ghost_cells);
          int patch_stride = component_stride * num_components;

          const Vector<1> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

          const double* view = &vec.getPatchView(0)(-num_ghost_cells, 0);
          for (int i = 0; i < num_local_patches; i++) {
            PatchView<const double, 1> ld = vec.getPatchView(i);
            CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i);
            CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * (i + 1) - 1);
          }
          if (ENABLE_DEBUG) {
            CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
            CHECK_THROWS_AS(vec.getPatchView(num_local_patches), RuntimeError);
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getPatchView managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 13 }) {
            Communicator comm(MPI_COMM_WORLD);

            array<int, 2> ns = { nx, ny };
            int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
            int patch_stride = component_stride * num_components;

            Vector<2> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

            double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            for (int i = 0; i < num_local_patches; i++) {
              PatchView<double, 2> ld = vec.getPatchView(i);
              CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i);
              CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * (i + 1) - 1);
            }
            if (ENABLE_DEBUG) {
              CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
              CHECK_THROWS_AS(vec.getPatchView(num_local_patches), RuntimeError);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<2> getPatchView const managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int num_local_patches : { 1, 13 }) {
            Communicator comm(MPI_COMM_WORLD);

            array<int, 2> ns = { nx, ny };
            int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells);
            int patch_stride = component_stride * num_components;

            const Vector<2> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

            const double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, 0);
            for (int i = 0; i < num_local_patches; i++) {
              PatchView<const double, 2> ld = vec.getPatchView(i);
              CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i);
              CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * (i + 1) - 1);
            }

            if (ENABLE_DEBUG) {
              CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
              CHECK_THROWS_AS(vec.getPatchView(num_local_patches), RuntimeError);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getPatchView managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int nz : { 1, 4, 5 }) {
            for (int num_local_patches : { 1, 13 }) {
              Communicator comm(MPI_COMM_WORLD);

              array<int, 3> ns = { nx, ny, nz };
              int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells);
              int patch_stride = component_stride * num_components;

              Vector<3> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

              double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (int i = 0; i < num_local_patches; i++) {
                PatchView<double, 3> ld = vec.getPatchView(i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * (i + 1) - 1);
              }
              if (ENABLE_DEBUG) {
                CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
                CHECK_THROWS_AS(vec.getPatchView(num_local_patches), RuntimeError);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getPatchView const managed constructor")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int nz : { 1, 4, 5 }) {
            for (int num_local_patches : { 1, 13 }) {
              Communicator comm(MPI_COMM_WORLD);

              array<int, 3> ns = { nx, ny, nz };
              int component_stride = (nx + 2 * num_ghost_cells) * (ny + 2 * num_ghost_cells) * (nz + 2 * num_ghost_cells);
              int patch_stride = component_stride * num_components;

              const Vector<3> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

              const double* view = &vec.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (int i = 0; i < num_local_patches; i++) {
                PatchView<const double, 3> ld = vec.getPatchView(i);
                CHECK_EQ(&ld[ld.getGhostStart()], view + patch_stride * i);
                CHECK_EQ(&ld[ld.getGhostEnd()], view + patch_stride * (i + 1) - 1);
              }
              if (ENABLE_DEBUG) {
                CHECK_THROWS_AS(vec.getPatchView(-1), RuntimeError);
                CHECK_THROWS_AS(vec.getPatchView(num_local_patches), RuntimeError);
              }
            }
          }
        }
      }
    }
  }
}
