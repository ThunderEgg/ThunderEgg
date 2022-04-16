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

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;
#define MESHES "mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";
TEST_CASE("Vector<3> getMPIComm")
{
  for (int num_components : { 1, 2, 3 }) {
    for (int num_ghost_cells : { 0, 1, 5 }) {
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

              // CHECK(vec.getMPIComm() == comm);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getNumComponents")
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

              CHECK(vec.getNumComponents() == num_components);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getNumLocalPatches")
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

              CHECK(vec.getNumLocalPatches() == num_local_patches);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getNumLocalCells")
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

              CHECK(vec.getNumLocalCells() == nx * ny * nz * num_local_patches);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getComponentViews")
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

              for (int i = 0; i < vec.getNumLocalPatches(); i++) {
                auto lds = vec.getPatchView(i);
                for (int c = 0; c < vec.getNumComponents(); c++) {
                  auto ld = vec.getComponentView(c, i);
                  CHECK(&ld[{ 0, 0, 0 }] == &lds[{ 0, 0, 0, c }]);
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> getComponentViews const")
{
  for (int num_components : { 1, 2, 3 }) {
    for (auto num_ghost_cells : { 0, 1, 5 }) {
      for (int nx : { 1, 4, 5 }) {
        for (int ny : { 1, 4, 5 }) {
          for (int nz : { 1, 4, 5 }) {
            for (int num_local_patches : { 1, 13 }) {
              Communicator comm(MPI_COMM_WORLD);
              array<int, 3> ns = { nx, ny, nz };

              const Vector<3> vec(comm, ns, num_components, num_local_patches, num_ghost_cells);

              INFO("num_ghost_cells:   " << num_ghost_cells);
              INFO("nx:                " << nx);
              INFO("ny:                " << ny);
              INFO("nz:                " << nz);
              INFO("num_local_patches: " << num_local_patches);
              INFO("num_components:    " << num_components);

              for (int i = 0; i < vec.getNumLocalPatches(); i++) {
                auto lds = vec.getPatchView(i);
                for (int c = 0; c < vec.getNumComponents(); c++) {
                  auto ld = vec.getComponentView(c, i);
                  CHECK(&ld[{ 0, 0, 0 }] == &lds[{ 0, 0, 0, c }]);
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> set")
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

              vec.set(28);

              for (int i = 0; i < vec.getNumLocalPatches(); i++) {
                for (int c = 0; c < vec.getNumComponents(); c++) {
                  auto ld = vec.getComponentView(c, i);
                  Loop::Nested<3>(ld.getGhostStart(), ld.getGhostEnd(), [&](std::array<int, 3>& coord) {
                    if (isGhost(coord, ns, num_ghost_cells)) {
                      CHECK(ld[coord] == 0);
                    } else {
                      CHECK(ld[coord] == 28);
                    }
                  });
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> setWithGhost")
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

              vec.setWithGhost(28);

              for (int i = 0; i < vec.getNumLocalPatches(); i++) {
                for (int c = 0; c < vec.getNumComponents(); c++) {
                  auto ld = vec.getComponentView(c, i);
                  Loop::Nested<3>(ld.getGhostStart(), ld.getGhostEnd(), [&](std::array<int, 3>& coord) { CHECK(ld[coord] == 28); });
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> scale")
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

              vec.setWithGhost(28);
              vec.scale(0.25);

              for (int i = 0; i < vec.getNumLocalPatches(); i++) {
                for (int c = 0; c < vec.getNumComponents(); c++) {
                  auto ld = vec.getComponentView(c, i);
                  Loop::Nested<3>(ld.getGhostStart(), ld.getGhostEnd(), [&](std::array<int, 3>& coord) {
                    if (isGhost(coord, ns, num_ghost_cells)) {
                      CHECK(ld[coord] == 28);
                    } else {
                      CHECK(ld[coord] == Catch::Approx(7));
                    }
                  });
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> shift")
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

              vec.setWithGhost(1);
              vec.shift(28);

              for (int i = 0; i < vec.getNumLocalPatches(); i++) {
                for (int c = 0; c < vec.getNumComponents(); c++) {
                  auto ld = vec.getComponentView(c, i);
                  Loop::Nested<3>(ld.getGhostStart(), ld.getGhostEnd(), [&](std::array<int, 3>& coord) {
                    if (isGhost(coord, ns, num_ghost_cells)) {
                      CHECK(ld[coord] == 1);
                    } else {
                      CHECK(ld[coord] == Catch::Approx(29));
                    }
                  });
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> copy")
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

              b.copy(a);

              for (int i = 0; i < a.getNumLocalPatches(); i++) {
                for (int c = 0; c < a.getNumComponents(); c++) {
                  auto a_ld = a.getComponentView(c, i);
                  auto b_ld = b.getComponentView(c, i);
                  Loop::Nested<3>(a_ld.getGhostStart(), a_ld.getGhostEnd(), [&](std::array<int, 3>& coord) {
                    if (isGhost(coord, ns, num_ghost_cells)) {
                      CHECK(b_ld[coord] == 0);
                    } else {
                      CHECK(b_ld[coord] == a_ld[coord]);
                    }
                  });
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> copyWithGhost")
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

              b.copyWithGhost(a);

              for (int i = 0; i < a.getNumLocalPatches(); i++) {
                for (int c = 0; c < a.getNumComponents(); c++) {
                  auto a_ld = a.getComponentView(c, i);
                  auto b_ld = b.getComponentView(c, i);
                  Loop::Nested<3>(a_ld.getGhostStart(), a_ld.getGhostEnd(), [&](std::array<int, 3>& coord) { CHECK(b_ld[coord] == a_ld[coord]); });
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> add")
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
              Vector<3> b_copy(comm, ns, num_components, num_local_patches, num_ghost_cells);
              Vector<3> expected(comm, ns, num_components, num_local_patches, num_ghost_cells);

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
              double* b_copy_data = &b_copy.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                b_data[i] = (x - 0.5) * (x - 0.5);
                b_copy_data[i] = (x - 0.5) * (x - 0.5);
              }

              double* expected_data = &expected.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                expected_data[i] = 10 - (x - 0.75) * (x - 0.75) + (x - 0.5) * (x - 0.5);
              }
              b.add(a);

              for (int i = 0; i < a.getNumLocalPatches(); i++) {
                for (int c = 0; c < a.getNumComponents(); c++) {
                  auto expected_ld = expected.getComponentView(c, i);
                  auto b_ld = b.getComponentView(c, i);
                  auto b_copy_ld = b_copy.getComponentView(c, i);
                  Loop::Nested<3>(b_ld.getGhostStart(), b_ld.getGhostEnd(), [&](std::array<int, 3>& coord) {
                    if (isGhost(coord, ns, num_ghost_cells)) {
                      CHECK(b_ld[coord] == b_copy_ld[coord]);
                    } else {
                      CHECK(b_ld[coord] == Catch::Approx(expected_ld[coord]));
                    }
                  });
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> addScaled")
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
              Vector<3> b_copy(comm, ns, num_components, num_local_patches, num_ghost_cells);
              Vector<3> expected(comm, ns, num_components, num_local_patches, num_ghost_cells);

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
              double* b_copy_data = &b_copy.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                b_data[i] = (x - 0.5) * (x - 0.5);
                b_copy_data[i] = (x - 0.5) * (x - 0.5);
              }

              double* expected_data = &expected.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                expected_data[i] = 0.7 * (10 - (x - 0.75) * (x - 0.75)) + (x - 0.5) * (x - 0.5);
              }
              b.addScaled(0.7, a);

              for (int i = 0; i < a.getNumLocalPatches(); i++) {
                for (int c = 0; c < a.getNumComponents(); c++) {
                  auto expected_ld = expected.getComponentView(c, i);
                  auto b_ld = b.getComponentView(c, i);
                  auto b_copy_ld = b_copy.getComponentView(c, i);
                  Loop::Nested<3>(b_ld.getGhostStart(), b_ld.getGhostEnd(), [&](std::array<int, 3>& coord) {
                    if (isGhost(coord, ns, num_ghost_cells)) {
                      CHECK(b_ld[coord] == b_copy_ld[coord]);
                    } else {
                      CHECK(b_ld[coord] == Catch::Approx(expected_ld[coord]));
                    }
                  });
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> scaleThenAdd")
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
              Vector<3> b_copy(comm, ns, num_components, num_local_patches, num_ghost_cells);
              Vector<3> expected(comm, ns, num_components, num_local_patches, num_ghost_cells);

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
              double* b_copy_data = &b_copy.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                b_data[i] = (x - 0.5) * (x - 0.5);
                b_copy_data[i] = (x - 0.5) * (x - 0.5);
              }

              double* expected_data = &expected.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                expected_data[i] = (10 - (x - 0.75) * (x - 0.75)) + 0.7 * (x - 0.5) * (x - 0.5);
              }
              b.scaleThenAdd(0.7, a);

              for (int i = 0; i < a.getNumLocalPatches(); i++) {
                for (int c = 0; c < a.getNumComponents(); c++) {
                  auto expected_ld = expected.getComponentView(c, i);
                  auto b_ld = b.getComponentView(c, i);
                  auto b_copy_ld = b_copy.getComponentView(c, i);
                  Loop::Nested<3>(b_ld.getGhostStart(), b_ld.getGhostEnd(), [&](std::array<int, 3>& coord) {
                    if (isGhost(coord, ns, num_ghost_cells)) {
                      CHECK(b_ld[coord] == b_copy_ld[coord]);
                    } else {
                      CHECK(b_ld[coord] == Catch::Approx(expected_ld[coord]));
                    }
                  });
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> scaleThenAddScaled")
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
              Vector<3> b_copy(comm, ns, num_components, num_local_patches, num_ghost_cells);
              Vector<3> expected(comm, ns, num_components, num_local_patches, num_ghost_cells);

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
              double* b_copy_data = &b_copy.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                b_data[i] = (x - 0.5) * (x - 0.5);
                b_copy_data[i] = (x - 0.5) * (x - 0.5);
              }

              double* expected_data = &expected.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                expected_data[i] = -2 * (10 - (x - 0.75) * (x - 0.75)) + 0.7 * (x - 0.5) * (x - 0.5);
              }
              b.scaleThenAddScaled(0.7, -2, a);

              for (int i = 0; i < a.getNumLocalPatches(); i++) {
                for (int c = 0; c < a.getNumComponents(); c++) {
                  auto expected_ld = expected.getComponentView(c, i);
                  auto b_ld = b.getComponentView(c, i);
                  auto b_copy_ld = b_copy.getComponentView(c, i);
                  Loop::Nested<3>(b_ld.getGhostStart(), b_ld.getGhostEnd(), [&](std::array<int, 3>& coord) {
                    if (isGhost(coord, ns, num_ghost_cells)) {
                      CHECK(b_ld[coord] == b_copy_ld[coord]);
                    } else {
                      CHECK(b_ld[coord] == Catch::Approx(expected_ld[coord]));
                    }
                  });
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("Vector<3> scaleThenAddScaled two vectors")
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
              Vector<3> c(comm, ns, num_components, num_local_patches, num_ghost_cells);
              Vector<3> b_copy(comm, ns, num_components, num_local_patches, num_ghost_cells);
              Vector<3> expected(comm, ns, num_components, num_local_patches, num_ghost_cells);

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
              double* b_copy_data = &b_copy.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                b_data[i] = (x - 0.5) * (x - 0.5);
                b_copy_data[i] = (x - 0.5) * (x - 0.5);
              }

              double* c_data = &c.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                c_data[i] = 1 + (x - 0.25) * (x - 0.25);
              }

              double* expected_data = &expected.getPatchView(0)(-num_ghost_cells, -num_ghost_cells, -num_ghost_cells, 0);
              for (size_t i = 0; i < size; i++) {
                double x = (i + 0.5) / size;
                expected_data[i] = -2 * (10 - (x - 0.75) * (x - 0.75)) + 0.7 * (x - 0.5) * (x - 0.5) + 9 * (1 + (x - 0.25) * (x - 0.25));
              }
              b.scaleThenAddScaled(0.7, -2, a, 9, c);

              for (int i = 0; i < a.getNumLocalPatches(); i++) {
                for (int c = 0; c < a.getNumComponents(); c++) {
                  auto expected_ld = expected.getComponentView(c, i);
                  auto b_ld = b.getComponentView(c, i);
                  auto b_copy_ld = b_copy.getComponentView(c, i);
                  Loop::Nested<3>(b_ld.getGhostStart(), b_ld.getGhostEnd(), [&](std::array<int, 3>& coord) {
                    if (isGhost(coord, ns, num_ghost_cells)) {
                      CHECK(b_ld[coord] == b_copy_ld[coord]);
                    } else {
                      CHECK(b_ld[coord] == Catch::Approx(expected_ld[coord]));
                    }
                  });
                }
              }
            }
          }
        }
      }
    }
  }
}
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
              expected_norm = sqrt(expected_norm);

              CHECK(vec.twoNorm() == Catch::Approx(expected_norm));
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

              CHECK(vec.infNorm() == expected_norm);
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

              CHECK(a.dot(b) == Catch::Approx(expected_value));
            }
          }
        }
      }
    }
  }
}
