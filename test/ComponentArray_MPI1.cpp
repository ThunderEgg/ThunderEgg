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
#include <ThunderEgg/ComponentArray.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
using namespace std;
using namespace ThunderEgg;
TEST_CASE("ComponentArray getStart", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        CHECK(ca.getStart()[0] == 0);
        CHECK(ca.getStart()[1] == 0);
      }
    }
  }
}
TEST_CASE("ComponentArray getEnd", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        CHECK(ca.getEnd()[0] == nx - 1);
        CHECK(ca.getEnd()[1] == ny - 1);
      }
    }
  }
}
TEST_CASE("ComponentArray getGhostStart", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        CHECK(ca.getGhostStart()[0] == -num_ghost);
        CHECK(ca.getGhostStart()[1] == -num_ghost);
      }
    }
  }
}
TEST_CASE("ComponentArray getGhostEnd", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        CHECK(ca.getGhostEnd()[0] == nx - 1 + num_ghost);
        CHECK(ca.getGhostEnd()[1] == ny - 1 + num_ghost);
      }
    }
  }
}
TEST_CASE("ComponentArray getStrides", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        CHECK(ca.getStrides()[0] == 1);
        CHECK(ca.getStrides()[1] == nx + 2 * num_ghost);
      }
    }
  }
}
TEST_CASE("ComponentArray squarebracket operator", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        double* start = &ca[{ 0, 0 }];
        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            CHECK(&ca[{ xi, yi }] == start + xi + yi * (nx + 2 * num_ghost));
          }
        }
        CHECK(ca.getStrides()[0] == 1);
        CHECK(ca.getStrides()[1] == nx + 2 * num_ghost);
      }
    }
  }
}
TEST_CASE("ComponentArray squarebracket operator const", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        const ComponentArray<2> ca({ nx, ny }, num_ghost);

        const double* start = &ca[{ 0, 0 }];
        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            CHECK(&ca[{ xi, yi }] == start + xi + yi * (nx + 2 * num_ghost));
          }
        }
        CHECK(ca.getStrides()[0] == 1);
        CHECK(ca.getStrides()[1] == nx + 2 * num_ghost);
      }
    }
  }
}
TEST_CASE("ComponentArray default is zero", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            CHECK(ca[{ xi, yi }] == 0.0);
          }
        }
        CHECK(ca.getStrides()[0] == 1);
        CHECK(ca.getStrides()[1] == nx + 2 * num_ghost);
      }
    }
  }
}
TEST_CASE("ComponentArray<2> getSliceOn<1>", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
          INFO("xi: " << xi);
          View<double, 1> slice = ca.getSliceOn(Side<2>::west(), { xi });
          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            CHECK(&ca[{ xi, yi }] == &slice[{ yi }]);
          }
        }

        for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
          INFO("xi: " << xi);
          View<double, 1> slice = ca.getSliceOn(Side<2>::east(), { xi });
          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            CHECK(&ca[{ nx - 1 - xi, yi }] == &slice[{ yi }]);
          }
        }

        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          INFO("yi: " << yi);
          View<double, 1> slice = ca.getSliceOn(Side<2>::south(), { yi });
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            CHECK(&ca[{ xi, yi }] == &slice[{ xi }]);
          }
        }

        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          INFO("yi: " << yi);
          View<double, 1> slice = ca.getSliceOn(Side<2>::north(), { yi });
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            CHECK(&ca[{ xi, ny - 1 - yi }] == &slice[{ xi }]);
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<2> getSliceOn<1> const", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        const ComponentArray<2> ca({ nx, ny }, num_ghost);

        for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
          INFO("xi: " << xi);
          View<const double, 1> slice = ca.getSliceOn(Side<2>::west(), { xi });
          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            CHECK(&ca[{ xi, yi }] == &slice[{ yi }]);
          }
        }

        for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
          INFO("xi: " << xi);
          View<const double, 1> slice = ca.getSliceOn(Side<2>::east(), { xi });
          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            CHECK(&ca[{ nx - 1 - xi, yi }] == &slice[{ yi }]);
          }
        }

        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          INFO("yi: " << yi);
          View<const double, 1> slice = ca.getSliceOn(Side<2>::south(), { yi });
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            CHECK(&ca[{ xi, yi }] == &slice[{ xi }]);
          }
        }

        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          INFO("yi: " << yi);
          View<const double, 1> slice = ca.getSliceOn(Side<2>::north(), { yi });
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            CHECK(&ca[{ xi, ny - 1 - yi }] == &slice[{ xi }]);
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<3> getSliceOn<2>", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          ComponentArray<3> ca({ nx, ny, nz }, num_ghost);

          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            View<double, 2> slice = ca.getSliceOn(Side<3>::west(), { xi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ xi, yi, zi }] == &slice[{ yi, zi }]);
              }
            }
          }

          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            View<double, 2> slice = ca.getSliceOn(Side<3>::east(), { xi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ nx - 1 - xi, yi, zi }] == &slice[{ yi, zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            View<double, 2> slice = ca.getSliceOn(Side<3>::south(), { yi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, zi }] == &slice[{ xi, zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            View<double, 2> slice = ca.getSliceOn(Side<3>::north(), { yi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, ny - 1 - yi, zi }] == &slice[{ xi, zi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            View<double, 2> slice = ca.getSliceOn(Side<3>::bottom(), { zi });
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, zi }] == &slice[{ xi, yi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            View<double, 2> slice = ca.getSliceOn(Side<3>::top(), { zi });
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, nz - 1 - zi }] == &slice[{ xi, yi }]);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<3> getSliceOn<2> const", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          const ComponentArray<3> ca({ nx, ny, nz }, num_ghost);

          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            View<const double, 2> slice = ca.getSliceOn(Side<3>::west(), { xi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ xi, yi, zi }] == &slice[{ yi, zi }]);
              }
            }
          }

          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            View<const double, 2> slice = ca.getSliceOn(Side<3>::east(), { xi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ nx - 1 - xi, yi, zi }] == &slice[{ yi, zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            View<const double, 2> slice = ca.getSliceOn(Side<3>::south(), { yi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, zi }] == &slice[{ xi, zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            View<const double, 2> slice = ca.getSliceOn(Side<3>::north(), { yi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, ny - 1 - yi, zi }] == &slice[{ xi, zi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            View<const double, 2> slice = ca.getSliceOn(Side<3>::bottom(), { zi });
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, zi }] == &slice[{ xi, yi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            View<const double, 2> slice = ca.getSliceOn(Side<3>::top(), { zi });
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, nz - 1 - zi }] == &slice[{ xi, yi }]);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<3> getSliceOn<1>", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          ComponentArray<3> ca({ nx, ny, nz }, num_ghost);

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              View<double, 1> slice = ca.getSliceOn(Edge::bs(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, zi }] == &slice[{ xi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              View<double, 1> slice = ca.getSliceOn(Edge::tn(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, ny - 1 - yi, nz - 1 - zi }] == &slice[{ xi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              View<double, 1> slice = ca.getSliceOn(Edge::bn(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, ny - 1 - yi, zi }] == &slice[{ xi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              View<double, 1> slice = ca.getSliceOn(Edge::ts(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, nz - 1 - zi }] == &slice[{ xi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 1> slice = ca.getSliceOn(Edge::bw(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ xi, yi, zi }] == &slice[{ yi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 1> slice = ca.getSliceOn(Edge::te(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ nx - 1 - xi, yi, nz - 1 - zi }] == &slice[{ yi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 1> slice = ca.getSliceOn(Edge::be(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ nx - 1 - xi, yi, zi }] == &slice[{ yi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 1> slice = ca.getSliceOn(Edge::tw(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ xi, yi, nz - 1 - zi }] == &slice[{ yi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 1> slice = ca.getSliceOn(Edge::sw(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                INFO("zi: " << zi);
                CHECK(&ca[{ xi, yi, zi }] == &slice[{ zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 1> slice = ca.getSliceOn(Edge::ne(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                INFO("zi: " << zi);
                CHECK(&ca[{ nx - 1 - xi, ny - 1 - yi, zi }] == &slice[{ zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 1> slice = ca.getSliceOn(Edge::se(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                INFO("zi: " << zi);
                CHECK(&ca[{ nx - 1 - xi, yi, zi }] == &slice[{ zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 1> slice = ca.getSliceOn(Edge::nw(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                INFO("zi: " << zi);
                CHECK(&ca[{ xi, ny - 1 - yi, zi }] == &slice[{ zi }]);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<3> getSliceOn<1> const", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          const ComponentArray<3> ca({ nx, ny, nz }, num_ghost);

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              View<const double, 1> slice = ca.getSliceOn(Edge::bs(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, zi }] == &slice[{ xi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              View<const double, 1> slice = ca.getSliceOn(Edge::tn(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, ny - 1 - yi, nz - 1 - zi }] == &slice[{ xi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              View<const double, 1> slice = ca.getSliceOn(Edge::bn(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, ny - 1 - yi, zi }] == &slice[{ xi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              View<const double, 1> slice = ca.getSliceOn(Edge::ts(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, nz - 1 - zi }] == &slice[{ xi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<const double, 1> slice = ca.getSliceOn(Edge::bw(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ xi, yi, zi }] == &slice[{ yi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<const double, 1> slice = ca.getSliceOn(Edge::te(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ nx - 1 - xi, yi, nz - 1 - zi }] == &slice[{ yi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<const double, 1> slice = ca.getSliceOn(Edge::be(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ nx - 1 - xi, yi, zi }] == &slice[{ yi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<const double, 1> slice = ca.getSliceOn(Edge::tw(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ xi, yi, nz - 1 - zi }] == &slice[{ yi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<const double, 1> slice = ca.getSliceOn(Edge::sw(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                INFO("zi: " << zi);
                CHECK(&ca[{ xi, yi, zi }] == &slice[{ zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<const double, 1> slice = ca.getSliceOn(Edge::ne(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                INFO("zi: " << zi);
                CHECK(&ca[{ nx - 1 - xi, ny - 1 - yi, zi }] == &slice[{ zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<const double, 1> slice = ca.getSliceOn(Edge::se(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                INFO("zi: " << zi);
                CHECK(&ca[{ nx - 1 - xi, yi, zi }] == &slice[{ zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<const double, 1> slice = ca.getSliceOn(Edge::nw(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                INFO("zi: " << zi);
                CHECK(&ca[{ xi, ny - 1 - yi, zi }] == &slice[{ zi }]);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<3> getSliceOn<0>", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          ComponentArray<3> ca({ nx, ny, nz }, num_ghost);

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, zi }] == &(ca.getSliceOn(Corner<3>::bsw(), { xi, yi, zi }))[{}]);
                CHECK(&ca[{ nx - 1 - xi, yi, zi }] == &(ca.getSliceOn(Corner<3>::bse(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ xi, ny - 1 - yi, zi }] == &(ca.getSliceOn(Corner<3>::bnw(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ nx - 1 - xi, ny - 1 - yi, zi }] == &(ca.getSliceOn(Corner<3>::bne(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ xi, yi, nz - 1 - zi }] == &(ca.getSliceOn(Corner<3>::tsw(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ nx - 1 - xi, yi, nz - 1 - zi }] == &(ca.getSliceOn(Corner<3>::tse(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ xi, ny - 1 - yi, nz - 1 - zi }] == &(ca.getSliceOn(Corner<3>::tnw(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ nx - 1 - xi, ny - 1 - yi, nz - 1 - zi }] == &(ca.getSliceOn(Corner<3>::tne(), { xi, yi, zi })[{}]));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<3> getSliceOn<0> const", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          const ComponentArray<3> ca({ nx, ny, nz }, num_ghost);

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            INFO("zi: " << zi);
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, zi }] == &(ca.getSliceOn(Corner<3>::bsw(), { xi, yi, zi }))[{}]);
                CHECK(&ca[{ nx - 1 - xi, yi, zi }] == &(ca.getSliceOn(Corner<3>::bse(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ xi, ny - 1 - yi, zi }] == &(ca.getSliceOn(Corner<3>::bnw(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ nx - 1 - xi, ny - 1 - yi, zi }] == &(ca.getSliceOn(Corner<3>::bne(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ xi, yi, nz - 1 - zi }] == &(ca.getSliceOn(Corner<3>::tsw(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ nx - 1 - xi, yi, nz - 1 - zi }] == &(ca.getSliceOn(Corner<3>::tse(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ xi, ny - 1 - yi, nz - 1 - zi }] == &(ca.getSliceOn(Corner<3>::tnw(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ nx - 1 - xi, ny - 1 - yi, nz - 1 - zi }] == &(ca.getSliceOn(Corner<3>::tne(), { xi, yi, zi })[{}]));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<2> getSliceOn<0>", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          INFO("yi: " << yi);
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            CHECK(&ca[{ xi, yi }] == &(ca.getSliceOn(Corner<2>::sw(), { xi, yi })[{}]));
            CHECK(&ca[{ nx - 1 - xi, yi }] == &(ca.getSliceOn(Corner<2>::se(), { xi, yi })[{}]));
            CHECK(&ca[{ xi, ny - 1 - yi }] == &(ca.getSliceOn(Corner<2>::nw(), { xi, yi })[{}]));
            CHECK(&ca[{ nx - 1 - xi, ny - 1 - yi }] == &(ca.getSliceOn(Corner<2>::ne(), { xi, yi })[{}]));
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<2> getSliceOn<0> const", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        const ComponentArray<2> ca({ nx, ny }, num_ghost);

        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          INFO("yi: " << yi);
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            CHECK(&ca[{ xi, yi }] == &(ca.getSliceOn(Corner<2>::sw(), { xi, yi })[{}]));
            CHECK(&ca[{ nx - 1 - xi, yi }] == &(ca.getSliceOn(Corner<2>::se(), { xi, yi })[{}]));
            CHECK(&ca[{ xi, ny - 1 - yi }] == &(ca.getSliceOn(Corner<2>::nw(), { xi, yi })[{}]));
            CHECK(&ca[{ nx - 1 - xi, ny - 1 - yi }] == &(ca.getSliceOn(Corner<2>::ne(), { xi, yi })[{}]));
          }
        }
      }
    }
  }
}

TEST_CASE("ComponentArray<2> getGhostSliceOn<1>", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        for (unsigned char xi = 0; xi < num_ghost; xi++) {
          INFO("xi: " << xi);
          View<double, 1> slice_w = ca.getGhostSliceOn(Side<2>::west(), { xi });
          View<double, 1> slice_e = ca.getGhostSliceOn(Side<2>::east(), { xi });
          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            CHECK(&ca[{ -1 - xi, yi }] == &slice_w[{ yi }]);
            CHECK(&ca[{ nx + xi, yi }] == &slice_e[{ yi }]);
          }
        }

        for (unsigned char yi = 0; yi < num_ghost; yi++) {
          INFO("yi: " << yi);
          View<double, 1> slice_s = ca.getGhostSliceOn(Side<2>::south(), { yi });
          View<double, 1> slice_n = ca.getGhostSliceOn(Side<2>::north(), { yi });
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            CHECK(&ca[{ xi, -1 - yi }] == &slice_s[{ xi }]);
            CHECK(&ca[{ xi, ny + yi }] == &slice_n[{ xi }]);
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<3> getGhostSliceOn<2>", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          const ComponentArray<3> ca({ nx, ny, nz }, num_ghost);

          for (unsigned char xi = 0; xi < num_ghost; xi++) {
            INFO("xi: " << xi);
            View<double, 2> slice_w = ca.getGhostSliceOn(Side<3>::west(), { xi });
            View<double, 2> slice_e = ca.getGhostSliceOn(Side<3>::east(), { xi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ -1 - xi, yi, zi }] == &slice_w[{ yi, zi }]);
                CHECK(&ca[{ nx + xi, yi, zi }] == &slice_e[{ yi, zi }]);
              }
            }
          }

          for (unsigned char yi = 0; yi < num_ghost; yi++) {
            INFO("yi: " << yi);
            View<double, 2> slice_s = ca.getGhostSliceOn(Side<3>::south(), { yi });
            View<double, 2> slice_n = ca.getGhostSliceOn(Side<3>::north(), { yi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, -1 - yi, zi }] == &slice_s[{ xi, zi }]);
                CHECK(&ca[{ xi, ny + yi, zi }] == &slice_n[{ xi, zi }]);
              }
            }
          }

          for (unsigned char zi = 0; zi < num_ghost; zi++) {
            INFO("zi: " << zi);
            View<double, 2> slice_b = ca.getGhostSliceOn(Side<3>::bottom(), { zi });
            View<double, 2> slice_t = ca.getGhostSliceOn(Side<3>::top(), { zi });
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, yi, -1 - zi }] == &slice_b[{ xi, yi }]);
                CHECK(&ca[{ xi, yi, nz + zi }] == &slice_t[{ xi, yi }]);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<3> getGhostSliceOn<1>", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          ComponentArray<3> ca({ nx, ny, nz }, num_ghost);

          for (unsigned char zi = 0; zi < num_ghost; zi++) {
            INFO("zi: " << zi);
            for (unsigned char yi = 0; yi < num_ghost; yi++) {
              INFO("yi: " << yi);
              View<double, 1> slice_bs = ca.getGhostSliceOn(Edge::bs(), { yi, zi });
              View<double, 1> slice_tn = ca.getGhostSliceOn(Edge::tn(), { yi, zi });
              View<double, 1> slice_bn = ca.getGhostSliceOn(Edge::bn(), { yi, zi });
              View<double, 1> slice_ts = ca.getGhostSliceOn(Edge::ts(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ xi, -1 - yi, -1 - zi }] == &slice_bs[{ xi }]);
                CHECK(&ca[{ xi, ny + yi, nz + zi }] == &slice_tn[{ xi }]);
                CHECK(&ca[{ xi, ny + yi, -1 - zi }] == &slice_bn[{ xi }]);
                CHECK(&ca[{ xi, -1 - yi, nz + zi }] == &slice_ts[{ xi }]);
              }
            }
          }

          for (unsigned char zi = 0; zi < num_ghost; zi++) {
            INFO("zi: " << zi);
            for (unsigned char xi = 0; xi < num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 1> slice_bw = ca.getGhostSliceOn(Edge::bw(), { xi, zi });
              View<double, 1> slice_te = ca.getGhostSliceOn(Edge::te(), { xi, zi });
              View<double, 1> slice_be = ca.getGhostSliceOn(Edge::be(), { xi, zi });
              View<double, 1> slice_tw = ca.getGhostSliceOn(Edge::tw(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&ca[{ -1 - xi, yi, -1 - zi }] == &slice_bw[{ yi }]);
                CHECK(&ca[{ nx + xi, yi, nz + zi }] == &slice_te[{ yi }]);
                CHECK(&ca[{ nx + xi, yi, -1 - zi }] == &slice_be[{ yi }]);
                CHECK(&ca[{ -1 - xi, yi, nz + zi }] == &slice_tw[{ yi }]);
              }
            }
          }

          for (unsigned char yi = 0; yi < num_ghost; yi++) {
            INFO("yi: " << yi);
            for (unsigned char xi = 0; xi < num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 1> slice_sw = ca.getGhostSliceOn(Edge::sw(), { xi, yi });
              View<double, 1> slice_ne = ca.getGhostSliceOn(Edge::ne(), { xi, yi });
              View<double, 1> slice_se = ca.getGhostSliceOn(Edge::se(), { xi, yi });
              View<double, 1> slice_nw = ca.getGhostSliceOn(Edge::nw(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                INFO("zi: " << zi);
                CHECK(&ca[{ -1 - xi, -1 - yi, zi }] == &slice_sw[{ zi }]);
                CHECK(&ca[{ nx + xi, ny + yi, zi }] == &slice_ne[{ zi }]);
                CHECK(&ca[{ nx + xi, -1 - yi, zi }] == &slice_se[{ zi }]);
                CHECK(&ca[{ -1 - xi, ny + yi, zi }] == &slice_nw[{ zi }]);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<3> getGhostSliceOn<0>", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          ComponentArray<3> ca({ nx, ny, nz }, num_ghost);

          for (unsigned char zi = 0; zi < num_ghost; zi++) {
            INFO("zi: " << zi);
            for (unsigned char yi = 0; yi < num_ghost; yi++) {
              INFO("yi: " << yi);
              for (unsigned char xi = 0; xi < num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&ca[{ -1 - xi, -1 - yi, -1 - zi }] == &(ca.getGhostSliceOn(Corner<3>::bsw(), { xi, yi, zi }))[{}]);
                CHECK(&ca[{ nx + xi, -1 - yi, -1 - zi }] == &(ca.getGhostSliceOn(Corner<3>::bse(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ -1 - xi, ny + yi, -1 - zi }] == &(ca.getGhostSliceOn(Corner<3>::bnw(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ nx + xi, ny + yi, -1 - zi }] == &(ca.getGhostSliceOn(Corner<3>::bne(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ -1 - xi, -1 - yi, nz + zi }] == &(ca.getGhostSliceOn(Corner<3>::tsw(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ nx + xi, -1 - yi, nz + zi }] == &(ca.getGhostSliceOn(Corner<3>::tse(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ -1 - xi, ny + yi, nz + zi }] == &(ca.getGhostSliceOn(Corner<3>::tnw(), { xi, yi, zi })[{}]));
                CHECK(&ca[{ nx + xi, ny + yi, nz + zi }] == &(ca.getGhostSliceOn(Corner<3>::tne(), { xi, yi, zi })[{}]));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray<2> getGhostSliceOn<0>", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        for (unsigned char yi = 0; yi < num_ghost; yi++) {
          INFO("yi: " << yi);
          for (unsigned char xi = 0; xi < num_ghost; xi++) {
            INFO("xi: " << xi);
            CHECK(&ca[{ -xi - 1, -yi - 1 }] == &(ca.getGhostSliceOn(Corner<2>::sw(), { xi, yi })[{}]));
            CHECK(&ca[{ nx + xi, -yi - 1 }] == &(ca.getGhostSliceOn(Corner<2>::se(), { xi, yi })[{}]));
            CHECK(&ca[{ -xi - 1, ny + yi }] == &(ca.getGhostSliceOn(Corner<2>::nw(), { xi, yi })[{}]));
            CHECK(&ca[{ nx + xi, ny + yi }] == &(ca.getGhostSliceOn(Corner<2>::ne(), { xi, yi })[{}]));
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentArray copy constructor", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        for (int iy = -num_ghost; iy < ny + num_ghost; iy++) {
          for (int ix = -num_ghost; ix < nx + num_ghost; ix++) {
            ca(ix, iy) = ix + iy;
          }
        }

        ComponentArray<2> ca_copy(ca);

        for (int iy = -num_ghost; iy < ny + num_ghost; iy++) {
          for (int ix = -num_ghost; ix < nx + num_ghost; ix++) {
            CHECK(ca(ix, iy) == ca_copy(ix, iy));
          }
        }

        CHECK(ca.getStrides() == ca_copy.getStrides());
        CHECK(ca.getGhostStart() == ca_copy.getGhostStart());
        CHECK(ca.getStart() == ca_copy.getStart());
        CHECK(ca.getEnd() == ca_copy.getEnd());
        CHECK(ca.getGhostEnd() == ca_copy.getGhostEnd());
        CHECK(&ca(0, 0) != &ca_copy(0, 0));
      }
    }
  }
}
TEST_CASE("ComponentArray copy assignment", "[ComponentArray]")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        ComponentArray<2> ca({ nx, ny }, num_ghost);

        for (int iy = -num_ghost; iy < ny + num_ghost; iy++) {
          for (int ix = -num_ghost; ix < nx + num_ghost; ix++) {
            ca(ix, iy) = ix + iy;
          }
        }

        ComponentArray<2> ca_copy({ nx, ny }, num_ghost);
        ca_copy = ca;

        for (int iy = -num_ghost; iy < ny + num_ghost; iy++) {
          for (int ix = -num_ghost; ix < nx + num_ghost; ix++) {
            CHECK(ca(ix, iy) == ca_copy(ix, iy));
          }
        }

        CHECK(ca.getStrides() == ca_copy.getStrides());
        CHECK(ca.getGhostStart() == ca_copy.getGhostStart());
        CHECK(ca.getStart() == ca_copy.getStart());
        CHECK(ca.getEnd() == ca_copy.getEnd());
        CHECK(ca.getGhostEnd() == ca_copy.getGhostEnd());
        CHECK(&ca(0, 0) != &ca_copy(0, 0));
      }
    }
  }
}
