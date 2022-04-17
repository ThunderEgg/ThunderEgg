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
#include <ThunderEgg/PatchArray.h>
#include <doctest.h>

using namespace std;
using namespace ThunderEgg;
TEST_CASE("PatchArray default constructor getStart")
{
  PatchArray<2> pa;

  CHECK(pa.getStart()[0] == 0);
  CHECK(pa.getStart()[1] == 0);
  CHECK(pa.getStart()[2] == 0);
}
TEST_CASE("PatchArray default constructor getEnd")
{
  PatchArray<2> pa;

  CHECK(pa.getEnd()[0] == -1);
  CHECK(pa.getEnd()[1] == -1);
  CHECK(pa.getEnd()[2] == -1);
}
TEST_CASE("PatchArray default constructor getGhostStart")
{
  PatchArray<2> pa;

  CHECK(pa.getGhostStart()[0] == 0);
  CHECK(pa.getGhostStart()[1] == 0);
  CHECK(pa.getGhostStart()[2] == 0);
}
TEST_CASE("PatchArray default constructor getGhostEnd")
{
  PatchArray<2> pa;

  CHECK(pa.getGhostEnd()[0] == -1);
  CHECK(pa.getGhostEnd()[1] == -1);
  CHECK(pa.getGhostEnd()[2] == -1);
}
TEST_CASE("PatchArray getStart")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          CHECK(pa.getStart()[0] == 0);
          CHECK(pa.getStart()[1] == 0);
          CHECK(pa.getStart()[2] == 0);
        }
      }
    }
  }
}
TEST_CASE("PatchArray getEnd")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          CHECK(pa.getEnd()[0] == nx - 1);
          CHECK(pa.getEnd()[1] == ny - 1);
          CHECK(pa.getEnd()[2] == num_components - 1);
        }
      }
    }
  }
}
TEST_CASE("PatchArray getGhostStart")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          CHECK(pa.getGhostStart()[0] == -num_ghost);
          CHECK(pa.getGhostStart()[1] == -num_ghost);
          CHECK(pa.getGhostStart()[2] == 0);
        }
      }
    }
  }
}
TEST_CASE("PatchArray getGhostEnd")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          CHECK(pa.getGhostEnd()[0] == nx - 1 + num_ghost);
          CHECK(pa.getGhostEnd()[1] == ny - 1 + num_ghost);
          CHECK(pa.getGhostEnd()[2] == num_components - 1);
        }
      }
    }
  }
}
TEST_CASE("PatchArray getStrides")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          CHECK(pa.getStrides()[0] == 1);
          CHECK(pa.getStrides()[1] == nx + 2 * num_ghost);
          CHECK(pa.getStrides()[2] == (nx + 2 * num_ghost) * (ny + 2 * num_ghost));
        }
      }
    }
  }
}
TEST_CASE("PatchArray squarebracket operator")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          double* start = &pa(0, 0, 0);
          for (unsigned char ci = 0; ci < num_components; ci++) {
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK(&pa[{ xi, yi, ci }] == start + xi + yi * (nx + 2 * num_ghost) + ci * (nx + 2 * num_ghost) * (ny + 2 * num_ghost));
              }
            }
          }
        }
      }
    }
  }
}

TEST_CASE("PatchArray squarebracket operator const")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          const PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          const double* start = &pa(0, 0, 0);
          for (unsigned char ci = 0; ci < num_components; ci++) {
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK(&pa[{ xi, yi, ci }] == start + xi + yi * (nx + 2 * num_ghost) + ci * (nx + 2 * num_ghost) * (ny + 2 * num_ghost));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray parens operator")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          double* start = &pa(0, 0, 0);
          for (unsigned char ci = 0; ci < num_components; ci++) {
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK(&pa(xi, yi, ci) == start + xi + yi * (nx + 2 * num_ghost) + ci * (nx + 2 * num_ghost) * (ny + 2 * num_ghost));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray parens operator const")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          const PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          const double* start = &pa(0, 0, 0);
          for (unsigned char ci = 0; ci < num_components; ci++) {
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK(&pa(xi, yi, ci) == start + xi + yi * (nx + 2 * num_ghost) + ci * (nx + 2 * num_ghost) * (ny + 2 * num_ghost));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray default is zero")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          for (unsigned char ci = 0; ci < num_components; ci++) {
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK(pa(xi, yi, ci) == 0.0);
              }
            }
          }
          CHECK(pa.getStrides()[0] == 1);
          CHECK(pa.getStrides()[1] == nx + 2 * num_ghost);
        }
      }
    }
  }
}
TEST_CASE("PatchArray<2> getSliceOn<1>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            View<double, 2> slice = pa.getSliceOn(Side<2>::west(), { xi });
            for (unsigned char ci = 0; ci < num_components; ci++) {
              INFO("ci: " << ci);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&pa(xi, yi, ci) == &slice(yi, ci));
              }
            }
          }

          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            View<double, 2> slice = pa.getSliceOn(Side<2>::east(), { xi });
            for (unsigned char ci = 0; ci < num_components; ci++) {
              INFO("ci: " << ci);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&pa(nx - 1 - xi, yi, ci) == &slice(yi, ci));
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            View<double, 2> slice = pa.getSliceOn(Side<2>::south(), { yi });
            for (unsigned char ci = 0; ci < num_components; ci++) {
              INFO("ci: " << ci);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&pa(xi, yi, ci) == &slice(xi, ci));
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            View<double, 2> slice = pa.getSliceOn(Side<2>::north(), { yi });
            for (unsigned char ci = 0; ci < num_components; ci++) {
              INFO("ci: " << ci);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&pa(xi, ny - 1 - yi, ci) == &slice(xi, ci));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<2> getSliceOn<1> const")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          const PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            View<const double, 2> slice = pa.getSliceOn(Side<2>::west(), { xi });
            for (unsigned char ci = 0; ci < num_components; ci++) {
              INFO("ci: " << ci);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&pa(xi, yi, ci) == &slice(yi, ci));
              }
            }
          }

          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            INFO("xi: " << xi);
            View<const double, 2> slice = pa.getSliceOn(Side<2>::east(), { xi });
            for (unsigned char ci = 0; ci < num_components; ci++) {
              INFO("ci: " << ci);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&pa(nx - 1 - xi, yi, ci) == &slice(yi, ci));
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            View<const double, 2> slice = pa.getSliceOn(Side<2>::south(), { yi });
            for (unsigned char ci = 0; ci < num_components; ci++) {
              INFO("ci: " << ci);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&pa(xi, yi, ci) == &slice(xi, ci));
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            View<const double, 2> slice = pa.getSliceOn(Side<2>::north(), { yi });
            for (unsigned char ci = 0; ci < num_components; ci++) {
              INFO("ci: " << ci);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&pa(xi, ny - 1 - yi, ci) == &slice(xi, ci));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<3> getSliceOn<2>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_components : { 1, 2 }) {
          for (auto num_ghost : { 0, 1, 2 }) {

            PatchArray<3> pa({ nx, ny, nz }, num_components, num_ghost);

            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 3> slice = pa.getSliceOn(Side<3>::west(), { xi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                  INFO("zi: " << zi);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(xi, yi, zi, ci) == &slice(yi, zi, ci));
                  }
                }
              }
            }

            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 3> slice = pa.getSliceOn(Side<3>::east(), { xi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                  INFO("zi: " << zi);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(nx - 1 - xi, yi, zi, ci) == &slice(yi, zi, ci));
                  }
                }
              }
            }

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              View<double, 3> slice = pa.getSliceOn(Side<3>::south(), { yi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                  INFO("zi: " << zi);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, yi, zi, ci) == &slice(xi, zi, ci));
                  }
                }
              }
            }

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              View<double, 3> slice = pa.getSliceOn(Side<3>::north(), { yi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                  INFO("zi: " << zi);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, ny - 1 - yi, zi, ci) == &slice(xi, zi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              View<double, 3> slice = pa.getSliceOn(Side<3>::bottom(), { zi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                  INFO("yi: " << yi);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, yi, zi, ci) == &slice(xi, yi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              View<double, 3> slice = pa.getSliceOn(Side<3>::top(), { zi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                  INFO("yi: " << yi);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, yi, nz - 1 - zi, ci) == &slice(xi, yi, ci));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<3> getSliceOn<2> const")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_components : { 1, 2 }) {
          for (auto num_ghost : { 0, 1, 2 }) {

            const PatchArray<3> pa({ nx, ny, nz }, num_components, num_ghost);

            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<const double, 3> slice = pa.getSliceOn(Side<3>::west(), { xi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                  INFO("zi: " << zi);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(xi, yi, zi, ci) == &slice(yi, zi, ci));
                  }
                }
              }
            }

            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              View<const double, 3> slice = pa.getSliceOn(Side<3>::east(), { xi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                  INFO("zi: " << zi);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(nx - 1 - xi, yi, zi, ci) == &slice(yi, zi, ci));
                  }
                }
              }
            }

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              View<const double, 3> slice = pa.getSliceOn(Side<3>::south(), { yi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                  INFO("zi: " << zi);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, yi, zi, ci) == &slice(xi, zi, ci));
                  }
                }
              }
            }

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              View<const double, 3> slice = pa.getSliceOn(Side<3>::north(), { yi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                  INFO("zi: " << zi);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, ny - 1 - yi, zi, ci) == &slice(xi, zi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              View<const double, 3> slice = pa.getSliceOn(Side<3>::bottom(), { zi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                  INFO("yi: " << yi);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, yi, zi, ci) == &slice(xi, yi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              View<const double, 3> slice = pa.getSliceOn(Side<3>::top(), { zi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                  INFO("yi: " << yi);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, yi, nz - 1 - zi, ci) == &slice(xi, yi, ci));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<3> getSliceOn<1>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_components : { 1, 2 }) {
          for (auto num_ghost : { 0, 1, 2 }) {

            PatchArray<3> pa({ nx, ny, nz }, num_components, num_ghost);

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                View<double, 2> slice = pa.getSliceOn(Edge::bs(), { yi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, yi, zi, ci) == &slice(xi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                View<double, 2> slice = pa.getSliceOn(Edge::tn(), { yi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, ny - 1 - yi, nz - 1 - zi, ci) == &slice(xi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                View<double, 2> slice = pa.getSliceOn(Edge::bn(), { yi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, ny - 1 - yi, zi, ci) == &slice(xi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                View<double, 2> slice = pa.getSliceOn(Edge::ts(), { yi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, yi, nz - 1 - zi, ci) == &slice(xi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<double, 2> slice = pa.getSliceOn(Edge::bw(), { xi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(xi, yi, zi, ci) == &slice(yi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<double, 2> slice = pa.getSliceOn(Edge::te(), { xi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(nx - 1 - xi, yi, nz - 1 - zi, ci) == &slice(yi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<double, 2> slice = pa.getSliceOn(Edge::be(), { xi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(nx - 1 - xi, yi, zi, ci) == &slice(yi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<double, 2> slice = pa.getSliceOn(Edge::tw(), { xi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(xi, yi, nz - 1 - zi, ci) == &slice(yi, ci));
                  }
                }
              }
            }

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<double, 2> slice = pa.getSliceOn(Edge::sw(), { xi, yi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                    INFO("zi: " << zi);
                    CHECK(&pa(xi, yi, zi, ci) == &slice(zi, ci));
                  }
                }
              }
            }

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<double, 2> slice = pa.getSliceOn(Edge::ne(), { xi, yi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                    INFO("zi: " << zi);
                    CHECK(&pa(nx - 1 - xi, ny - 1 - yi, zi, ci) == &slice(zi, ci));
                  }
                }
              }
            }

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<double, 2> slice = pa.getSliceOn(Edge::se(), { xi, yi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                    INFO("zi: " << zi);
                    CHECK(&pa(nx - 1 - xi, yi, zi, ci) == &slice(zi, ci));
                  }
                }
              }
            }

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<double, 2> slice = pa.getSliceOn(Edge::nw(), { xi, yi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                    INFO("zi: " << zi);
                    CHECK(&pa(xi, ny - 1 - yi, zi, ci) == &slice(zi, ci));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<3> getSliceOn<1> const")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_components : { 1, 2 }) {
          for (auto num_ghost : { 0, 1, 2 }) {

            const PatchArray<3> pa({ nx, ny, nz }, num_components, num_ghost);

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                View<const double, 2> slice = pa.getSliceOn(Edge::bs(), { yi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, yi, zi, ci) == &slice(xi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                View<const double, 2> slice = pa.getSliceOn(Edge::tn(), { yi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, ny - 1 - yi, nz - 1 - zi, ci) == &slice(xi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                View<const double, 2> slice = pa.getSliceOn(Edge::bn(), { yi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, ny - 1 - yi, zi, ci) == &slice(xi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                View<const double, 2> slice = pa.getSliceOn(Edge::ts(), { yi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, yi, nz - 1 - zi, ci) == &slice(xi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<const double, 2> slice = pa.getSliceOn(Edge::bw(), { xi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(xi, yi, zi, ci) == &slice(yi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<const double, 2> slice = pa.getSliceOn(Edge::te(), { xi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(nx - 1 - xi, yi, nz - 1 - zi, ci) == &slice(yi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<const double, 2> slice = pa.getSliceOn(Edge::be(), { xi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(nx - 1 - xi, yi, zi, ci) == &slice(yi, ci));
                  }
                }
              }
            }

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<const double, 2> slice = pa.getSliceOn(Edge::tw(), { xi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(xi, yi, nz - 1 - zi, ci) == &slice(yi, ci));
                  }
                }
              }
            }

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<const double, 2> slice = pa.getSliceOn(Edge::sw(), { xi, yi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                    INFO("zi: " << zi);
                    CHECK(&pa(xi, yi, zi, ci) == &slice(zi, ci));
                  }
                }
              }
            }

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<const double, 2> slice = pa.getSliceOn(Edge::ne(), { xi, yi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                    INFO("zi: " << zi);
                    CHECK(&pa(nx - 1 - xi, ny - 1 - yi, zi, ci) == &slice(zi, ci));
                  }
                }
              }
            }

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<const double, 2> slice = pa.getSliceOn(Edge::se(), { xi, yi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                    INFO("zi: " << zi);
                    CHECK(&pa(nx - 1 - xi, yi, zi, ci) == &slice(zi, ci));
                  }
                }
              }
            }

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              INFO("yi: " << yi);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                View<const double, 2> slice = pa.getSliceOn(Edge::nw(), { xi, yi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                    INFO("zi: " << zi);
                    CHECK(&pa(xi, ny - 1 - yi, zi, ci) == &slice(zi, ci));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<3> getSliceOn<0>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_components : { 1, 2 }) {
          for (auto num_ghost : { 0, 1, 2 }) {

            PatchArray<3> pa({ nx, ny, nz }, num_components, num_ghost);

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                  INFO("xi: " << xi);
                  for (unsigned char ci = 0; ci < num_components; ci++) {
                    INFO("ci: " << ci);
                    CHECK(&pa(xi, yi, zi, ci) == &(pa.getSliceOn(Corner<3>::bsw(), { xi, yi, zi }))(ci));
                    CHECK(&pa(nx - 1 - xi, yi, zi, ci) == &(pa.getSliceOn(Corner<3>::bse(), { xi, yi, zi })(ci)));
                    CHECK(&pa(xi, ny - 1 - yi, zi, ci) == &(pa.getSliceOn(Corner<3>::bnw(), { xi, yi, zi })(ci)));
                    CHECK(&pa(nx - 1 - xi, ny - 1 - yi, zi, ci) == &(pa.getSliceOn(Corner<3>::bne(), { xi, yi, zi })(ci)));
                    CHECK(&pa(xi, yi, nz - 1 - zi, ci) == &(pa.getSliceOn(Corner<3>::tsw(), { xi, yi, zi })(ci)));
                    CHECK(&pa(nx - 1 - xi, yi, nz - 1 - zi, ci) == &(pa.getSliceOn(Corner<3>::tse(), { xi, yi, zi })(ci)));
                    CHECK(&pa(xi, ny - 1 - yi, nz - 1 - zi, ci) == &(pa.getSliceOn(Corner<3>::tnw(), { xi, yi, zi })(ci)));
                    CHECK(&pa(nx - 1 - xi, ny - 1 - yi, nz - 1 - zi, ci) == &(pa.getSliceOn(Corner<3>::tne(), { xi, yi, zi })(ci)));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<3> getSliceOn<0> const")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_components : { 1, 2 }) {
          for (auto num_ghost : { 0, 1, 2 }) {

            const PatchArray<3> pa({ nx, ny, nz }, num_components, num_ghost);

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              INFO("zi: " << zi);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                  INFO("xi: " << xi);
                  for (unsigned char ci = 0; ci < num_components; ci++) {
                    INFO("ci: " << ci);
                    CHECK(&pa(xi, yi, zi, ci) == &(pa.getSliceOn(Corner<3>::bsw(), { xi, yi, zi }))(ci));
                    CHECK(&pa(nx - 1 - xi, yi, zi, ci) == &(pa.getSliceOn(Corner<3>::bse(), { xi, yi, zi })(ci)));
                    CHECK(&pa(xi, ny - 1 - yi, zi, ci) == &(pa.getSliceOn(Corner<3>::bnw(), { xi, yi, zi })(ci)));
                    CHECK(&pa(nx - 1 - xi, ny - 1 - yi, zi, ci) == &(pa.getSliceOn(Corner<3>::bne(), { xi, yi, zi })(ci)));
                    CHECK(&pa(xi, yi, nz - 1 - zi, ci) == &(pa.getSliceOn(Corner<3>::tsw(), { xi, yi, zi })(ci)));
                    CHECK(&pa(nx - 1 - xi, yi, nz - 1 - zi, ci) == &(pa.getSliceOn(Corner<3>::tse(), { xi, yi, zi })(ci)));
                    CHECK(&pa(xi, ny - 1 - yi, nz - 1 - zi, ci) == &(pa.getSliceOn(Corner<3>::tnw(), { xi, yi, zi })(ci)));
                    CHECK(&pa(nx - 1 - xi, ny - 1 - yi, nz - 1 - zi, ci) == &(pa.getSliceOn(Corner<3>::tne(), { xi, yi, zi })(ci)));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<2> getSliceOn<0>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                CHECK(&pa(xi, yi, ci) == &(pa.getSliceOn(Corner<2>::sw(), { xi, yi })(ci)));
                CHECK(&pa(nx - 1 - xi, yi, ci) == &(pa.getSliceOn(Corner<2>::se(), { xi, yi })(ci)));
                CHECK(&pa(xi, ny - 1 - yi, ci) == &(pa.getSliceOn(Corner<2>::nw(), { xi, yi })(ci)));
                CHECK(&pa(nx - 1 - xi, ny - 1 - yi, ci) == &(pa.getSliceOn(Corner<2>::ne(), { xi, yi })(ci)));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<2> getSliceOn<0> const")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          const PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            INFO("yi: " << yi);
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              INFO("xi: " << xi);
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                CHECK(&pa(xi, yi, ci) == &(pa.getSliceOn(Corner<2>::sw(), { xi, yi })(ci)));
                CHECK(&pa(nx - 1 - xi, yi, ci) == &(pa.getSliceOn(Corner<2>::se(), { xi, yi })(ci)));
                CHECK(&pa(xi, ny - 1 - yi, ci) == &(pa.getSliceOn(Corner<2>::nw(), { xi, yi })(ci)));
                CHECK(&pa(nx - 1 - xi, ny - 1 - yi, ci) == &(pa.getSliceOn(Corner<2>::ne(), { xi, yi })(ci)));
              }
            }
          }
        }
      }
    }
  }
}

TEST_CASE("PatchArray<2> getGhostSliceOn<1>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          for (unsigned char xi = 0; xi < num_ghost; xi++) {
            INFO("xi: " << xi);
            View<double, 2> slice_w = pa.getGhostSliceOn(Side<2>::west(), { xi });
            View<double, 2> slice_e = pa.getGhostSliceOn(Side<2>::east(), { xi });
            for (unsigned char ci = 0; ci < num_components; ci++) {
              INFO("ci: " << ci);
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                INFO("yi: " << yi);
                CHECK(&pa(-1 - xi, yi, ci) == &slice_w(yi, ci));
                CHECK(&pa(nx + xi, yi, ci) == &slice_e(yi, ci));
              }
            }
          }

          for (unsigned char yi = 0; yi < num_ghost; yi++) {
            INFO("yi: " << yi);
            View<double, 2> slice_s = pa.getGhostSliceOn(Side<2>::south(), { yi });
            View<double, 2> slice_n = pa.getGhostSliceOn(Side<2>::north(), { yi });
            for (unsigned char ci = 0; ci < num_components; ci++) {
              INFO("ci: " << ci);
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                INFO("xi: " << xi);
                CHECK(&pa(xi, -1 - yi, ci) == &slice_s(xi, ci));
                CHECK(&pa(xi, ny + yi, ci) == &slice_n(xi, ci));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<3> getGhostSliceOn<2>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_components : { 1, 2 }) {
          for (auto num_ghost : { 0, 1, 2 }) {

            const PatchArray<3> pa({ nx, ny, nz }, num_components, num_ghost);

            for (unsigned char xi = 0; xi < num_ghost; xi++) {
              INFO("xi: " << xi);
              View<double, 3> slice_w = pa.getGhostSliceOn(Side<3>::west(), { xi });
              View<double, 3> slice_e = pa.getGhostSliceOn(Side<3>::east(), { xi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                  INFO("zi: " << zi);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(-1 - xi, yi, zi, ci) == &slice_w(yi, zi, ci));
                    CHECK(&pa(nx + xi, yi, zi, ci) == &slice_e(yi, zi, ci));
                  }
                }
              }
            }

            for (unsigned char yi = 0; yi < num_ghost; yi++) {
              INFO("yi: " << yi);
              View<double, 3> slice_s = pa.getGhostSliceOn(Side<3>::south(), { yi });
              View<double, 3> slice_n = pa.getGhostSliceOn(Side<3>::north(), { yi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                  INFO("zi: " << zi);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, -1 - yi, zi, ci) == &slice_s(xi, zi, ci));
                    CHECK(&pa(xi, ny + yi, zi, ci) == &slice_n(xi, zi, ci));
                  }
                }
              }
            }

            for (unsigned char zi = 0; zi < num_ghost; zi++) {
              INFO("zi: " << zi);
              View<double, 3> slice_b = pa.getGhostSliceOn(Side<3>::bottom(), { zi });
              View<double, 3> slice_t = pa.getGhostSliceOn(Side<3>::top(), { zi });
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                  INFO("yi: " << yi);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, yi, -1 - zi, ci) == &slice_b(xi, yi, ci));
                    CHECK(&pa(xi, yi, nz + zi, ci) == &slice_t(xi, yi, ci));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<3> getGhostSliceOn<1>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_components : { 1, 2 }) {
          for (auto num_ghost : { 0, 1, 2 }) {

            PatchArray<3> pa({ nx, ny, nz }, num_components, num_ghost);

            for (unsigned char zi = 0; zi < num_ghost; zi++) {
              INFO("zi: " << zi);
              for (unsigned char yi = 0; yi < num_ghost; yi++) {
                INFO("yi: " << yi);
                View<double, 2> slice_bs = pa.getGhostSliceOn(Edge::bs(), { yi, zi });
                View<double, 2> slice_tn = pa.getGhostSliceOn(Edge::tn(), { yi, zi });
                View<double, 2> slice_bn = pa.getGhostSliceOn(Edge::bn(), { yi, zi });
                View<double, 2> slice_ts = pa.getGhostSliceOn(Edge::ts(), { yi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                    INFO("xi: " << xi);
                    CHECK(&pa(xi, -1 - yi, -1 - zi, ci) == &slice_bs(xi, ci));
                    CHECK(&pa(xi, ny + yi, nz + zi, ci) == &slice_tn(xi, ci));
                    CHECK(&pa(xi, ny + yi, -1 - zi, ci) == &slice_bn(xi, ci));
                    CHECK(&pa(xi, -1 - yi, nz + zi, ci) == &slice_ts(xi, ci));
                  }
                }
              }
            }

            for (unsigned char zi = 0; zi < num_ghost; zi++) {
              INFO("zi: " << zi);
              for (unsigned char xi = 0; xi < num_ghost; xi++) {
                INFO("xi: " << xi);
                View<double, 2> slice_bw = pa.getGhostSliceOn(Edge::bw(), { xi, zi });
                View<double, 2> slice_te = pa.getGhostSliceOn(Edge::te(), { xi, zi });
                View<double, 2> slice_be = pa.getGhostSliceOn(Edge::be(), { xi, zi });
                View<double, 2> slice_tw = pa.getGhostSliceOn(Edge::tw(), { xi, zi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                    INFO("yi: " << yi);
                    CHECK(&pa(-1 - xi, yi, -1 - zi, ci) == &slice_bw(yi, ci));
                    CHECK(&pa(nx + xi, yi, nz + zi, ci) == &slice_te(yi, ci));
                    CHECK(&pa(nx + xi, yi, -1 - zi, ci) == &slice_be(yi, ci));
                    CHECK(&pa(-1 - xi, yi, nz + zi, ci) == &slice_tw(yi, ci));
                  }
                }
              }
            }

            for (unsigned char yi = 0; yi < num_ghost; yi++) {
              INFO("yi: " << yi);
              for (unsigned char xi = 0; xi < num_ghost; xi++) {
                INFO("xi: " << xi);
                View<double, 2> slice_sw = pa.getGhostSliceOn(Edge::sw(), { xi, yi });
                View<double, 2> slice_ne = pa.getGhostSliceOn(Edge::ne(), { xi, yi });
                View<double, 2> slice_se = pa.getGhostSliceOn(Edge::se(), { xi, yi });
                View<double, 2> slice_nw = pa.getGhostSliceOn(Edge::nw(), { xi, yi });
                for (unsigned char ci = 0; ci < num_components; ci++) {
                  INFO("ci: " << ci);
                  for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                    INFO("zi: " << zi);
                    CHECK(&pa(-1 - xi, -1 - yi, zi, ci) == &slice_sw(zi, ci));
                    CHECK(&pa(nx + xi, ny + yi, zi, ci) == &slice_ne(zi, ci));
                    CHECK(&pa(nx + xi, -1 - yi, zi, ci) == &slice_se(zi, ci));
                    CHECK(&pa(-1 - xi, ny + yi, zi, ci) == &slice_nw(zi, ci));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<3> getGhostSliceOn<0>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_components : { 1, 2 }) {
          for (auto num_ghost : { 0, 1, 2 }) {

            PatchArray<3> pa({ nx, ny, nz }, num_components, num_ghost);

            for (unsigned char zi = 0; zi < num_ghost; zi++) {
              INFO("zi: " << zi);
              for (unsigned char yi = 0; yi < num_ghost; yi++) {
                INFO("yi: " << yi);
                for (unsigned char xi = 0; xi < num_ghost; xi++) {
                  INFO("xi: " << xi);
                  for (unsigned char ci = 0; ci < num_components; ci++) {
                    INFO("ci: " << ci);
                    CHECK(&pa(-1 - xi, -1 - yi, -1 - zi, ci) == &(pa.getGhostSliceOn(Corner<3>::bsw(), { xi, yi, zi })(ci)));
                    CHECK(&pa(nx + xi, -1 - yi, -1 - zi, ci) == &(pa.getGhostSliceOn(Corner<3>::bse(), { xi, yi, zi })(ci)));
                    CHECK(&pa(-1 - xi, ny + yi, -1 - zi, ci) == &(pa.getGhostSliceOn(Corner<3>::bnw(), { xi, yi, zi })(ci)));
                    CHECK(&pa(nx + xi, ny + yi, -1 - zi, ci) == &(pa.getGhostSliceOn(Corner<3>::bne(), { xi, yi, zi })(ci)));
                    CHECK(&pa(-1 - xi, -1 - yi, nz + zi, ci) == &(pa.getGhostSliceOn(Corner<3>::tsw(), { xi, yi, zi })(ci)));
                    CHECK(&pa(nx + xi, -1 - yi, nz + zi, ci) == &(pa.getGhostSliceOn(Corner<3>::tse(), { xi, yi, zi })(ci)));
                    CHECK(&pa(-1 - xi, ny + yi, nz + zi, ci) == &(pa.getGhostSliceOn(Corner<3>::tnw(), { xi, yi, zi })(ci)));
                    CHECK(&pa(nx + xi, ny + yi, nz + zi, ci) == &(pa.getGhostSliceOn(Corner<3>::tne(), { xi, yi, zi })(ci)));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray<2> getGhostSliceOn<0>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          for (unsigned char yi = 0; yi < num_ghost; yi++) {
            INFO("yi: " << yi);
            for (unsigned char xi = 0; xi < num_ghost; xi++) {
              INFO("xi: " << xi);
              for (unsigned char ci = 0; ci < num_components; ci++) {
                INFO("ci: " << ci);
                CHECK(&pa(-xi - 1, -yi - 1, ci) == &(pa.getGhostSliceOn(Corner<2>::sw(), { xi, yi })(ci)));
                CHECK(&pa(nx + xi, -yi - 1, ci) == &(pa.getGhostSliceOn(Corner<2>::se(), { xi, yi })(ci)));
                CHECK(&pa(-xi - 1, ny + yi, ci) == &(pa.getGhostSliceOn(Corner<2>::nw(), { xi, yi })(ci)));
                CHECK(&pa(nx + xi, ny + yi, ci) == &(pa.getGhostSliceOn(Corner<2>::ne(), { xi, yi })(ci)));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("PatchArray copy constructor")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          for (int c = 0; c < num_components; c++) {
            for (int iy = -num_ghost; iy < ny + num_ghost; iy++) {
              for (int ix = -num_ghost; ix < nx + num_ghost; ix++) {
                pa(ix, iy, c) = ix + iy + c;
              }
            }
          }

          PatchArray<2> pa_copy(pa);

          for (int c = 0; c < num_components; c++) {
            for (int iy = -num_ghost; iy < ny + num_ghost; iy++) {
              for (int ix = -num_ghost; ix < nx + num_ghost; ix++) {
                CHECK(pa(ix, iy, c) == pa_copy(ix, iy, c));
              }
            }
          }

          CHECK(pa.getStrides() == pa_copy.getStrides());
          CHECK(pa.getGhostStart() == pa_copy.getGhostStart());
          CHECK(pa.getStart() == pa_copy.getStart());
          CHECK(pa.getEnd() == pa_copy.getEnd());
          CHECK(pa.getGhostEnd() == pa_copy.getGhostEnd());
          CHECK(&pa(0, 0, 0) != &pa_copy(0, 0, 0));
        }
      }
    }
  }
}
TEST_CASE("PatchArray copy assignment")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_components : { 1, 2 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          PatchArray<2> pa({ nx, ny }, num_components, num_ghost);

          for (int c = 0; c < num_components; c++) {
            for (int iy = -num_ghost; iy < ny + num_ghost; iy++) {
              for (int ix = -num_ghost; ix < nx + num_ghost; ix++) {
                pa(ix, iy, c) = ix + iy + c;
              }
            }
          }

          PatchArray<2> pa_copy({ nx, ny }, num_components, num_ghost);
          pa_copy = pa;

          for (int c = 0; c < num_components; c++) {
            for (int iy = -num_ghost; iy < ny + num_ghost; iy++) {
              for (int ix = -num_ghost; ix < nx + num_ghost; ix++) {
                CHECK(pa(ix, iy, c) == pa_copy(ix, iy, c));
              }
            }
          }

          CHECK(pa.getStrides() == pa_copy.getStrides());
          CHECK(pa.getGhostStart() == pa_copy.getGhostStart());
          CHECK(pa.getStart() == pa_copy.getStart());
          CHECK(pa.getEnd() == pa_copy.getEnd());
          CHECK(pa.getGhostEnd() == pa_copy.getGhostEnd());
          CHECK(&pa(0, 0, 0) != &pa_copy(0, 0, 0));
        }
      }
    }
  }
}
