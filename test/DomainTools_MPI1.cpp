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
#include <ThunderEgg/DomainTools.h>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;
using namespace doctest;

TEST_CASE("DomainTools::GetRealCoord 1D")
{
  for (auto nx : { 1, 2, 3, 10, 16, 17 }) {
    for (auto startx : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto lengthx : { 1.0, 0.23, 3.0 }) {
        PatchInfo<1> pinfo;

        pinfo.ns = { nx };
        pinfo.spacings = { lengthx / nx };
        pinfo.starts = { startx };

        for (int coordx = -1; coordx <= nx; coordx++) {
          array<int, 1> coord = { coordx };
          array<double, 1> expected;
          if (coordx == -1) {
            expected[0] = startx;
          } else if (coordx == nx) {
            expected[0] = startx + lengthx;
          } else {
            expected[0] = startx + (0.5 + coordx) * lengthx / nx;
          }

          array<double, 1> result;
          DomainTools::GetRealCoord<1>(pinfo, coord, result);

          CHECK(result[0] + 100 == Approx(expected[0] + 100));
        }
      }
    }
  }
}
TEST_CASE("DomainTools::GetRealCoord 2D")
{
  for (auto ny : { 1, 2, 3 }) {
    for (auto starty : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto lengthy : { 1.0, 0.23, 3.0 }) {
        PatchInfo<2> pinfo;

        auto nx = 2;
        auto startx = 0.0;
        auto lengthx = 1.0;

        pinfo.ns = { nx, ny };
        pinfo.spacings = { lengthx / nx, lengthy / ny };
        pinfo.starts = { startx, starty };

        for (int coordy = -1; coordy <= ny; coordy++) {
          for (int coordx = -1; coordx <= nx; coordx++) {
            array<double, 2> expected;
            if (coordx == -1) {
              expected[0] = startx;
            } else if (coordx == nx) {
              expected[0] = startx + lengthx;
            } else {
              expected[0] = startx + (0.5 + coordx) * lengthx / nx;
            }
            if (coordy == -1) {
              expected[1] = starty;
            } else if (coordy == ny) {
              expected[1] = starty + lengthy;
            } else {
              expected[1] = starty + (0.5 + coordy) * lengthy / ny;
            }

            array<double, 2> result;
            array<int, 2> coord = { coordx, coordy };
            DomainTools::GetRealCoord<2>(pinfo, coord, result);

            CHECK(result[0] + 100 == Approx(expected[0] + 100));
            CHECK(result[1] + 100 == Approx(expected[1] + 100));
          }
        }
      }
    }
  }
}
TEST_CASE("DomainTools::GetRealCoord 3D")
{
  for (auto nz : { 1, 2, 3 }) {
    for (auto startz : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto lengthz : { 1.0, 0.23, 3.0 }) {
        PatchInfo<3> pinfo;

        auto nx = 2;
        auto ny = 2;
        auto startx = 0.0;
        auto starty = 0.0;
        auto lengthx = 1.0;
        auto lengthy = 1.0;

        pinfo.ns = { nx, ny, nz };
        pinfo.spacings = { lengthx / nx, lengthy / ny, lengthz / nz };
        pinfo.starts = { startx, starty, startz };

        for (int coordz = -1; coordz <= nz; coordz++) {
          for (int coordy = -1; coordy <= ny; coordy++) {
            for (int coordx = -1; coordx <= nx; coordx++) {
              array<double, 3> expected;
              if (coordx == -1) {
                expected[0] = startx;
              } else if (coordx == nx) {
                expected[0] = startx + lengthx;
              } else {
                expected[0] = startx + (0.5 + coordx) * lengthx / nx;
              }
              if (coordy == -1) {
                expected[1] = starty;
              } else if (coordy == ny) {
                expected[1] = starty + lengthy;
              } else {
                expected[1] = starty + (0.5 + coordy) * lengthy / ny;
              }
              if (coordz == -1) {
                expected[2] = startz;
              } else if (coordz == nz) {
                expected[2] = startz + lengthz;
              } else {
                expected[2] = startz + (0.5 + coordz) * lengthz / nz;
              }

              array<double, 3> result;
              array<int, 3> coord = { coordx, coordy, coordz };
              DomainTools::GetRealCoord<3>(pinfo, coord, result);

              CHECK(result[0] + 100 == Approx(expected[0] + 100));
              CHECK(result[1] + 100 == Approx(expected[1] + 100));
              CHECK(result[2] + 100 == Approx(expected[2] + 100));
            }
          }
        }
      }
    }
  }
}
TEST_CASE("DomainTools::getRealCoordGhost 1D")
{
  for (auto nx : { 1, 2, 3, 10, 16, 17 }) {
    for (auto startx : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto lengthx : { 1.0, 0.23, 3.0 }) {
        PatchInfo<1> pinfo;

        pinfo.ns = { nx };
        pinfo.spacings = { lengthx / nx };
        pinfo.starts = { startx };

        for (int coordx = -1; coordx <= nx; coordx++) {
          array<int, 1> coord = { coordx };
          array<double, 1> expected;
          expected[0] = startx + (0.5 + coordx) * lengthx / nx;

          array<double, 1> result;
          DomainTools::GetRealCoordGhost<1>(pinfo, coord, result);

          CHECK(result[0] + 100 == Approx(expected[0] + 100));
        }
      }
    }
  }
}
TEST_CASE("DomainTools::getRealCoordGhost 2D")
{
  for (auto ny : { 1, 2, 3 }) {
    for (auto starty : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto lengthy : { 1.0, 0.23, 3.0 }) {
        PatchInfo<2> pinfo;

        auto nx = 2;
        auto startx = 0.0;
        auto lengthx = 1.0;

        pinfo.ns = { nx, ny };
        pinfo.spacings = { lengthx / nx, lengthy / ny };
        pinfo.starts = { startx, starty };

        for (int coordy = -1; coordy <= ny; coordy++) {
          for (int coordx = -1; coordx <= nx; coordx++) {
            array<double, 2> expected;
            expected[0] = startx + (0.5 + coordx) * lengthx / nx;
            expected[1] = starty + (0.5 + coordy) * lengthy / ny;

            array<double, 2> result;
            array<int, 2> coord = { coordx, coordy };
            DomainTools::GetRealCoordGhost<2>(pinfo, coord, result);

            CHECK(result[0] + 100 == Approx(expected[0] + 100));
            CHECK(result[1] + 100 == Approx(expected[1] + 100));
          }
        }
      }
    }
  }
}
TEST_CASE("DomainTools::getRealCoordGhost 3D")
{
  for (auto nz : { 1, 2, 3 }) {
    for (auto startz : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto lengthz : { 1.0, 0.23, 3.0 }) {
        PatchInfo<3> pinfo;

        auto nx = 2;
        auto ny = 2;
        auto startx = 0.0;
        auto starty = 0.0;
        auto lengthx = 1.0;
        auto lengthy = 1.0;

        pinfo.ns = { nx, ny, nz };
        pinfo.spacings = { lengthx / nx, lengthy / ny, lengthz / nz };
        pinfo.starts = { startx, starty, startz };

        for (int coordz = -1; coordz <= nz; coordz++) {
          for (int coordy = -1; coordy <= ny; coordy++) {
            for (int coordx = -1; coordx <= nx; coordx++) {
              array<double, 3> expected;
              expected[0] = startx + (0.5 + coordx) * lengthx / nx;
              expected[1] = starty + (0.5 + coordy) * lengthy / ny;
              expected[2] = startz + (0.5 + coordz) * lengthz / nz;

              array<double, 3> result;
              array<int, 3> coord = { coordx, coordy, coordz };
              DomainTools::GetRealCoordGhost<3>(pinfo, coord, result);

              CHECK(result[0] + 100 == Approx(expected[0] + 100));
              CHECK(result[1] + 100 == Approx(expected[1] + 100));
              CHECK(result[2] + 100 == Approx(expected[2] + 100));
            }
          }
        }
      }
    }
  }
}
TEST_CASE("DomainTools::GetRealCoordBound 1D")
{
  for (auto nx : { 1, 2, 3, 10, 16, 17 }) {
    for (auto startx : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto lengthx : { 1.0, 0.23, 3.0 }) {
        PatchInfo<1> pinfo;

        pinfo.ns = { nx };
        pinfo.spacings = { lengthx / nx };
        pinfo.starts = { startx };

        std::array<int, 0> coord;
        std::array<double, 1> result;
        DomainTools::GetRealCoordBound<1>(pinfo, coord, Side<1>::west(), result);
        CHECK(result[0] + 100 == Approx(startx + 100));
        DomainTools::GetRealCoordBound<1>(pinfo, coord, Side<1>::east(), result);
        CHECK(result[0] + 100 == Approx(startx + lengthx + 100));
      }
    }
  }
}
TEST_CASE("DomainTools::GetRealCoordBound 2D")
{
  for (auto ny : { 1, 2, 3 }) {
    for (auto starty : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto lengthy : { 1.0, 0.23, 3.0 }) {
        PatchInfo<2> pinfo;

        auto nx = 2;
        auto startx = 0.0;
        auto lengthx = 1.0;

        pinfo.ns = { nx, ny };
        pinfo.spacings = { lengthx / nx, lengthy / ny };
        pinfo.starts = { startx, starty };

        for (int coordx = -1; coordx <= nx; coordx++) {
          array<double, 2> expected;
          if (coordx == -1) {
            expected[0] = startx;
          } else if (coordx == nx) {
            expected[0] = startx + lengthx;
          } else {
            expected[0] = startx + (0.5 + coordx) * lengthx / nx;
          }

          array<double, 2> result;
          array<int, 1> coord = { coordx };
          DomainTools::GetRealCoordBound<2>(pinfo, coord, Side<2>::south(), result);

          CHECK(result[0] + 100 == Approx(expected[0] + 100));
          CHECK(result[1] + 100 == Approx(starty + 100));

          DomainTools::GetRealCoordBound<2>(pinfo, coord, Side<2>::north(), result);
          CHECK(result[0] + 100 == Approx(expected[0] + 100));
          CHECK(result[1] + 100 == Approx(starty + lengthy + 100));
        }
        for (int coordy = -1; coordy <= ny; coordy++) {
          array<double, 2> expected;
          if (coordy == -1) {
            expected[1] = starty;
          } else if (coordy == ny) {
            expected[1] = starty + lengthy;
          } else {
            expected[1] = starty + (0.5 + coordy) * lengthy / ny;
          }

          array<double, 2> result;
          array<int, 1> coord = { coordy };
          DomainTools::GetRealCoordBound<2>(pinfo, coord, Side<2>::west(), result);

          CHECK(result[0] + 100 == Approx(startx + 100));
          CHECK(result[1] + 100 == Approx(expected[1] + 100));

          DomainTools::GetRealCoordBound<2>(pinfo, coord, Side<2>::east(), result);

          CHECK(result[0] + 100 == Approx(startx + lengthx + 100));
          CHECK(result[1] + 100 == Approx(expected[1] + 100));
        }
      }
    }
  }
}
TEST_CASE("DomainTools::setValues 1D g=x")
{
  for (auto nx : { 1, 2, 10, 13 }) {
    for (auto startx : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto spacingx : { 0.01, 1.0, 3.14 }) {
        for (auto num_ghost : { 0, 1, 2, 3, 4, 5 }) {
          Communicator comm(MPI_COMM_WORLD);

          auto f = [](const std::array<double, 1> coord) { return coord[0]; };

          vector<PatchInfo<1>> pinfos(1);

          pinfos[0].id = 0;
          pinfos[0].ns = { nx };
          pinfos[0].spacings = { spacingx };
          pinfos[0].starts = { startx };
          pinfos[0].num_ghost_cells = num_ghost;
          Domain<1> d(comm, 0, { nx }, num_ghost, pinfos.begin(), pinfos.end());

          Vector<1> vec(d, 1);

          DomainTools::SetValues<1>(d, vec, f);
          auto ld = vec.getComponentView(0, 0);
          Loop::Nested<1>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 1> coord) {
            if (coord[0] < 0 || coord[0] >= nx) {
              CHECK(ld[coord] + 100 == Approx(0.0 + 100));
            } else {
              std::array<double, 1> real_coord;
              DomainTools::GetRealCoord<1>(pinfos[0], coord, real_coord);
              CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
            }
          });
        }
      }
    }
  }
}
TEST_CASE("DomainTools::setValues 1D f=x**2")
{
  for (auto nx : { 1, 2, 10, 13 }) {
    for (auto startx : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto spacingx : { 0.01, 1.0, 3.14 }) {
        for (auto num_ghost : { 0, 1, 2, 3, 4, 5 }) {
          Communicator comm(MPI_COMM_WORLD);

          auto f = [](const std::array<double, 1> coord) { return coord[0] * coord[0]; };

          vector<PatchInfo<1>> pinfos(1);

          pinfos[0].id = 0;
          pinfos[0].ns = { nx };
          pinfos[0].spacings = { spacingx };
          pinfos[0].starts = { startx };
          pinfos[0].num_ghost_cells = num_ghost;
          Domain<1> d(comm, 0, { nx }, num_ghost, pinfos.begin(), pinfos.end());

          Vector<1> vec(d, 1);

          DomainTools::SetValues<1>(d, vec, f);
          auto ld = vec.getComponentView(0, 0);
          Loop::Nested<1>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 1> coord) {
            if (coord[0] < 0 || coord[0] >= nx) {
              CHECK(ld[coord] + 100 == Approx(0.0 + 100));
            } else {
              std::array<double, 1> real_coord;
              DomainTools::GetRealCoord<1>(pinfos[0], coord, real_coord);
              CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
            }
          });
        }
      }
    }
  }
}
TEST_CASE("DomainTools::setValues 2D f=x+y")
{
  for (auto ny : { 1, 2, 10, 13 }) {
    for (auto starty : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto spacingy : { 0.01, 1.0, 3.14 }) {
        for (auto num_ghost : { 0, 1, 2, 3, 4, 5 }) {
          Communicator comm(MPI_COMM_WORLD);

          auto f = [](const std::array<double, 2> coord) { return coord[0] + coord[1]; };

          vector<PatchInfo<2>> pinfos(1);

          int nx = 3;
          double startx = 0.0;
          double spacingx = 0.1;

          pinfos[0].id = 0;
          pinfos[0].ns = { nx, ny };
          pinfos[0].spacings = { spacingx, spacingy };
          pinfos[0].starts = { startx, starty };
          pinfos[0].num_ghost_cells = num_ghost;
          Domain<2> d(comm, 0, { nx, ny }, num_ghost, pinfos.begin(), pinfos.end());

          Vector<2> vec(d, 1);

          DomainTools::SetValues<2>(d, vec, f);
          auto ld = vec.getComponentView(0, 0);
          Loop::Nested<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
            if (coord[0] < 0 || coord[0] >= nx || coord[1] < 0 || coord[1] >= ny) {
              CHECK(ld[coord] + 100 == Approx(0.0 + 100));
            } else {
              std::array<double, 2> real_coord;
              DomainTools::GetRealCoord<2>(pinfos[0], coord, real_coord);
              CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
            }
          });
        }
      }
    }
  }
}
TEST_CASE("DomainTools::setValues 2D f=x+y,g=x*y")
{
  for (auto ny : { 1, 2, 10, 13 }) {
    for (auto starty : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto spacingy : { 0.01, 1.0, 3.14 }) {
        for (auto num_ghost : { 0, 1, 2, 3, 4, 5 }) {
          Communicator comm(MPI_COMM_WORLD);

          auto f = [](const std::array<double, 2> coord) { return coord[0] + coord[1]; };
          auto g = [](const std::array<double, 2> coord) { return coord[0] * coord[1]; };

          vector<PatchInfo<2>> pinfos(1);

          int nx = 3;
          double startx = 0.0;
          double spacingx = 0.1;

          pinfos[0].id = 0;
          pinfos[0].ns = { nx, ny };
          pinfos[0].spacings = { spacingx, spacingy };
          pinfos[0].starts = { startx, starty };
          pinfos[0].num_ghost_cells = num_ghost;
          Domain<2> d(comm, 0, { nx, ny }, num_ghost, pinfos.begin(), pinfos.end());

          Vector<2> vec(d, 2);

          DomainTools::SetValues<2>(d, vec, f, g);
          auto ld = vec.getComponentView(0, 0);
          Loop::Nested<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
            if (coord[0] < 0 || coord[0] >= nx || coord[1] < 0 || coord[1] >= ny) {
              CHECK(ld[coord] + 100 == Approx(0.0 + 100));
            } else {
              std::array<double, 2> real_coord;
              DomainTools::GetRealCoord<2>(pinfos[0], coord, real_coord);
              CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
            }
          });
          auto ld2 = vec.getComponentView(1, 0);
          Loop::Nested<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
            if (coord[0] < 0 || coord[0] >= nx || coord[1] < 0 || coord[1] >= ny) {
              CHECK(ld2[coord] + 100 == Approx(0.0 + 100));
            } else {
              std::array<double, 2> real_coord;
              DomainTools::GetRealCoord<2>(pinfos[0], coord, real_coord);
              CHECK(ld2[coord] + 100 == Approx(g(real_coord) + 100));
            }
          });
        }
      }
    }
  }
}
TEST_CASE("DomainTools::setValues 2D f=x*y")
{
  for (auto ny : { 1, 2, 10, 13 }) {
    for (auto starty : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto spacingy : { 0.01, 1.0, 3.14 }) {
        for (auto num_ghost : { 0, 1, 2, 3, 4, 5 }) {
          Communicator comm(MPI_COMM_WORLD);

          auto f = [](const std::array<double, 2> coord) { return coord[0] * coord[1]; };

          vector<PatchInfo<2>> pinfos(1);

          int nx = 3;
          double startx = 0.0;
          double spacingx = 0.1;

          pinfos[0].id = 0;
          pinfos[0].ns = { nx, ny };
          pinfos[0].spacings = { spacingx, spacingy };
          pinfos[0].starts = { startx, starty };
          pinfos[0].num_ghost_cells = num_ghost;
          Domain<2> d(comm, 0, { nx, ny }, num_ghost, pinfos.begin(), pinfos.end());

          Vector<2> vec(d, 1);

          DomainTools::SetValues<2>(d, vec, f);
          auto ld = vec.getComponentView(0, 0);
          Loop::Nested<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
            if (coord[0] < 0 || coord[0] >= nx || coord[1] < 0 || coord[1] >= ny) {
              CHECK(ld[coord] + 100 == Approx(0.0 + 100));
            } else {
              std::array<double, 2> real_coord;
              DomainTools::GetRealCoord<2>(pinfos[0], coord, real_coord);
              CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
            }
          });
        }
      }
    }
  }
}
TEST_CASE("DomainTools::setValuesWithGhost 1D f=x")
{
  for (auto nx : { 1, 2, 10, 13 }) {
    for (auto startx : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto spacingx : { 0.01, 1.0, 3.14 }) {
        for (auto num_ghost : { 0, 1, 2, 3, 4, 5 }) {
          Communicator comm(MPI_COMM_WORLD);

          auto f = [](const std::array<double, 1> coord) { return coord[0]; };

          vector<PatchInfo<1>> pinfos(1);

          pinfos[0].id = 0;
          pinfos[0].ns = { nx };
          pinfos[0].spacings = { spacingx };
          pinfos[0].starts = { startx };
          pinfos[0].num_ghost_cells = num_ghost;
          Domain<1> d(comm, 0, { nx }, num_ghost, pinfos.begin(), pinfos.end());

          Vector<1> vec(d, 1);

          DomainTools::SetValuesWithGhost<1>(d, vec, f);
          auto ld = vec.getComponentView(0, 0);
          Loop::Nested<1>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 1> coord) {
            std::array<double, 1> real_coord;
            DomainTools::GetRealCoordGhost<1>(pinfos[0], coord, real_coord);
            CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
          });
        }
      }
    }
  }
}
TEST_CASE("DomainTools::setValuesWithGhost 1D f=x**2")
{
  for (auto nx : { 1, 2, 10, 13 }) {
    for (auto startx : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto spacingx : { 0.01, 1.0, 3.14 }) {
        for (auto num_ghost : { 0, 1, 2, 3, 4, 5 }) {
          Communicator comm(MPI_COMM_WORLD);

          auto f = [](const std::array<double, 1> coord) { return coord[0] * coord[0]; };

          vector<PatchInfo<1>> pinfos(1);

          pinfos[0].id = 0;
          pinfos[0].ns = { nx };
          pinfos[0].spacings = { spacingx };
          pinfos[0].starts = { startx };
          pinfos[0].num_ghost_cells = num_ghost;
          Domain<1> d(comm, 0, { nx }, num_ghost, pinfos.begin(), pinfos.end());

          Vector<1> vec(d, 1);

          DomainTools::SetValuesWithGhost<1>(d, vec, f);
          auto ld = vec.getComponentView(0, 0);
          Loop::Nested<1>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 1> coord) {
            std::array<double, 1> real_coord;
            DomainTools::GetRealCoordGhost<1>(pinfos[0], coord, real_coord);
            CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
          });
        }
      }
    }
  }
}
TEST_CASE("DomainTools::setValuesWithGhost 2D f=x+y")
{
  for (auto ny : { 1, 2, 10, 13 }) {
    for (auto starty : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto spacingy : { 0.01, 1.0, 3.14 }) {
        for (auto num_ghost : { 0, 1, 2, 3, 4, 5 }) {
          Communicator comm(MPI_COMM_WORLD);

          auto f = [](const std::array<double, 2> coord) { return coord[0] + coord[1]; };

          vector<PatchInfo<2>> pinfos(1);

          int nx = 3;
          double startx = 0.0;
          double spacingx = 0.1;

          pinfos[0].id = 0;
          pinfos[0].ns = { nx, ny };
          pinfos[0].spacings = { spacingx, spacingy };
          pinfos[0].starts = { startx, starty };
          pinfos[0].num_ghost_cells = num_ghost;
          Domain<2> d(comm, 0, { nx, ny }, num_ghost, pinfos.begin(), pinfos.end());

          Vector<2> vec(d, 1);

          DomainTools::SetValuesWithGhost<2>(d, vec, f);
          auto ld = vec.getComponentView(0, 0);
          Loop::Nested<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
            std::array<double, 2> real_coord;
            DomainTools::GetRealCoordGhost<2>(pinfos[0], coord, real_coord);
            CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
          });
        }
      }
    }
  }
}
TEST_CASE("DomainTools::setValuesWithGhost 2D f=x+y,g=x*y")
{
  for (auto ny : { 1, 2, 10, 13 }) {
    for (auto starty : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto spacingy : { 0.01, 1.0, 3.14 }) {
        for (auto num_ghost : { 0, 1, 2, 3, 4, 5 }) {
          Communicator comm(MPI_COMM_WORLD);

          auto f = [](const std::array<double, 2> coord) { return coord[0] + coord[1]; };
          auto g = [](const std::array<double, 2> coord) { return coord[0] * coord[1]; };

          vector<PatchInfo<2>> pinfos(1);

          int nx = 3;
          double startx = 0.0;
          double spacingx = 0.1;

          pinfos[0].id = 0;
          pinfos[0].ns = { nx, ny };
          pinfos[0].spacings = { spacingx, spacingy };
          pinfos[0].starts = { startx, starty };
          pinfos[0].num_ghost_cells = num_ghost;
          Domain<2> d(comm, 0, { nx, ny }, num_ghost, pinfos.begin(), pinfos.end());

          Vector<2> vec(d, 2);

          DomainTools::SetValuesWithGhost<2>(d, vec, f, g);
          auto ld = vec.getComponentView(0, 0);
          Loop::Nested<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
            std::array<double, 2> real_coord;
            DomainTools::GetRealCoordGhost<2>(pinfos[0], coord, real_coord);
            CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
          });
          auto ld2 = vec.getComponentView(1, 0);
          Loop::Nested<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
            std::array<double, 2> real_coord;
            DomainTools::GetRealCoordGhost<2>(pinfos[0], coord, real_coord);
            CHECK(ld2[coord] + 100 == Approx(g(real_coord) + 100));
          });
        }
      }
    }
  }
}
TEST_CASE("DomainTools::setValuesWithGhost 2D throws when too many functions are given")
{
  for (auto ny : { 1, 2, 10, 13 }) {
    for (auto starty : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto spacingy : { 0.01, 1.0, 3.14 }) {
        for (auto num_ghost : { 0, 1, 2, 3, 4, 5 }) {
          Communicator comm(MPI_COMM_WORLD);

          auto f = [](const std::array<double, 2> coord) { return coord[0] + coord[1]; };
          auto g = [](const std::array<double, 2> coord) { return coord[0] * coord[1]; };

          vector<PatchInfo<2>> pinfos(1);

          int nx = 3;
          double startx = 0.0;
          double spacingx = 0.1;

          pinfos[0].id = 0;
          pinfos[0].ns = { nx, ny };
          pinfos[0].spacings = { spacingx, spacingy };
          pinfos[0].starts = { startx, starty };
          pinfos[0].num_ghost_cells = num_ghost;
          Domain<2> d(comm, 0, { nx, ny }, num_ghost, pinfos.begin(), pinfos.end());

          Vector<2> vec(d, 1);

          CHECK_THROWS_AS(DomainTools::SetValuesWithGhost<2>(d, vec, f, g), RuntimeError);
        }
      }
    }
  }
}
TEST_CASE("DomainTools::setValues 2D throws when too many functions are given")
{
  for (auto ny : { 1, 2, 10, 13 }) {
    for (auto starty : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto spacingy : { 0.01, 1.0, 3.14 }) {
        for (auto num_ghost : { 0, 1, 2, 3, 4, 5 }) {
          Communicator comm(MPI_COMM_WORLD);

          auto f = [](const std::array<double, 2> coord) { return coord[0] + coord[1]; };
          auto g = [](const std::array<double, 2> coord) { return coord[0] * coord[1]; };

          vector<PatchInfo<2>> pinfos(1);

          int nx = 3;
          double startx = 0.0;
          double spacingx = 0.1;

          pinfos[0].id = 0;
          pinfos[0].ns = { nx, ny };
          pinfos[0].spacings = { spacingx, spacingy };
          pinfos[0].starts = { startx, starty };
          pinfos[0].num_ghost_cells = num_ghost;
          Domain<2> d(comm, 0, { nx, ny }, num_ghost, pinfos.begin(), pinfos.end());

          Vector<2> vec(d, 1);

          CHECK_THROWS_AS(DomainTools::SetValues<2>(d, vec, f, g), RuntimeError);
        }
      }
    }
  }
}
TEST_CASE("DomainTools::setValuesWithGhost 2D f=x*y")
{
  for (auto ny : { 1, 2, 10, 13 }) {
    for (auto starty : { 0.0, -1.0, 1.0, 0.23, -0.23, 3.0, -3.0 }) {
      for (auto spacingy : { 0.01, 1.0, 3.14 }) {
        for (auto num_ghost : { 0, 1, 2, 3, 4, 5 }) {
          Communicator comm(MPI_COMM_WORLD);

          auto f = [](const std::array<double, 2> coord) { return coord[0] * coord[1]; };

          vector<PatchInfo<2>> pinfos(1);

          int nx = 3;
          double startx = 0.0;
          double spacingx = 0.1;

          pinfos[0].id = 0;
          pinfos[0].ns = { nx, ny };
          pinfos[0].spacings = { spacingx, spacingy };
          pinfos[0].starts = { startx, starty };
          pinfos[0].num_ghost_cells = num_ghost;
          Domain<2> d(comm, 0, { nx, ny }, num_ghost, pinfos.begin(), pinfos.end());

          Vector<2> vec(d, 1);

          DomainTools::SetValuesWithGhost<2>(d, vec, f);
          auto ld = vec.getComponentView(0, 0);
          Loop::Nested<2>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, 2> coord) {
            std::array<double, 2> real_coord;
            DomainTools::GetRealCoordGhost<2>(pinfos[0], coord, real_coord);
            CHECK(ld[coord] + 100 == Approx(f(real_coord) + 100));
          });
        }
      }
    }
  }
}
