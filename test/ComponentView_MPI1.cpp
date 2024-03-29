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
#include <ThunderEgg/Vector.h>
#include <doctest.h>

using namespace std;
using namespace ThunderEgg;
TEST_CASE("ComponentView constructor")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        array<int, 2> lengths = { nx, ny };
        array<int, 2> strides = { 1, nx + 2 * num_ghost };
        int size = 1;
        for (size_t i = 0; i < 2; i++) {
          size *= (lengths[i] + 2 * num_ghost);
        }
        vector<double> vec(size);
        iota(vec.begin(), vec.end(), 0);

        ComponentView<double, 2> ld(vec.data(), strides, lengths, num_ghost);

        for (int i = 0; i < 2; i++) {
          CHECK_EQ(ld.getStrides()[i], strides[i]);
          CHECK_EQ(ld.getStart()[i], 0);
          CHECK_EQ(ld.getEnd()[i], lengths[i] - 1);
          CHECK_EQ(ld.getGhostStart()[i], -num_ghost);
          CHECK_EQ(ld.getGhostEnd()[i], lengths[i] - 1 + num_ghost);
        }
      }
    }
  }
}
TEST_CASE("ComponentView<double,2> getSliceOn<1>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        array<int, 2> lengths = { nx, ny };
        array<int, 2> strides = { 1, nx + 2 * num_ghost };
        int size = 1;
        for (size_t i = 0; i < 2; i++) {
          size *= (lengths[i] + 2 * num_ghost);
        }
        double data[size];

        ComponentView<double, 2> ld(data, strides, lengths, num_ghost);

        for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
          View<double, 1> slice = ld.getSliceOn(Side<2>::west(), { xi });
          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            CHECK_EQ(&ld[{ xi, yi }], &slice[{ yi }]);
          }
        }

        for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
          View<double, 1> slice = ld.getSliceOn(Side<2>::east(), { xi });
          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            CHECK_EQ(&ld[{ nx - 1 - xi, yi }], &slice[{ yi }]);
          }
        }

        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          View<double, 1> slice = ld.getSliceOn(Side<2>::south(), { yi });
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            CHECK_EQ(&ld[{ xi, yi }], &slice[{ xi }]);
          }
        }

        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          View<double, 1> slice = ld.getSliceOn(Side<2>::north(), { yi });
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            CHECK_EQ(&ld[{ xi, ny - 1 - yi }], &slice[{ xi }]);
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView<double,3> getSliceOn<2>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          array<int, 3> lengths = { nx, ny, nz };
          array<int, 3> strides = { 1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost) };
          int size = 1;
          for (size_t i = 0; i < 3; i++) {
            size *= (lengths[i] + 2 * num_ghost);
          }
          double data[size];

          ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            View<double, 2> slice = ld.getSliceOn(Side<3>::west(), { xi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                CHECK_EQ(&ld[{ xi, yi, zi }], &slice[{ yi, zi }]);
              }
            }
          }

          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            View<double, 2> slice = ld.getSliceOn(Side<3>::east(), { xi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                CHECK_EQ(&ld[{ nx - 1 - xi, yi, zi }], &slice[{ yi, zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            View<double, 2> slice = ld.getSliceOn(Side<3>::south(), { yi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, yi, zi }], &slice[{ xi, zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            View<double, 2> slice = ld.getSliceOn(Side<3>::north(), { yi });
            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, ny - 1 - yi, zi }], &slice[{ xi, zi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            View<double, 2> slice = ld.getSliceOn(Side<3>::bottom(), { zi });
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, yi, zi }], &slice[{ xi, yi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            View<double, 2> slice = ld.getSliceOn(Side<3>::top(), { zi });
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, yi, nz - 1 - zi }], &slice[{ xi, yi }]);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView<double,3> getSliceOn<1>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          array<int, 3> lengths = { nx, ny, nz };
          array<int, 3> strides = { 1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost) };
          int size = 1;
          for (size_t i = 0; i < 3; i++) {
            size *= (lengths[i] + 2 * num_ghost);
          }
          double data[size];

          ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              View<double, 1> slice = ld.getSliceOn(Edge::bs(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, yi, zi }], &slice[{ xi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              View<double, 1> slice = ld.getSliceOn(Edge::tn(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, ny - 1 - yi, nz - 1 - zi }], &slice[{ xi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              View<double, 1> slice = ld.getSliceOn(Edge::bn(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, ny - 1 - yi, zi }], &slice[{ xi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              View<double, 1> slice = ld.getSliceOn(Edge::ts(), { yi, zi });
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, yi, nz - 1 - zi }], &slice[{ xi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              View<double, 1> slice = ld.getSliceOn(Edge::bw(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                CHECK_EQ(&ld[{ xi, yi, zi }], &slice[{ yi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              View<double, 1> slice = ld.getSliceOn(Edge::te(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                CHECK_EQ(&ld[{ nx - 1 - xi, yi, nz - 1 - zi }], &slice[{ yi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              View<double, 1> slice = ld.getSliceOn(Edge::be(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                CHECK_EQ(&ld[{ nx - 1 - xi, yi, zi }], &slice[{ yi }]);
              }
            }
          }

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              View<double, 1> slice = ld.getSliceOn(Edge::tw(), { xi, zi });
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                CHECK_EQ(&ld[{ xi, yi, nz - 1 - zi }], &slice[{ yi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              View<double, 1> slice = ld.getSliceOn(Edge::sw(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                CHECK_EQ(&ld[{ xi, yi, zi }], &slice[{ zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              View<double, 1> slice = ld.getSliceOn(Edge::ne(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                CHECK_EQ(&ld[{ nx - 1 - xi, ny - 1 - yi, zi }], &slice[{ zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              View<double, 1> slice = ld.getSliceOn(Edge::se(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                CHECK_EQ(&ld[{ nx - 1 - xi, yi, zi }], &slice[{ zi }]);
              }
            }
          }

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
              View<double, 1> slice = ld.getSliceOn(Edge::nw(), { xi, yi });
              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                CHECK_EQ(&ld[{ xi, ny - 1 - yi, zi }], &slice[{ zi }]);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView<double,3> getSliceOn<0>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          array<int, 3> lengths = { nx, ny, nz };
          array<int, 3> strides = { 1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost) };
          int size = 1;
          for (size_t i = 0; i < 3; i++) {
            size *= (lengths[i] + 2 * num_ghost);
          }
          double data[size];

          ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, yi, zi }], &(ld.getSliceOn(Corner<3>::bsw(), { xi, yi, zi }))[{}]);
                CHECK_EQ(&ld[{ nx - 1 - xi, yi, zi }], &(ld.getSliceOn(Corner<3>::bse(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ xi, ny - 1 - yi, zi }], &(ld.getSliceOn(Corner<3>::bnw(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ nx - 1 - xi, ny - 1 - yi, zi }], &(ld.getSliceOn(Corner<3>::bne(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ xi, yi, nz - 1 - zi }], &(ld.getSliceOn(Corner<3>::tsw(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ nx - 1 - xi, yi, nz - 1 - zi }], &(ld.getSliceOn(Corner<3>::tse(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ xi, ny - 1 - yi, nz - 1 - zi }], &(ld.getSliceOn(Corner<3>::tnw(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ nx - 1 - xi, ny - 1 - yi, nz - 1 - zi }], &(ld.getSliceOn(Corner<3>::tne(), { xi, yi, zi })[{}]));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView<double,3> getSliceOn<0> const")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          array<int, 3> lengths = { nx, ny, nz };
          array<int, 3> strides = { 1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost) };
          int size = 1;
          for (size_t i = 0; i < 3; i++) {
            size *= (lengths[i] + 2 * num_ghost);
          }
          double data[size];

          const ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

          for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, yi, zi }], &(ld.getSliceOn(Corner<3>::bsw(), { xi, yi, zi }))[{}]);
                CHECK_EQ(&ld[{ nx - 1 - xi, yi, zi }], &(ld.getSliceOn(Corner<3>::bse(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ xi, ny - 1 - yi, zi }], &(ld.getSliceOn(Corner<3>::bnw(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ nx - 1 - xi, ny - 1 - yi, zi }], &(ld.getSliceOn(Corner<3>::bne(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ xi, yi, nz - 1 - zi }], &(ld.getSliceOn(Corner<3>::tsw(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ nx - 1 - xi, yi, nz - 1 - zi }], &(ld.getSliceOn(Corner<3>::tse(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ xi, ny - 1 - yi, nz - 1 - zi }], &(ld.getSliceOn(Corner<3>::tnw(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ nx - 1 - xi, ny - 1 - yi, nz - 1 - zi }], &(ld.getSliceOn(Corner<3>::tne(), { xi, yi, zi })[{}]));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView<double,2> getSliceOn<0>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        array<int, 2> lengths = { nx, ny };
        array<int, 2> strides = { 1, nx + 2 * num_ghost };
        int size = 1;
        for (size_t i = 0; i < 2; i++) {
          size *= (lengths[i] + 2 * num_ghost);
        }
        double data[size];

        ComponentView<double, 2> ld(data, strides, lengths, num_ghost);

        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            CHECK_EQ(&ld[{ xi, yi }], &(ld.getSliceOn(Corner<2>::sw(), { xi, yi })[{}]));
            CHECK_EQ(&ld[{ nx - 1 - xi, yi }], &(ld.getSliceOn(Corner<2>::se(), { xi, yi })[{}]));
            CHECK_EQ(&ld[{ xi, ny - 1 - yi }], &(ld.getSliceOn(Corner<2>::nw(), { xi, yi })[{}]));
            CHECK_EQ(&ld[{ nx - 1 - xi, ny - 1 - yi }], &(ld.getSliceOn(Corner<2>::ne(), { xi, yi })[{}]));
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView<double,2> getSliceOn<0> const")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        array<int, 2> lengths = { nx, ny };
        array<int, 2> strides = { 1, nx + 2 * num_ghost };
        int size = 1;
        for (size_t i = 0; i < 2; i++) {
          size *= (lengths[i] + 2 * num_ghost);
        }
        double data[size];

        const ComponentView<double, 2> ld(data, strides, lengths, num_ghost);

        for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            CHECK_EQ(&ld[{ xi, yi }], &(ld.getSliceOn(Corner<2>::sw(), { xi, yi })[{}]));
            CHECK_EQ(&ld[{ nx - 1 - xi, yi }], &(ld.getSliceOn(Corner<2>::se(), { xi, yi })[{}]));
            CHECK_EQ(&ld[{ xi, ny - 1 - yi }], &(ld.getSliceOn(Corner<2>::nw(), { xi, yi })[{}]));
            CHECK_EQ(&ld[{ nx - 1 - xi, ny - 1 - yi }], &(ld.getSliceOn(Corner<2>::ne(), { xi, yi })[{}]));
          }
        }
      }
    }
  }
}

TEST_CASE("ComponentView<double,2> getGhostSliceOn<1>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        array<int, 2> lengths = { nx, ny };
        array<int, 2> strides = { 1, nx + 2 * num_ghost };
        int size = 1;
        for (size_t i = 0; i < 2; i++) {
          size *= (lengths[i] + 2 * num_ghost);
        }
        double data[size];

        const ComponentView<double, 2> ld(data, strides, lengths, num_ghost);

        for (unsigned char xi = 0; xi < num_ghost; xi++) {
          View<double, 1> slice_w = ld.getGhostSliceOn(Side<2>::west(), { xi });
          CHECK_EQ(slice_w.getGhostStart(), std::array<int, 1>({ -num_ghost }));
          CHECK_EQ(slice_w.getStart(), std::array<int, 1>({ 0 }));
          CHECK_EQ(slice_w.getEnd(), std::array<int, 1>({ ny - 1 }));
          CHECK_EQ(slice_w.getGhostEnd(), std::array<int, 1>({ ny - 1 + num_ghost }));

          View<double, 1> slice_e = ld.getGhostSliceOn(Side<2>::east(), { xi });
          CHECK_EQ(slice_e.getGhostStart(), std::array<int, 1>({ -num_ghost }));
          CHECK_EQ(slice_e.getStart(), std::array<int, 1>({ 0 }));
          CHECK_EQ(slice_e.getEnd(), std::array<int, 1>({ ny - 1 }));
          CHECK_EQ(slice_e.getGhostEnd(), std::array<int, 1>({ ny - 1 + num_ghost }));

          for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
            CHECK_EQ(&ld[{ -1 - xi, yi }], &slice_w[{ yi }]);
            CHECK_EQ(&ld[{ nx + xi, yi }], &slice_e[{ yi }]);
          }
        }

        for (unsigned char yi = 0; yi < num_ghost; yi++) {
          View<double, 1> slice_s = ld.getGhostSliceOn(Side<2>::south(), { yi });
          CHECK_EQ(slice_s.getGhostStart(), std::array<int, 1>({ -num_ghost }));
          CHECK_EQ(slice_s.getStart(), std::array<int, 1>({ 0 }));
          CHECK_EQ(slice_s.getEnd(), std::array<int, 1>({ nx - 1 }));
          CHECK_EQ(slice_s.getGhostEnd(), std::array<int, 1>({ nx - 1 + num_ghost }));

          View<double, 1> slice_n = ld.getGhostSliceOn(Side<2>::north(), { yi });
          CHECK_EQ(slice_s.getGhostStart(), std::array<int, 1>({ -num_ghost }));
          CHECK_EQ(slice_s.getStart(), std::array<int, 1>({ 0 }));
          CHECK_EQ(slice_s.getEnd(), std::array<int, 1>({ nx - 1 }));
          CHECK_EQ(slice_s.getGhostEnd(), std::array<int, 1>({ nx - 1 + num_ghost }));

          for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
            CHECK_EQ(&ld[{ xi, -1 - yi }], &slice_s[{ xi }]);
            CHECK_EQ(&ld[{ xi, ny + yi }], &slice_n[{ xi }]);
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView<double,3> getGhostSliceOn<2>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          array<int, 3> lengths = { nx, ny, nz };
          array<int, 3> strides = { 1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost) };
          int size = 1;
          for (size_t i = 0; i < 3; i++) {
            size *= (lengths[i] + 2 * num_ghost);
          }
          double data[size];

          const ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

          for (unsigned char xi = 0; xi < num_ghost; xi++) {
            View<double, 2> slice_w = ld.getGhostSliceOn(Side<3>::west(), { xi });
            CHECK_EQ(slice_w.getGhostStart(), std::array<int, 2>({ -num_ghost, -num_ghost }));
            CHECK_EQ(slice_w.getStart(), std::array<int, 2>({ 0, 0 }));
            CHECK_EQ(slice_w.getEnd(), std::array<int, 2>({ ny - 1, nz - 1 }));
            CHECK_EQ(slice_w.getGhostEnd(), std::array<int, 2>({ ny - 1 + num_ghost, nz - 1 + num_ghost }));

            View<double, 2> slice_e = ld.getGhostSliceOn(Side<3>::east(), { xi });
            CHECK_EQ(slice_e.getGhostStart(), std::array<int, 2>({ -num_ghost, -num_ghost }));
            CHECK_EQ(slice_e.getStart(), std::array<int, 2>({ 0, 0 }));
            CHECK_EQ(slice_e.getEnd(), std::array<int, 2>({ ny - 1, nz - 1 }));
            CHECK_EQ(slice_e.getGhostEnd(), std::array<int, 2>({ ny - 1 + num_ghost, nz - 1 + num_ghost }));

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                CHECK_EQ(&ld[{ -1 - xi, yi, zi }], &slice_w[{ yi, zi }]);
                CHECK_EQ(&ld[{ nx + xi, yi, zi }], &slice_e[{ yi, zi }]);
              }
            }
          }

          for (unsigned char yi = 0; yi < num_ghost; yi++) {
            View<double, 2> slice_s = ld.getGhostSliceOn(Side<3>::south(), { yi });
            CHECK_EQ(slice_s.getGhostStart(), std::array<int, 2>({ -num_ghost, -num_ghost }));
            CHECK_EQ(slice_s.getStart(), std::array<int, 2>({ 0, 0 }));
            CHECK_EQ(slice_s.getEnd(), std::array<int, 2>({ nx - 1, nz - 1 }));
            CHECK_EQ(slice_s.getGhostEnd(), std::array<int, 2>({ nx - 1 + num_ghost, nz - 1 + num_ghost }));

            View<double, 2> slice_n = ld.getGhostSliceOn(Side<3>::north(), { yi });
            CHECK_EQ(slice_n.getGhostStart(), std::array<int, 2>({ -num_ghost, -num_ghost }));
            CHECK_EQ(slice_n.getStart(), std::array<int, 2>({ 0, 0 }));
            CHECK_EQ(slice_n.getEnd(), std::array<int, 2>({ nx - 1, nz - 1 }));
            CHECK_EQ(slice_n.getGhostEnd(), std::array<int, 2>({ nx - 1 + num_ghost, nz - 1 + num_ghost }));

            for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, -1 - yi, zi }], &slice_s[{ xi, zi }]);
                CHECK_EQ(&ld[{ xi, ny + yi, zi }], &slice_n[{ xi, zi }]);
              }
            }
          }

          for (unsigned char zi = 0; zi < num_ghost; zi++) {
            View<double, 2> slice_b = ld.getGhostSliceOn(Side<3>::bottom(), { zi });
            CHECK_EQ(slice_b.getGhostStart(), std::array<int, 2>({ -num_ghost, -num_ghost }));
            CHECK_EQ(slice_b.getStart(), std::array<int, 2>({ 0, 0 }));
            CHECK_EQ(slice_b.getEnd(), std::array<int, 2>({ nx - 1, ny - 1 }));
            CHECK_EQ(slice_b.getGhostEnd(), std::array<int, 2>({ nx - 1 + num_ghost, ny - 1 + num_ghost }));

            View<double, 2> slice_t = ld.getGhostSliceOn(Side<3>::top(), { zi });
            CHECK_EQ(slice_t.getGhostStart(), std::array<int, 2>({ -num_ghost, -num_ghost }));
            CHECK_EQ(slice_t.getStart(), std::array<int, 2>({ 0, 0 }));
            CHECK_EQ(slice_t.getEnd(), std::array<int, 2>({ nx - 1, ny - 1 }));
            CHECK_EQ(slice_t.getGhostEnd(), std::array<int, 2>({ nx - 1 + num_ghost, ny - 1 + num_ghost }));

            for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, yi, -1 - zi }], &slice_b[{ xi, yi }]);
                CHECK_EQ(&ld[{ xi, yi, nz + zi }], &slice_t[{ xi, yi }]);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView<double,3> getGhostSliceOn<1>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          array<int, 3> lengths = { nx, ny, nz };
          array<int, 3> strides = { 1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost) };
          int size = 1;
          for (size_t i = 0; i < 3; i++) {
            size *= (lengths[i] + 2 * num_ghost);
          }
          double data[size];

          const ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

          for (unsigned char zi = 0; zi < num_ghost; zi++) {
            for (unsigned char yi = 0; yi < num_ghost; yi++) {
              View<double, 1> slice_bs = ld.getGhostSliceOn(Edge::bs(), { yi, zi });
              CHECK_EQ(slice_bs.getGhostStart(), std::array<int, 1>({ -num_ghost }));
              CHECK_EQ(slice_bs.getStart(), std::array<int, 1>({ 0 }));
              CHECK_EQ(slice_bs.getEnd(), std::array<int, 1>({ nx - 1 }));
              CHECK_EQ(slice_bs.getGhostEnd(), std::array<int, 1>({ nx - 1 + num_ghost }));

              View<double, 1> slice_tn = ld.getGhostSliceOn(Edge::tn(), { yi, zi });
              CHECK_EQ(slice_tn.getGhostStart(), std::array<int, 1>({ -num_ghost }));
              CHECK_EQ(slice_tn.getStart(), std::array<int, 1>({ 0 }));
              CHECK_EQ(slice_tn.getEnd(), std::array<int, 1>({ nx - 1 }));
              CHECK_EQ(slice_tn.getGhostEnd(), std::array<int, 1>({ nx - 1 + num_ghost }));

              View<double, 1> slice_bn = ld.getGhostSliceOn(Edge::bn(), { yi, zi });
              CHECK_EQ(slice_bn.getGhostStart(), std::array<int, 1>({ -num_ghost }));
              CHECK_EQ(slice_bn.getStart(), std::array<int, 1>({ 0 }));
              CHECK_EQ(slice_bn.getEnd(), std::array<int, 1>({ nx - 1 }));
              CHECK_EQ(slice_bn.getGhostEnd(), std::array<int, 1>({ nx - 1 + num_ghost }));

              View<double, 1> slice_ts = ld.getGhostSliceOn(Edge::ts(), { yi, zi });
              CHECK_EQ(slice_ts.getGhostStart(), std::array<int, 1>({ -num_ghost }));
              CHECK_EQ(slice_ts.getStart(), std::array<int, 1>({ 0 }));
              CHECK_EQ(slice_ts.getEnd(), std::array<int, 1>({ nx - 1 }));
              CHECK_EQ(slice_ts.getGhostEnd(), std::array<int, 1>({ nx - 1 + num_ghost }));

              for (int xi = -num_ghost; xi < nx + num_ghost; xi++) {
                CHECK_EQ(&ld[{ xi, -1 - yi, -1 - zi }], &slice_bs[{ xi }]);
                CHECK_EQ(&ld[{ xi, ny + yi, nz + zi }], &slice_tn[{ xi }]);
                CHECK_EQ(&ld[{ xi, ny + yi, -1 - zi }], &slice_bn[{ xi }]);
                CHECK_EQ(&ld[{ xi, -1 - yi, nz + zi }], &slice_ts[{ xi }]);
              }
            }
          }

          for (unsigned char zi = 0; zi < num_ghost; zi++) {
            for (unsigned char xi = 0; xi < num_ghost; xi++) {
              View<double, 1> slice_bw = ld.getGhostSliceOn(Edge::bw(), { xi, zi });
              CHECK_EQ(slice_bw.getGhostStart(), std::array<int, 1>({ -num_ghost }));
              CHECK_EQ(slice_bw.getStart(), std::array<int, 1>({ 0 }));
              CHECK_EQ(slice_bw.getEnd(), std::array<int, 1>({ ny - 1 }));
              CHECK_EQ(slice_bw.getGhostEnd(), std::array<int, 1>({ ny - 1 + num_ghost }));

              View<double, 1> slice_te = ld.getGhostSliceOn(Edge::te(), { xi, zi });
              CHECK_EQ(slice_te.getGhostStart(), std::array<int, 1>({ -num_ghost }));
              CHECK_EQ(slice_te.getStart(), std::array<int, 1>({ 0 }));
              CHECK_EQ(slice_te.getEnd(), std::array<int, 1>({ ny - 1 }));
              CHECK_EQ(slice_te.getGhostEnd(), std::array<int, 1>({ ny - 1 + num_ghost }));

              View<double, 1> slice_be = ld.getGhostSliceOn(Edge::be(), { xi, zi });
              CHECK_EQ(slice_be.getGhostStart(), std::array<int, 1>({ -num_ghost }));
              CHECK_EQ(slice_be.getStart(), std::array<int, 1>({ 0 }));
              CHECK_EQ(slice_be.getEnd(), std::array<int, 1>({ ny - 1 }));
              CHECK_EQ(slice_be.getGhostEnd(), std::array<int, 1>({ ny - 1 + num_ghost }));

              View<double, 1> slice_tw = ld.getGhostSliceOn(Edge::tw(), { xi, zi });
              CHECK_EQ(slice_tw.getGhostStart(), std::array<int, 1>({ -num_ghost }));
              CHECK_EQ(slice_tw.getStart(), std::array<int, 1>({ 0 }));
              CHECK_EQ(slice_tw.getEnd(), std::array<int, 1>({ ny - 1 }));
              CHECK_EQ(slice_tw.getGhostEnd(), std::array<int, 1>({ ny - 1 + num_ghost }));

              for (int yi = -num_ghost; yi < ny + num_ghost; yi++) {
                CHECK_EQ(&ld[{ -1 - xi, yi, -1 - zi }], &slice_bw[{ yi }]);
                CHECK_EQ(&ld[{ nx + xi, yi, nz + zi }], &slice_te[{ yi }]);
                CHECK_EQ(&ld[{ nx + xi, yi, -1 - zi }], &slice_be[{ yi }]);
                CHECK_EQ(&ld[{ -1 - xi, yi, nz + zi }], &slice_tw[{ yi }]);
              }
            }
          }

          for (unsigned char yi = 0; yi < num_ghost; yi++) {
            for (unsigned char xi = 0; xi < num_ghost; xi++) {
              View<double, 1> slice_sw = ld.getGhostSliceOn(Edge::sw(), { xi, yi });
              CHECK_EQ(slice_sw.getGhostStart(), std::array<int, 1>({ -num_ghost }));
              CHECK_EQ(slice_sw.getStart(), std::array<int, 1>({ 0 }));
              CHECK_EQ(slice_sw.getEnd(), std::array<int, 1>({ nz - 1 }));
              CHECK_EQ(slice_sw.getGhostEnd(), std::array<int, 1>({ nz - 1 + num_ghost }));

              View<double, 1> slice_ne = ld.getGhostSliceOn(Edge::ne(), { xi, yi });
              CHECK_EQ(slice_ne.getGhostStart(), std::array<int, 1>({ -num_ghost }));
              CHECK_EQ(slice_ne.getStart(), std::array<int, 1>({ 0 }));
              CHECK_EQ(slice_ne.getEnd(), std::array<int, 1>({ nz - 1 }));
              CHECK_EQ(slice_ne.getGhostEnd(), std::array<int, 1>({ nz - 1 + num_ghost }));

              View<double, 1> slice_se = ld.getGhostSliceOn(Edge::se(), { xi, yi });
              CHECK_EQ(slice_se.getGhostStart(), std::array<int, 1>({ -num_ghost }));
              CHECK_EQ(slice_se.getStart(), std::array<int, 1>({ 0 }));
              CHECK_EQ(slice_se.getEnd(), std::array<int, 1>({ nz - 1 }));
              CHECK_EQ(slice_se.getGhostEnd(), std::array<int, 1>({ nz - 1 + num_ghost }));

              View<double, 1> slice_nw = ld.getGhostSliceOn(Edge::nw(), { xi, yi });
              CHECK_EQ(slice_nw.getGhostStart(), std::array<int, 1>({ -num_ghost }));
              CHECK_EQ(slice_nw.getStart(), std::array<int, 1>({ 0 }));
              CHECK_EQ(slice_nw.getEnd(), std::array<int, 1>({ nz - 1 }));
              CHECK_EQ(slice_nw.getGhostEnd(), std::array<int, 1>({ nz - 1 + num_ghost }));

              for (int zi = -num_ghost; zi < nz + num_ghost; zi++) {
                CHECK_EQ(&ld[{ -1 - xi, -1 - yi, zi }], &slice_sw[{ zi }]);
                CHECK_EQ(&ld[{ nx + xi, ny + yi, zi }], &slice_ne[{ zi }]);
                CHECK_EQ(&ld[{ nx + xi, -1 - yi, zi }], &slice_se[{ zi }]);
                CHECK_EQ(&ld[{ -1 - xi, ny + yi, zi }], &slice_nw[{ zi }]);
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView<double,3> getGhostSliceOn<0>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto nz : { 2, 3 }) {
        for (auto num_ghost : { 0, 1, 2 }) {

          array<int, 3> lengths = { nx, ny, nz };
          array<int, 3> strides = { 1, nx + 2 * num_ghost, (nx + 2 * num_ghost) * (ny + 2 * num_ghost) };
          int size = 1;
          for (size_t i = 0; i < 3; i++) {
            size *= (lengths[i] + 2 * num_ghost);
          }
          double data[size];

          const ComponentView<double, 3> ld(data, strides, lengths, num_ghost);

          for (unsigned char zi = 0; zi < num_ghost; zi++) {
            for (unsigned char yi = 0; yi < num_ghost; yi++) {
              for (unsigned char xi = 0; xi < num_ghost; xi++) {
                CHECK_EQ(&ld[{ -1 - xi, -1 - yi, -1 - zi }], &(ld.getGhostSliceOn(Corner<3>::bsw(), { xi, yi, zi }))[{}]);
                CHECK_EQ(&ld[{ nx + xi, -1 - yi, -1 - zi }], &(ld.getGhostSliceOn(Corner<3>::bse(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ -1 - xi, ny + yi, -1 - zi }], &(ld.getGhostSliceOn(Corner<3>::bnw(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ nx + xi, ny + yi, -1 - zi }], &(ld.getGhostSliceOn(Corner<3>::bne(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ -1 - xi, -1 - yi, nz + zi }], &(ld.getGhostSliceOn(Corner<3>::tsw(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ nx + xi, -1 - yi, nz + zi }], &(ld.getGhostSliceOn(Corner<3>::tse(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ -1 - xi, ny + yi, nz + zi }], &(ld.getGhostSliceOn(Corner<3>::tnw(), { xi, yi, zi })[{}]));
                CHECK_EQ(&ld[{ nx + xi, ny + yi, nz + zi }], &(ld.getGhostSliceOn(Corner<3>::tne(), { xi, yi, zi })[{}]));
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView<double,2> getGhostSliceOn<0>")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        array<int, 2> lengths = { nx, ny };
        array<int, 2> strides = { 1, nx + 2 * num_ghost };
        int size = 1;
        for (size_t i = 0; i < 2; i++) {
          size *= (lengths[i] + 2 * num_ghost);
        }
        double data[size];

        const ComponentView<double, 2> ld(data, strides, lengths, num_ghost);

        for (unsigned char yi = 0; yi < num_ghost; yi++) {
          for (unsigned char xi = 0; xi < num_ghost; xi++) {
            CHECK_EQ(&ld[{ -xi - 1, -yi - 1 }], &(ld.getGhostSliceOn(Corner<2>::sw(), { xi, yi })[{}]));
            CHECK_EQ(&ld[{ nx + xi, -yi - 1 }], &(ld.getGhostSliceOn(Corner<2>::se(), { xi, yi })[{}]));
            CHECK_EQ(&ld[{ -xi - 1, ny + yi }], &(ld.getGhostSliceOn(Corner<2>::nw(), { xi, yi })[{}]));
            CHECK_EQ(&ld[{ nx + xi, ny + yi }], &(ld.getGhostSliceOn(Corner<2>::ne(), { xi, yi })[{}]));
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView squarebracket operator")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        array<int, 2> lengths = { nx, ny };
        array<int, 2> strides = { 1, nx + 2 * num_ghost };
        int size = 1;
        for (size_t i = 0; i < 2; i++) {
          size *= (lengths[i] + 2 * num_ghost);
        }
        double data[size];

        ComponentView<double, 2> v(data, strides, lengths, num_ghost);

        double* start = &v[{ 0, 0 }];

        for (int yi = -num_ghost - 1; yi < ny + num_ghost + 1; yi++) {
          for (int xi = -num_ghost - 1; xi < nx + num_ghost + 1; xi++) {
            if (xi < -num_ghost || xi >= nx + num_ghost || yi < -num_ghost || yi >= ny + num_ghost) {
              // oob coord
              if constexpr (ENABLE_DEBUG) {
                CHECK_THROWS_AS((v[{ xi, yi }]), RuntimeError);
              }
            } else {
              CHECK_EQ(&v[{ xi, yi }], start + xi + yi * (nx + 2 * num_ghost));
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView squarebracket operator const")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        array<int, 2> lengths = { nx, ny };
        array<int, 2> strides = { 1, nx + 2 * num_ghost };
        int size = 1;
        for (size_t i = 0; i < 2; i++) {
          size *= (lengths[i] + 2 * num_ghost);
        }
        double data[size];

        const ComponentView<double, 2> v(data, strides, lengths, num_ghost);

        const double* start = &v[{ 0, 0 }];
        for (int yi = -num_ghost - 1; yi < ny + num_ghost + 1; yi++) {
          for (int xi = -num_ghost - 1; xi < nx + num_ghost + 1; xi++) {
            if (xi < -num_ghost || xi >= nx + num_ghost || yi < -num_ghost || yi >= ny + num_ghost) {
              // oob coord
              if constexpr (ENABLE_DEBUG) {
                CHECK_THROWS_AS((v[{ xi, yi }]), RuntimeError);
              }
            } else {
              CHECK_EQ(&v[{ xi, yi }], start + xi + yi * (nx + 2 * num_ghost));
            }
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView set")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        array<int, 2> lengths = { nx, ny };
        array<int, 2> strides = { 1, nx + 2 * num_ghost };
        int size = 1;
        for (size_t i = 0; i < 2; i++) {
          size *= (lengths[i] + 2 * num_ghost);
        }
        double data[size];

        ComponentView<double, 2> v(data, strides, lengths, num_ghost);

        double value = 0;
        for (int yi = -num_ghost - 1; yi < ny + num_ghost + 1; yi++) {
          for (int xi = -num_ghost - 1; xi < nx + num_ghost + 1; xi++) {
            if ((xi < -num_ghost) || (xi >= nx + num_ghost) || (yi < -num_ghost) || (yi >= ny + num_ghost)) {
              // oob coord
              if constexpr (ENABLE_DEBUG) {
                CHECK_THROWS_AS(v.set({ xi, yi }, value), RuntimeError);
              }
            } else {
              v.set({ xi, yi }, value);
              CHECK_EQ(v[{ xi, yi }], value);
            }
            value++;
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView set const")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        array<int, 2> lengths = { nx, ny };
        array<int, 2> strides = { 1, nx + 2 * num_ghost };
        int size = 1;
        for (size_t i = 0; i < 2; i++) {
          size *= (lengths[i] + 2 * num_ghost);
        }
        double data[size];

        ComponentView<const double, 2> v(data, strides, lengths, num_ghost);

        double value = 0;
        for (int yi = -num_ghost - 1; yi < ny + num_ghost + 1; yi++) {
          for (int xi = -num_ghost - 1; xi < nx + num_ghost + 1; xi++) {
            if (xi >= 0 && xi < nx && yi >= 0 && yi < ny) {
              // intertior coord
              if constexpr (ENABLE_DEBUG) {
                CHECK_THROWS_AS(v.set({ xi, yi }, value), RuntimeError);
              }
            } else if (xi < -num_ghost || xi >= nx + num_ghost || yi < -num_ghost || yi >= ny + num_ghost) {
              // oob coord
              if constexpr (ENABLE_DEBUG) {
                CHECK_THROWS_AS(v.set({ xi, yi }, value), RuntimeError);
              }
            } else {
              v.set({ xi, yi }, value);
              CHECK_EQ(v[{ xi, yi }], value);
            }
            value++;
          }
        }
      }
    }
  }
}
TEST_CASE("ComponentView implicit conversion to const type")
{
  for (auto nx : { 2, 3 }) {
    for (auto ny : { 2, 3 }) {
      for (auto num_ghost : { 0, 1, 2 }) {

        array<int, 2> lengths = { nx, ny };
        array<int, 2> strides = { 1, nx + 2 * num_ghost };
        int size = 1;
        for (size_t i = 0; i < 2; i++) {
          size *= (lengths[i] + 2 * num_ghost);
        }
        double data[size];

        ComponentView<double, 2> v(data, strides, lengths, num_ghost);
        ComponentView<const double, 2> vc = v;

        CHECK_EQ(vc.getGhostStart(), v.getGhostStart());
        CHECK_EQ(vc.getStart(), v.getStart());
        CHECK_EQ(vc.getEnd(), v.getEnd());
        CHECK_EQ(vc.getGhostEnd(), v.getGhostEnd());
        CHECK_EQ(vc.getStrides(), v.getStrides());
        CHECK_EQ(&vc[vc.getGhostStart()], &v[v.getGhostStart()]);
      }
    }
  }
}
