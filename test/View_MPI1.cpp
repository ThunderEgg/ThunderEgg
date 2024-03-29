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
#include <ThunderEgg/View.h>
#include <doctest.h>

using namespace std;
using namespace ThunderEgg;
TEST_CASE("View default constructor")
{
  View<double, 2> v;

  for (int i = 0; i < 2; i++) {
    CHECK_EQ(v.getStrides()[i], 0);
    CHECK_EQ(v.getStart()[i], 0);
    CHECK_EQ(v.getEnd()[i], -1);
    CHECK_EQ(v.getGhostStart()[i], 0);
    CHECK_EQ(v.getGhostEnd()[i], -1);
  }
  if constexpr (ENABLE_DEBUG) {
    CHECK_THROWS_AS((v[{ 0, 0 }]), RuntimeError);
    CHECK_THROWS_AS(v(0, 0), RuntimeError);
  }
}
TEST_CASE("View constructor")
{
  for (auto x_stride : { 1, 2 }) {
    for (auto y_stride : { 1, 2 }) {
      for (auto x_ghost_start : { -1, 1 }) {
        for (auto x_start : { 2, 3 }) {
          for (auto x_end : { 3, 4 }) {
            for (auto x_ghost_end : { 5, 6 }) {
              for (auto y_ghost_start : { -1, 1 }) {
                for (auto y_start : { 2, 3 }) {
                  for (auto y_end : { 3, 4 }) {
                    for (auto y_ghost_end : { 5, 6 }) {

                      double data;

                      View<double, 2> v(&data, { x_stride, y_stride }, { x_ghost_start, y_ghost_start }, { x_start, y_start }, { x_end, y_end }, { x_ghost_end, y_ghost_end });

                      CHECK_EQ(v.getStrides()[0], x_stride);
                      CHECK_EQ(v.getStrides()[1], y_stride);
                      CHECK_EQ(v.getGhostStart()[0], x_ghost_start);
                      CHECK_EQ(v.getGhostStart()[1], y_ghost_start);
                      CHECK_EQ(v.getStart()[0], x_start);
                      CHECK_EQ(v.getStart()[1], y_start);
                      CHECK_EQ(v.getEnd()[0], x_end);
                      CHECK_EQ(v.getEnd()[1], y_end);
                      CHECK_EQ(v.getGhostEnd()[0], x_ghost_end);
                      CHECK_EQ(v.getGhostEnd()[1], y_ghost_end);
                      CHECK_EQ(&v[v.getGhostStart()], &data);
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
}
TEST_CASE("View squarebracket operator")
{
  for (auto x_stride : { 1, 2 }) {
    for (auto y_stride : { 4, 5 }) {
      for (auto x_ghost_start : { -1, 1 }) {
        for (auto x_start : { 2, 3 }) {
          for (auto x_end : { 3, 4 }) {
            for (auto x_ghost_end : { 5, 6 }) {
              for (auto y_ghost_start : { -1, 1 }) {
                for (auto y_start : { 2, 3 }) {
                  for (auto y_end : { 3, 4 }) {
                    for (auto y_ghost_end : { 5, 6 }) {

                      double data[100];
                      for (int i = 0; i < 100; i++) {
                        data[i] = 0;
                      }

                      View<double, 2> v(data, { x_stride, y_stride }, { x_ghost_start, y_ghost_start }, { x_start, y_start }, { x_end, y_end }, { x_ghost_end, y_ghost_end });

                      double* origin = data - (x_ghost_start * x_stride + y_ghost_start * y_stride);
                      for (int yi = y_ghost_start - 1; yi <= y_ghost_end + 1; yi++) {
                        for (int xi = x_ghost_start - 1; xi <= x_ghost_end + 1; xi++) {
                          if (xi < x_ghost_start || xi > x_ghost_end || yi < y_ghost_start || yi > y_ghost_end) {
                            // oob coord
                            if constexpr (ENABLE_DEBUG) {
                              CHECK_THROWS_AS((v[{ xi, yi }]), RuntimeError);
                            }
                          } else {
                            CHECK_EQ(&v[{ xi, yi }], origin + xi * x_stride + yi * y_stride);
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
      }
    }
  }
}
TEST_CASE("View parens operator")
{
  for (auto x_stride : { 1, 2 }) {
    for (auto y_stride : { 4, 5 }) {
      for (auto x_ghost_start : { -1, 1 }) {
        for (auto x_start : { 2, 3 }) {
          for (auto x_end : { 3, 4 }) {
            for (auto x_ghost_end : { 5, 6 }) {
              for (auto y_ghost_start : { -1, 1 }) {
                for (auto y_start : { 2, 3 }) {
                  for (auto y_end : { 3, 4 }) {
                    for (auto y_ghost_end : { 5, 6 }) {

                      double data[100];
                      for (int i = 0; i < 100; i++) {
                        data[i] = 0;
                      }

                      View<double, 2> v(data, { x_stride, y_stride }, { x_ghost_start, y_ghost_start }, { x_start, y_start }, { x_end, y_end }, { x_ghost_end, y_ghost_end });

                      double* origin = data - (x_ghost_start * x_stride + y_ghost_start * y_stride);
                      for (int yi = y_ghost_start - 1; yi <= y_ghost_end + 1; yi++) {
                        for (int xi = x_ghost_start - 1; xi <= x_ghost_end + 1; xi++) {
                          if (xi < x_ghost_start || xi > x_ghost_end || yi < y_ghost_start || yi > y_ghost_end) {
                            // oob coord
                            if constexpr (ENABLE_DEBUG) {
                              CHECK_THROWS_AS((v(xi, yi)), RuntimeError);
                            }
                          } else {
                            CHECK_EQ(&v(xi, yi), origin + xi * x_stride + yi * y_stride);
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
      }
    }
  }
}
TEST_CASE("View set")
{
  for (auto x_stride : { 1, 2 }) {
    for (auto y_stride : { 4, 5 }) {
      for (auto x_ghost_start : { -1, 1 }) {
        for (auto x_start : { 2, 3 }) {
          for (auto x_end : { 3, 4 }) {
            for (auto x_ghost_end : { 5, 6 }) {
              for (auto y_ghost_start : { -1, 1 }) {
                for (auto y_start : { 2, 3 }) {
                  for (auto y_end : { 3, 4 }) {
                    for (auto y_ghost_end : { 5, 6 }) {

                      double data[1000];
                      for (int i = 0; i < 1000; i++) {
                        data[i] = 0;
                      }

                      View<double, 2> v(data, { x_stride, y_stride }, { x_ghost_start, y_ghost_start }, { x_start, y_start }, { x_end, y_end }, { x_ghost_end, y_ghost_end });

                      double value = 0;
                      for (int yi = y_ghost_start - 1; yi <= y_ghost_end + 1; yi++) {
                        for (int xi = x_ghost_start - 1; xi <= x_ghost_end + 1; xi++) {
                          if (xi < x_ghost_start || xi > x_ghost_end || yi < y_ghost_start || yi > y_ghost_end) {
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
            }
          }
        }
      }
    }
  }
}
TEST_CASE("View set const")
{
  for (auto x_stride : { 1, 2 }) {
    for (auto y_stride : { 4, 5 }) {
      for (auto x_ghost_start : { -1, 1 }) {
        for (auto x_start : { 2, 3 }) {
          for (auto x_end : { 3, 4 }) {
            for (auto x_ghost_end : { 5, 6 }) {
              for (auto y_ghost_start : { -1, 1 }) {
                for (auto y_start : { 2, 3 }) {
                  for (auto y_end : { 3, 4 }) {
                    for (auto y_ghost_end : { 5, 6 }) {

                      double data[1000];
                      for (int i = 0; i < 1000; i++) {
                        data[i] = 0;
                      }

                      View<const double, 2> v(data, { x_stride, y_stride }, { x_ghost_start, y_ghost_start }, { x_start, y_start }, { x_end, y_end }, { x_ghost_end, y_ghost_end });

                      double value = 0;
                      for (int yi = y_ghost_start - 1; yi <= y_ghost_end + 1; yi++) {
                        for (int xi = x_ghost_start - 1; xi <= x_ghost_end + 1; xi++) {
                          if (xi >= x_start && xi <= x_end && yi >= y_start && yi <= y_end) {
                            // intertior coord
                            if constexpr (ENABLE_DEBUG) {
                              CHECK_THROWS_AS(v.set({ xi, yi }, value), RuntimeError);
                            }
                          } else if (xi < x_ghost_start || xi > x_ghost_end || yi < y_ghost_start || yi > y_ghost_end) {
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
            }
          }
        }
      }
    }
  }
}
TEST_CASE("View implicit conversion to const type")
{
  for (auto x_stride : { 1, 2 }) {
    for (auto y_stride : { 4, 5 }) {
      for (auto x_ghost_start : { -1, 1 }) {
        for (auto x_start : { 2, 3 }) {
          for (auto x_end : { 3, 4 }) {
            for (auto x_ghost_end : { 5, 6 }) {
              for (auto y_ghost_start : { -1, 1 }) {
                for (auto y_start : { 2, 3 }) {
                  for (auto y_end : { 3, 4 }) {
                    for (auto y_ghost_end : { 5, 6 }) {

                      double data[1000];
                      for (int i = 0; i < 1000; i++) {
                        data[i] = 0;
                      }

                      View<double, 2> v(data, { x_stride, y_stride }, { x_ghost_start, y_ghost_start }, { x_start, y_start }, { x_end, y_end }, { x_ghost_end, y_ghost_end });
                      View<const double, 2> vc = v;

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
            }
          }
        }
      }
    }
  }
}
