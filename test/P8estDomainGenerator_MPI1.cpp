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

#include "P8estDomainGenerator_SHARED.h"
#include <ThunderEgg/P8estDomainGenerator.h>

#include <ThunderEgg/Domain.h>
#include <mpi.h>
#include <p4est_base.h>
#include <p8est.h>
#include <p8est_connectivity.h>
#include <p8est_extended.h>
#include <p8est_geometry.h>
#include <vector>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

namespace {
struct FourTreeBSW
{
  p8est_connectivity_t* conn;
  p8est_geometry_t* geom;
  p8est_t* p8est;
  P8estDomainGenerator::BlockMapFunc bmf;
  double scale_x = 1.0;
  double scale_y = 1.0;
  double scale_z = 1.0;
  int n;
  int rank;

  FourTreeBSW(int base_level)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    n = 1 << base_level;
    conn = p8est_connectivity_new_brick(n, n, n, 0, 0, 0);

    geom = p8est_geometry_new_connectivity(conn);

    p8est = p8est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

    for (int i = base_level; i < 2; i++) {
      p8est_refine(
        p8est,
        false,
        [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int { return 1; },
        nullptr);
    }

    p8est_partition(p8est, true, nullptr);

    bmf = [&](int block_no, double unit_x, double unit_y, double unit_z, double& x, double& y, double& z) {
      const double coord[3] = { unit_x, unit_y, unit_z };
      double out_coord[3];
      p8est_geometry_connectivity_X(geom, block_no, coord, out_coord);
      x = scale_x * out_coord[0] / n;
      y = scale_y * out_coord[1] / n;
      z = scale_z * out_coord[2] / n;
    };
  }
  ~FourTreeBSW()
  {
    p8est_destroy(p8est);
    p8est_geometry_destroy(geom);
    p8est_connectivity_destroy(conn);
  }
};
} // namespace
TEST_CASE("P8estDomainGenerator 4x4x4 hasCoarserDomain")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    FourTreeBSW tree(base_level);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    for (int i = base_level; i < 3; i++) {
      CHECK_UNARY(dg.hasCoarserDomain());
      Domain<3> domain = dg.getCoarserDomain();
    }
    CHECK_UNARY_FALSE(dg.hasCoarserDomain());
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 Uniform Number of Patches")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    FourTreeBSW tree(base_level);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    for (int curr_level = 2; curr_level >= base_level; curr_level--) {
      Domain<3> domain = dg.getCoarserDomain();
      int n = 1 << curr_level; // 2^curr_level
      CHECK_EQ(domain.getNumGlobalPatches(), n * n * n);
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 RefineLevel")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    FourTreeBSW tree(base_level);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    for (int curr_level = 2 - base_level; curr_level >= 0; curr_level--) {
      Domain<3> domain = dg.getCoarserDomain();
      for (auto patch : domain.getPatchInfoVector()) {
        CHECK_EQ(patch.refine_level, curr_level);
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 rank")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    FourTreeBSW tree(base_level);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    for (int curr_level = 2 - base_level; curr_level >= 0; curr_level--) {
      Domain<3> domain = dg.getCoarserDomain();
      for (auto patch : domain.getPatchInfoVector()) {
        CHECK_EQ(patch.rank, tree.rank);
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 spacings")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    for (int nx : { 5, 10 }) {
      for (int ny : { 5, 10 }) {
        for (int nz : { 5, 10 }) {
          for (double scale_x : { 0.5, 1.0 }) {
            for (double scale_y : { 0.5, 1.0 }) {
              for (double scale_z : { 0.5, 1.0 }) {
                FourTreeBSW tree(base_level);

                tree.scale_x = scale_x;
                tree.scale_y = scale_y;
                tree.scale_z = scale_z;

                P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, tree.bmf);

                for (int curr_level = 2; curr_level >= base_level; curr_level--) {
                  int n = 1 << curr_level;
                  double patch_length = 1.0 / n;
                  Domain<3> domain = dg.getCoarserDomain();
                  for (auto patch : domain.getPatchInfoVector()) {
                    CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * patch_length / nx));
                    CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * patch_length / ny));
                    CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * patch_length / nz));
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
TEST_CASE("P8estDomainGenerator 4x4x4 ns")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    for (int nx : { 5, 10 }) {
      for (int ny : { 5, 10 }) {
        for (int nz : { 5, 10 }) {
          FourTreeBSW tree(base_level);

          P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, tree.bmf);

          for (int curr_level = 2 - base_level; curr_level >= 0; curr_level--) {
            Domain<3> domain = dg.getCoarserDomain();
            for (auto patch : domain.getPatchInfoVector()) {
              CHECK_EQ(patch.ns[0], nx);
              CHECK_EQ(patch.ns[1], ny);
              CHECK_EQ(patch.ns[2], nz);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 starts")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    for (int nx : { 5, 10 }) {
      for (int ny : { 5, 10 }) {
        for (int nz : { 5, 10 }) {
          FourTreeBSW tree(base_level);

          P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, tree.bmf);

          for (int curr_level = 2; curr_level >= base_level; curr_level--) {
            Domain<3> domain = dg.getCoarserDomain();

            int n = 1 << curr_level; // 2^curr_level
            vector<vector<vector<int>>> num_patches(n, vector<vector<int>>(n, vector<int>(4, 0)));
            for (auto patch : domain.getPatchInfoVector()) {
              int i = (int)(patch.starts[0] * n);
              int j = (int)(patch.starts[1] * n);
              int k = (int)(patch.starts[2] * n);
              CHECK_EQ(patch.starts[0] * n, doctest::Approx((double)i));
              CHECK_EQ(patch.starts[1] * n, doctest::Approx((double)j));
              CHECK_EQ(patch.starts[2] * n, doctest::Approx((double)k));
              num_patches[i][j][k]++;
            }
            // Check that each patch is unique
            for (int i = 0; i < n; i++) {
              for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                  CHECK_EQ(num_patches[i][j][k], 1);
                }
              }
            }
          }
        }
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 num_ghost_cells")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    for (int num_ghost_cells : { 0, 1, 2 }) {
      FourTreeBSW tree(base_level);

      P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, num_ghost_cells, tree.bmf);
      for (int curr_level = 2 - base_level; curr_level >= 0; curr_level--) {
        Domain<3> domain = dg.getCoarserDomain();
        for (auto patch : domain.getPatchInfoVector()) {
          CHECK_EQ(patch.num_ghost_cells, num_ghost_cells);
        }
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 Uniform neighbor nfos")
{
  {
    FourTreeBSW tree(0);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> domain_2 = dg.getCoarserDomain();
    Domain<3> domain_1 = dg.getCoarserDomain();
    Domain<3> domain_0 = dg.getCoarserDomain();

    CheckRootDomainNeighbors(domain_0);

    Check2x2x2DomainNeighbors(domain_1);

    Check4x4x4DomainNeighbors(domain_2);
  }
  {
    FourTreeBSW tree(1);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> domain_2 = dg.getCoarserDomain();
    Domain<3> domain_1 = dg.getCoarserDomain();

    Check2x2x2DomainNeighbors(domain_1);

    Check4x4x4DomainNeighbors(domain_2);
  }
  {
    FourTreeBSW tree(2);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> domain_2 = dg.getCoarserDomain();

    Check4x4x4DomainNeighbors(domain_2);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4 Uniform child/parent ids and ranks")
{
  {
    FourTreeBSW tree(0);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> domain_2 = dg.getCoarserDomain();
    Domain<3> domain_1 = dg.getCoarserDomain();
    Domain<3> domain_0 = dg.getCoarserDomain();

    CheckChildIdsAndRanksNull(domain_2);

    CheckParentAndChildIdsAndRanks(domain_1, 1, domain_2, 2);

    CheckParentAndChildIdsAndRanks(domain_0, 0, domain_1, 1);

    CheckParentIdsAndRanksNull(domain_0);
  }
  {
    FourTreeBSW tree(1);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> domain_2 = dg.getCoarserDomain();
    Domain<3> domain_1 = dg.getCoarserDomain();

    CheckChildIdsAndRanksNull(domain_2);

    CheckParentAndChildIdsAndRanks(domain_1, 1, domain_2, 2);

    CheckParentIdsAndRanksNull(domain_1);
  }
  {
    FourTreeBSW tree(2);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> domain_2 = dg.getCoarserDomain();

    CheckChildIdsAndRanksNull(domain_2);

    CheckParentIdsAndRanksNull(domain_2);
  }
}
namespace {
struct FourTreeRefineBSW
{
  p8est_connectivity_t* conn;
  p8est_geometry_t* geom;
  p8est_t* p8est;
  P8estDomainGenerator::BlockMapFunc bmf;
  double scale_x = 1.0;
  double scale_y = 1.0;
  double scale_z = 1.0;
  int n;
  int rank;

  FourTreeRefineBSW(int base_level)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    conn = p8est_connectivity_new_unitcube();

    p8est = p8est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    n = 1 << base_level;
    conn = p8est_connectivity_new_brick(n, n, n, 0, 0, 0);

    geom = p8est_geometry_new_connectivity(conn);

    p8est = p8est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

    for (int i = base_level; i < 2; i++) {
      p8est_refine(
        p8est,
        false,
        [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int { return 1; },
        nullptr);
    }
    p8est_refine(
      p8est,
      false,
      [](p8est_t* p8est, p4est_topidx_t witch_tree, p8est_quadrant_t* quadrant) -> int {
        return witch_tree == 0 && quadrant->x == 0 && quadrant->y == 0 && quadrant->z == 0;
      },
      nullptr);

    p8est_partition(p8est, true, nullptr);

    bmf = [&](int block_no, double unit_x, double unit_y, double unit_z, double& x, double& y, double& z) {
      const double coord[3] = { unit_x, unit_y, unit_z };
      double out_coord[3];
      p8est_geometry_connectivity_X(geom, block_no, coord, out_coord);
      x = scale_x * out_coord[0] / n;
      y = scale_y * out_coord[1] / n;
      z = scale_z * out_coord[2] / n;
    };
  }
  ~FourTreeRefineBSW()
  {
    p8est_destroy(p8est);
    p8est_geometry_destroy(geom);
    p8est_connectivity_destroy(conn);
  }
};
} // namespace
TEST_CASE("P8estDomainGenerator 4x4x4rbsw hasCoarserDomain")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    FourTreeRefineBSW tree(base_level);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    for (int i = base_level; i < 3 + 1; i++) {
      CHECK_UNARY(dg.hasCoarserDomain());
      Domain<3> domain = dg.getCoarserDomain();
    }
    CHECK_UNARY_FALSE(dg.hasCoarserDomain());
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw Uniform Number of Patches")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    FourTreeRefineBSW tree(base_level);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> finest_domain = dg.getCoarserDomain();
    CHECK_EQ(finest_domain.getNumGlobalPatches(), 71);

    for (int curr_level = 2; curr_level >= base_level; curr_level--) {
      Domain<3> domain = dg.getCoarserDomain();
      int n = 1 << curr_level; // 2^curr_level
      CHECK_EQ(domain.getNumGlobalPatches(), n * n * n);
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw RefineLevel")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    FourTreeRefineBSW tree(base_level);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> finest_domain = dg.getCoarserDomain();
    for (auto patch : finest_domain.getPatchInfoVector()) {
      if (patch.starts[0] < 0.24 && patch.starts[1] < 0.24 && patch.starts[2] < 0.24) {
        CHECK_EQ(patch.refine_level, 3 - base_level);
      } else {
        CHECK_EQ(patch.refine_level, 2 - base_level);
      }
    }

    for (int curr_level = 2 - base_level; curr_level >= 0; curr_level--) {
      Domain<3> domain = dg.getCoarserDomain();
      for (auto patch : domain.getPatchInfoVector()) {
        CHECK_EQ(patch.refine_level, curr_level);
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw rank")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    FourTreeRefineBSW tree(base_level);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    for (int curr_level = 2 - base_level + 1; curr_level >= 0; curr_level--) {
      Domain<3> domain = dg.getCoarserDomain();
      for (auto patch : domain.getPatchInfoVector()) {
        CHECK_EQ(patch.rank, tree.rank);
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw spacings")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    for (int nx : { 5, 10 }) {
      for (int ny : { 5, 10 }) {
        for (int nz : { 5, 10 }) {
          for (double scale_x : { 0.5, 1.0 }) {
            for (double scale_y : { 0.5, 1.0 }) {
              for (double scale_z : { 0.5, 1.0 }) {
                FourTreeRefineBSW tree(base_level);

                tree.scale_x = scale_x;
                tree.scale_y = scale_y;
                tree.scale_z = scale_z;

                P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, tree.bmf);

                Domain<3> finest_domain = dg.getCoarserDomain();

                for (auto patch : finest_domain.getPatchInfoVector()) {
                  if (patch.starts[0] < 0.24 * scale_x && patch.starts[1] < 0.24 * scale_y &&
                      patch.starts[2] < 0.24 * scale_z) {
                    CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 0.125 / nx));
                    CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 0.125 / ny));
                    CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 0.125 / nz));
                  } else {
                    CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * 0.25 / nx));
                    CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * 0.25 / ny));
                    CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * 0.25 / nz));
                  }
                }

                for (int curr_level = 2; curr_level >= base_level; curr_level--) {
                  int n = 1 << curr_level;
                  double patch_length = 1.0 / n;
                  Domain<3> domain = dg.getCoarserDomain();
                  for (auto patch : domain.getPatchInfoVector()) {
                    CHECK_EQ(patch.spacings[0], doctest::Approx(scale_x * patch_length / nx));
                    CHECK_EQ(patch.spacings[1], doctest::Approx(scale_y * patch_length / ny));
                    CHECK_EQ(patch.spacings[2], doctest::Approx(scale_z * patch_length / nz));
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
TEST_CASE("P8estDomainGenerator 4x4x4rbsw ns")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    for (int nx : { 5, 10 }) {
      for (int ny : { 5, 10 }) {
        for (int nz : { 5, 10 }) {
          FourTreeRefineBSW tree(base_level);

          P8estDomainGenerator dg(tree.p8est, { nx, ny, nz }, 1, tree.bmf);

          for (int curr_level = 3 - base_level; curr_level >= 0; curr_level--) {
            Domain<3> domain = dg.getCoarserDomain();
            for (auto patch : domain.getPatchInfoVector()) {
              CHECK_EQ(patch.ns[0], nx);
              CHECK_EQ(patch.ns[1], ny);
              CHECK_EQ(patch.ns[2], nz);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw num_ghost_cells")
{
  for (int base_level = 0; base_level < 3; base_level++) {
    for (int num_ghost_cells : { 0, 1, 2 }) {
      FourTreeRefineBSW tree(base_level);

      P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, num_ghost_cells, tree.bmf);

      for (int curr_level = 3 - base_level; curr_level >= 0; curr_level--) {
        Domain<3> domain = dg.getCoarserDomain();
        for (auto patch : domain.getPatchInfoVector()) {
          CHECK_EQ(patch.num_ghost_cells, num_ghost_cells);
        }
      }
    }
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw neighbor nfos")
{
  {
    FourTreeRefineBSW tree(0);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> domain_3 = dg.getCoarserDomain();
    Domain<3> domain_2 = dg.getCoarserDomain();
    Domain<3> domain_1 = dg.getCoarserDomain();
    Domain<3> domain_0 = dg.getCoarserDomain();

    CheckRootDomainNeighbors(domain_0);

    Check2x2x2DomainNeighbors(domain_1);

    Check4x4x4DomainNeighbors(domain_2);

    Check4x4x4RefinedBSWDomainNeighbors(domain_3);
  }
  {
    FourTreeRefineBSW tree(1);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> domain_3 = dg.getCoarserDomain();
    Domain<3> domain_2 = dg.getCoarserDomain();
    Domain<3> domain_1 = dg.getCoarserDomain();

    Check2x2x2DomainNeighbors(domain_1);

    Check4x4x4DomainNeighbors(domain_2);

    Check4x4x4RefinedBSWDomainNeighbors(domain_3);
  }
  {
    FourTreeRefineBSW tree(2);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> domain_3 = dg.getCoarserDomain();
    Domain<3> domain_2 = dg.getCoarserDomain();

    Check4x4x4DomainNeighbors(domain_2);

    Check4x4x4RefinedBSWDomainNeighbors(domain_3);
  }
}
TEST_CASE("P8estDomainGenerator 4x4x4rbsw child/parent ids and ranks")
{
  {
    FourTreeRefineBSW tree(0);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> domain_3 = dg.getCoarserDomain();
    Domain<3> domain_2 = dg.getCoarserDomain();
    Domain<3> domain_1 = dg.getCoarserDomain();
    Domain<3> domain_0 = dg.getCoarserDomain();

    CheckChildIdsAndRanksNull(domain_3);

    CheckParentAndChildIdsAndRanksRefined(domain_2, 2, domain_3, 3);

    CheckParentAndChildIdsAndRanks(domain_1, 1, domain_2, 2);

    CheckParentAndChildIdsAndRanks(domain_0, 0, domain_1, 1);

    CheckParentIdsAndRanksNull(domain_0);
  }
  {
    FourTreeRefineBSW tree(1);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> domain_3 = dg.getCoarserDomain();
    Domain<3> domain_2 = dg.getCoarserDomain();
    Domain<3> domain_1 = dg.getCoarserDomain();

    CheckChildIdsAndRanksNull(domain_3);

    CheckParentAndChildIdsAndRanksRefined(domain_2, 2, domain_3, 3);

    CheckParentAndChildIdsAndRanks(domain_1, 1, domain_2, 2);

    CheckParentIdsAndRanksNull(domain_1);
  }
  {
    FourTreeRefineBSW tree(2);

    P8estDomainGenerator dg(tree.p8est, { 10, 10, 10 }, 1, tree.bmf);

    Domain<3> domain_3 = dg.getCoarserDomain();
    Domain<3> domain_2 = dg.getCoarserDomain();

    CheckChildIdsAndRanksNull(domain_3);

    CheckParentAndChildIdsAndRanksRefined(domain_2, 2, domain_3, 3);

    CheckParentIdsAndRanksNull(domain_2);
  }
}