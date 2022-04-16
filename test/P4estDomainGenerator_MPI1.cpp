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

#include <ThunderEgg/P4estDomainGenerator.h>

#include <p4est.h>
#include <p4est_mesh.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("SinglePatch", "[p4estDomGen]")
{
  p4est_connectivity_t* conn = p4est_connectivity_new_unitsquare();
  p4est_t* p4est = p4est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);
  P4estDomainGenerator::BlockMapFunc bmf = [](int block_no, double unit_x, double unit_y, double& x, double& y) {
    x = unit_x;
    y = unit_y;
  };
  P4estDomainGenerator dg(p4est, { 10, 10 }, 1, bmf);
  Domain<2> domain = dg.getFinestDomain();
  CHECK(domain.getNumGlobalPatches() == 1);
  auto patch = domain.getPatchInfoVector()[0];
  CHECK_FALSE(patch.hasNbr(Side<2>::west()));
  CHECK_FALSE(patch.hasNbr(Side<2>::east()));
  CHECK_FALSE(patch.hasNbr(Side<2>::south()));
  CHECK_FALSE(patch.hasNbr(Side<2>::north()));
  CHECK_FALSE(patch.hasNbr(Corner<2>::sw()));
  CHECK_FALSE(patch.hasNbr(Corner<2>::se()));
  CHECK_FALSE(patch.hasNbr(Corner<2>::nw()));
  CHECK_FALSE(patch.hasNbr(Corner<2>::ne()));
  CHECK(patch.spacings[0] == Catch::Approx(1.0 / 10));
  CHECK(patch.spacings[1] == Catch::Approx(1.0 / 10));
  CHECK(patch.starts[0] == Catch::Approx(0));
  CHECK(patch.starts[1] == Catch::Approx(0));
  CHECK(patch.ns[0] == 10);
  CHECK(patch.ns[1] == 10);

  CHECK_FALSE(dg.hasCoarserDomain());
}
TEST_CASE("P4estDomainGenerator 2x2 Uniform", "[p4estDomGen]")
{
  for (int nx : { 5, 10 }) {
    for (int ny : { 5, 10 }) {
      for (double scale_x : { 0.5, 1.0 }) {
        for (double scale_y : { 0.5, 1.0 }) {
          for (int num_ghost_cells : { 0, 1, 2 }) {
            p4est_connectivity_t* conn = p4est_connectivity_new_unitsquare();

            p4est_t* p4est = p4est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

            p4est_refine(
              p4est, false, [](p4est_t* p4est, p4est_topidx_t witch_tree, p4est_quadrant_t* quadrant) -> int { return 1; }, nullptr);

            INFO("nx: " << nx);
            INFO("ny: " << ny);
            INFO("scale_x: " << scale_x);
            INFO("scale_x: " << scale_y);
            INFO("num_ghost_cells: " << num_ghost_cells);

            P4estDomainGenerator::BlockMapFunc bmf = [&](int block_no, double unit_x, double unit_y, double& x, double& y) {
              x = scale_x * unit_x;
              y = scale_y * unit_y;
            };

            P4estDomainGenerator dg(p4est, { nx, ny }, num_ghost_cells, bmf);

            CHECK(dg.hasCoarserDomain());
            Domain<2> domain_1 = dg.getFinestDomain();
            CHECK(dg.hasCoarserDomain());
            Domain<2> domain_0 = dg.getCoarserDomain();
            CHECK_FALSE(dg.hasCoarserDomain());

            // SECTION("correct number of patches")
            {
              CHECK(domain_1.getNumGlobalPatches() == 4);
              CHECK(domain_0.getNumGlobalPatches() == 1);
            }
            // SECTION("patches have correct spacings")
            {
              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
                CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
              }

              auto patch = domain_0.getPatchInfoVector()[0];

              CHECK(patch.spacings[0] == Catch::Approx(scale_x * 1.0 / nx));
              CHECK(patch.spacings[1] == Catch::Approx(scale_y * 1.0 / ny));
            }
            // SECTION("patches have correct ns")
            {
              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK(patch.ns[0] == nx);
                CHECK(patch.ns[1] == ny);
              }

              auto patch = domain_0.getPatchInfoVector()[0];

              CHECK(patch.ns[0] == nx);
              CHECK(patch.ns[1] == ny);
            }
            // SECTION("patches have ranks set")
            {
              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK(patch.rank == 0);
              }

              auto patch = domain_0.getPatchInfoVector()[0];

              CHECK(patch.rank == 0);
            }
            // SECTION("patches have refine_level set")
            {
              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK(patch.refine_level == 1);
              }

              auto patch = domain_0.getPatchInfoVector()[0];

              CHECK(patch.refine_level == 0);
            }
            // SECTION("patches have num_ghost_cells set")
            {
              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK(patch.num_ghost_cells == num_ghost_cells);
              }

              auto patch = domain_0.getPatchInfoVector()[0];

              CHECK(patch.num_ghost_cells == num_ghost_cells);
            }
            const PatchInfo<2>* domain_1_sw_patch = nullptr;
            const PatchInfo<2>* domain_1_se_patch = nullptr;
            const PatchInfo<2>* domain_1_nw_patch = nullptr;
            const PatchInfo<2>* domain_1_ne_patch = nullptr;
            const PatchInfo<2>* domain_0_coarser_patch = nullptr;

            for (const PatchInfo<2>& patch : domain_1.getPatchInfoVector()) {
              double x = patch.starts[0];
              double y = patch.starts[1];
              if (x == Catch::Approx(0) && y == Catch::Approx(0)) {
                domain_1_sw_patch = &patch;
              }
              if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0)) {
                domain_1_se_patch = &patch;
              }
              if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y)) {
                domain_1_nw_patch = &patch;
              }
              if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
                domain_1_ne_patch = &patch;
              }
            }

            domain_0_coarser_patch = &domain_0.getPatchInfoVector()[0];

            // SECTION("patches have starts set")
            {
              REQUIRE(domain_1_sw_patch != nullptr);
              REQUIRE(domain_1_se_patch != nullptr);
              REQUIRE(domain_1_nw_patch != nullptr);
              REQUIRE(domain_1_ne_patch != nullptr);

              CHECK(domain_0_coarser_patch->starts[0] == Catch::Approx(0.0));
              CHECK(domain_0_coarser_patch->starts[1] == Catch::Approx(0.0));
            }
            // SECTION("parent ids are set correctly")
            {
              CHECK(domain_1_sw_patch->parent_id == domain_0_coarser_patch->id);
              CHECK(domain_1_se_patch->parent_id == domain_0_coarser_patch->id);
              CHECK(domain_1_nw_patch->parent_id == domain_0_coarser_patch->id);
              CHECK(domain_1_ne_patch->parent_id == domain_0_coarser_patch->id);

              CHECK(domain_0_coarser_patch->parent_id == -1);
            }
            // SECTION("child ids are set correctly")
            {
              for (int i = 0; i < 4; i++) {
                CHECK(domain_1_sw_patch->child_ids[i] == -1);
                CHECK(domain_1_se_patch->child_ids[i] == -1);
                CHECK(domain_1_nw_patch->child_ids[i] == -1);
                CHECK(domain_1_ne_patch->child_ids[i] == -1);
              }

              CHECK(domain_0_coarser_patch->child_ids[0] == domain_1_sw_patch->id);
              CHECK(domain_0_coarser_patch->child_ids[1] == domain_1_se_patch->id);
              CHECK(domain_0_coarser_patch->child_ids[2] == domain_1_nw_patch->id);
              CHECK(domain_0_coarser_patch->child_ids[3] == domain_1_ne_patch->id);
            }
            // SECTION("orth on parent is set correctly")
            {
              CHECK(domain_1_sw_patch->orth_on_parent == Orthant<2>::sw());
              CHECK(domain_1_se_patch->orth_on_parent == Orthant<2>::se());
              CHECK(domain_1_nw_patch->orth_on_parent == Orthant<2>::nw());
              CHECK(domain_1_ne_patch->orth_on_parent == Orthant<2>::ne());

              CHECK(domain_0_coarser_patch->orth_on_parent == Orthant<2>::null());
            }
            // SECTION("parent ranks are set correctly")
            {
              CHECK(domain_1_sw_patch->parent_rank == domain_0_coarser_patch->rank);
              CHECK(domain_1_se_patch->parent_rank == domain_0_coarser_patch->rank);
              CHECK(domain_1_nw_patch->parent_rank == domain_0_coarser_patch->rank);
              CHECK(domain_1_ne_patch->parent_rank == domain_0_coarser_patch->rank);

              CHECK(domain_0_coarser_patch->parent_rank == -1);
            }
            // SECTION("child ranks are set correctly")
            {
              for (int i = 0; i < 4; i++) {
                CHECK(domain_1_sw_patch->child_ranks[i] == -1);
                CHECK(domain_1_se_patch->child_ranks[i] == -1);
                CHECK(domain_1_nw_patch->child_ranks[i] == -1);
                CHECK(domain_1_ne_patch->child_ranks[i] == -1);
              }

              CHECK(domain_0_coarser_patch->child_ranks[0] == domain_1_sw_patch->rank);
              CHECK(domain_0_coarser_patch->child_ranks[1] == domain_1_se_patch->rank);
              CHECK(domain_0_coarser_patch->child_ranks[2] == domain_1_nw_patch->rank);
              CHECK(domain_0_coarser_patch->child_ranks[3] == domain_1_ne_patch->rank);
            }
            // SECTION("nbr_info ids are correct")
            {
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Side<2>::west()));
              CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_1_se_patch->id);
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Side<2>::south()));
              CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_1_nw_patch->id);

              CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_1_sw_patch->id);
              CHECK_FALSE(domain_1_se_patch->hasNbr(Side<2>::east()));
              CHECK_FALSE(domain_1_se_patch->hasNbr(Side<2>::south()));
              CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_1_ne_patch->id);

              CHECK_FALSE(domain_1_nw_patch->hasNbr(Side<2>::west()));
              CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_1_ne_patch->id);
              CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_1_sw_patch->id);
              CHECK_FALSE(domain_1_nw_patch->hasNbr(Side<2>::north()));

              CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_1_nw_patch->id);
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Side<2>::east()));
              CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_1_se_patch->id);
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Side<2>::north()));

              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::west()));
              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::east()));
              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::south()));
              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::north()));
            }
            // SECTION("nbr_info ranks are correct")
            {
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Side<2>::west()));
              CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_se_patch->rank);
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Side<2>::south()));
              CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_nw_patch->rank);

              CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_sw_patch->rank);
              CHECK_FALSE(domain_1_se_patch->hasNbr(Side<2>::east()));
              CHECK_FALSE(domain_1_se_patch->hasNbr(Side<2>::south()));
              CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_ne_patch->rank);

              CHECK_FALSE(domain_1_nw_patch->hasNbr(Side<2>::west()));
              CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_ne_patch->rank);
              CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_sw_patch->rank);
              CHECK_FALSE(domain_1_nw_patch->hasNbr(Side<2>::north()));

              CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_nw_patch->rank);
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Side<2>::east()));
              CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_se_patch->rank);
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Side<2>::north()));

              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::west()));
              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::east()));
              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::south()));
              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Side<2>::north()));
            }
            // SECTION("corner_nbr_info ids are correct")
            {
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::sw()));
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::se()));
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_1_sw_patch->hasNbr(Corner<2>::ne()));
              CHECK(domain_1_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_1_ne_patch->id);

              CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::sw()));
              CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_1_se_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_1_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_1_nw_patch->id);
              CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::ne()));

              CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_1_nw_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_1_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_1_se_patch->id);
              CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::nw()));
              CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::ne()));

              CHECK(domain_1_ne_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_1_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_1_sw_patch->id);
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::se()));
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::nw()));
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::ne()));

              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::sw()));
              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::se()));
              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::nw()));
              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::ne()));
            }
            // SECTION("corner_nbr_info ranks are correct")
            {
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::sw()));
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::se()));
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_1_sw_patch->hasNbr(Corner<2>::ne()));
              CHECK(domain_1_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_1_ne_patch->rank);

              CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::sw()));
              CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_1_se_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_1_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_1_nw_patch->rank);
              CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::ne()));

              CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_1_nw_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_1_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_1_se_patch->rank);
              CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::nw()));
              CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::ne()));

              CHECK(domain_1_ne_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_1_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_1_sw_patch->rank);
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::se()));
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::nw()));
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::ne()));

              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::sw()));
              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::se()));
              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::nw()));
              CHECK_FALSE(domain_0_coarser_patch->hasNbr(Corner<2>::ne()));
            }
          }
        }
      }
    }
  }
}
TEST_CASE("P4estDomainGenerator 2x2 Refined SW", "[p4estDomGen]")
{
  for (int nx : { 5, 10 }) {
    for (int ny : { 5, 10 }) {
      for (double scale_x : { 0.5, 1.0 }) {
        for (double scale_y : { 0.5, 1.0 }) {
          for (int num_ghost_cells : { 0, 1, 2 }) {
            p4est_connectivity_t* conn = p4est_connectivity_new_unitsquare();

            p4est_t* p4est = p4est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);

            p4est_refine(
              p4est, false, [](p4est_t* p4est, p4est_topidx_t witch_tree, p4est_quadrant_t* quadrant) -> int { return 1; }, nullptr);
            p4est_refine(
              p4est, false, [](p4est_t* p4est, p4est_topidx_t witch_tree, p4est_quadrant_t* quadrant) -> int { return (quadrant->x == 0 && quadrant->y == 0); }, nullptr);

            INFO("nx: " << nx);
            INFO("ny: " << ny);
            INFO("scale_x: " << scale_x);
            INFO("scale_x: " << scale_y);
            INFO("num_ghost_cells: " << num_ghost_cells);

            P4estDomainGenerator::BlockMapFunc bmf = [&](int block_no, double unit_x, double unit_y, double& x, double& y) {
              x = scale_x * unit_x;
              y = scale_y * unit_y;
            };

            P4estDomainGenerator dg(p4est, { nx, ny }, num_ghost_cells, bmf);

            CHECK(dg.hasCoarserDomain());
            Domain<2> domain_2 = dg.getFinestDomain();
            CHECK(dg.hasCoarserDomain());
            Domain<2> domain_1 = dg.getCoarserDomain();
            CHECK(dg.hasCoarserDomain());
            Domain<2> domain_0 = dg.getCoarserDomain();
            CHECK_FALSE(dg.hasCoarserDomain());

            // SECTION("patches have correct ns")
            {
              for (auto patch : domain_2.getPatchInfoVector()) {
                CHECK(patch.ns[0] == nx);
                CHECK(patch.ns[1] == ny);
              }

              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK(patch.ns[0] == nx);
                CHECK(patch.ns[1] == ny);
              }

              auto patch = domain_0.getPatchInfoVector()[0];

              CHECK(patch.ns[0] == nx);
              CHECK(patch.ns[1] == ny);
            }
            // SECTION("patches have ranks set")
            {
              for (auto patch : domain_2.getPatchInfoVector()) {
                CHECK(patch.rank == 0);
              }

              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK(patch.rank == 0);
              }

              auto patch = domain_0.getPatchInfoVector()[0];

              CHECK(patch.rank == 0);
            }

            // SECTION("patches have num_ghost_cells set")
            {
              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK(patch.num_ghost_cells == num_ghost_cells);
              }

              auto patch = domain_0.getPatchInfoVector()[0];

              CHECK(patch.num_ghost_cells == num_ghost_cells);
            }

            const PatchInfo<2>* domain_2_sw_sw_patch = nullptr;
            const PatchInfo<2>* domain_2_sw_se_patch = nullptr;
            const PatchInfo<2>* domain_2_sw_nw_patch = nullptr;
            const PatchInfo<2>* domain_2_sw_ne_patch = nullptr;
            const PatchInfo<2>* domain_2_se_patch = nullptr;
            const PatchInfo<2>* domain_2_nw_patch = nullptr;
            const PatchInfo<2>* domain_2_ne_patch = nullptr;

            const PatchInfo<2>* domain_1_sw_patch = nullptr;
            const PatchInfo<2>* domain_1_se_patch = nullptr;
            const PatchInfo<2>* domain_1_nw_patch = nullptr;
            const PatchInfo<2>* domain_1_ne_patch = nullptr;

            const PatchInfo<2>* domain_0_patch = nullptr;

            for (const PatchInfo<2>& patch : domain_2.getPatchInfoVector()) {
              double x = patch.starts[0];
              double y = patch.starts[1];
              if (x == Catch::Approx(0) && y == Catch::Approx(0)) {
                domain_2_sw_sw_patch = &patch;
              }
              if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0)) {
                domain_2_sw_se_patch = &patch;
              }
              if (x == Catch::Approx(0) && y == Catch::Approx(0.25 * scale_y)) {
                domain_2_sw_nw_patch = &patch;
              }
              if (x == Catch::Approx(0.25 * scale_x) && y == Catch::Approx(0.25 * scale_y)) {
                domain_2_sw_ne_patch = &patch;
              }
              if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0)) {
                domain_2_se_patch = &patch;
              }
              if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y)) {
                domain_2_nw_patch = &patch;
              }
              if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
                domain_2_ne_patch = &patch;
              }
            }

            for (const PatchInfo<2>& patch : domain_1.getPatchInfoVector()) {
              double x = patch.starts[0];
              double y = patch.starts[1];
              if (x == Catch::Approx(0) && y == Catch::Approx(0)) {
                domain_1_sw_patch = &patch;
              }
              if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0)) {
                domain_1_se_patch = &patch;
              }
              if (x == Catch::Approx(0) && y == Catch::Approx(0.5 * scale_y)) {
                domain_1_nw_patch = &patch;
              }
              if (x == Catch::Approx(0.5 * scale_x) && y == Catch::Approx(0.5 * scale_y)) {
                domain_1_ne_patch = &patch;
              }
            }

            domain_0_patch = &domain_0.getPatchInfoVector()[0];

            // SECTION("patches have starts set")
            {
              REQUIRE(domain_2_sw_sw_patch != nullptr);
              REQUIRE(domain_2_sw_se_patch != nullptr);
              REQUIRE(domain_2_sw_nw_patch != nullptr);
              REQUIRE(domain_2_sw_ne_patch != nullptr);
              REQUIRE(domain_2_se_patch != nullptr);
              REQUIRE(domain_2_nw_patch != nullptr);
              REQUIRE(domain_2_ne_patch != nullptr);

              REQUIRE(domain_1_sw_patch != nullptr);
              REQUIRE(domain_1_se_patch != nullptr);
              REQUIRE(domain_1_nw_patch != nullptr);
              REQUIRE(domain_1_ne_patch != nullptr);

              CHECK(domain_0_patch->starts[0] == Catch::Approx(0.0));
              CHECK(domain_0_patch->starts[1] == Catch::Approx(0.0));
            }

            // SECTION("patches have correct spacings")
            {
              CHECK(domain_2_sw_sw_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
              CHECK(domain_2_sw_sw_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
              CHECK(domain_2_sw_se_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
              CHECK(domain_2_sw_se_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
              CHECK(domain_2_sw_nw_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
              CHECK(domain_2_sw_nw_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));
              CHECK(domain_2_sw_ne_patch->spacings[0] == Catch::Approx(scale_x * 0.25 / nx));
              CHECK(domain_2_sw_ne_patch->spacings[1] == Catch::Approx(scale_y * 0.25 / ny));

              CHECK(domain_2_se_patch->spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
              CHECK(domain_2_se_patch->spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
              CHECK(domain_2_nw_patch->spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
              CHECK(domain_2_nw_patch->spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
              CHECK(domain_2_ne_patch->spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
              CHECK(domain_2_ne_patch->spacings[1] == Catch::Approx(scale_y * 0.5 / ny));

              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK(patch.spacings[0] == Catch::Approx(scale_x * 0.5 / nx));
                CHECK(patch.spacings[1] == Catch::Approx(scale_y * 0.5 / ny));
              }

              auto patch = domain_0.getPatchInfoVector()[0];

              CHECK(patch.spacings[0] == Catch::Approx(scale_x * 1.0 / nx));
              CHECK(patch.spacings[1] == Catch::Approx(scale_y * 1.0 / ny));
            }

            // SECTION("patches have refine_level set")
            {
              CHECK(domain_2_sw_sw_patch->refine_level == 2);
              CHECK(domain_2_sw_se_patch->refine_level == 2);
              CHECK(domain_2_sw_nw_patch->refine_level == 2);
              CHECK(domain_2_sw_ne_patch->refine_level == 2);

              CHECK(domain_2_se_patch->refine_level == 1);
              CHECK(domain_2_nw_patch->refine_level == 1);
              CHECK(domain_2_ne_patch->refine_level == 1);

              for (auto patch : domain_1.getPatchInfoVector()) {
                CHECK(patch.refine_level == 1);
              }

              auto patch = domain_0.getPatchInfoVector()[0];

              CHECK(patch.refine_level == 0);
            }

            // SECTION("parent ids are set correctly")
            {
              CHECK(domain_2_sw_sw_patch->parent_id == domain_1_sw_patch->id);
              CHECK(domain_2_sw_se_patch->parent_id == domain_1_sw_patch->id);
              CHECK(domain_2_sw_nw_patch->parent_id == domain_1_sw_patch->id);
              CHECK(domain_2_sw_ne_patch->parent_id == domain_1_sw_patch->id);

              CHECK(domain_2_se_patch->parent_id == domain_1_se_patch->id);
              CHECK(domain_2_nw_patch->parent_id == domain_1_nw_patch->id);
              CHECK(domain_2_ne_patch->parent_id == domain_1_ne_patch->id);

              CHECK(domain_1_sw_patch->parent_id == domain_0_patch->id);
              CHECK(domain_1_se_patch->parent_id == domain_0_patch->id);
              CHECK(domain_1_nw_patch->parent_id == domain_0_patch->id);
              CHECK(domain_1_ne_patch->parent_id == domain_0_patch->id);

              CHECK(domain_0_patch->parent_id == -1);
            }
            // SECTION("child ids are set correctly")
            {
              for (int i = 0; i < 4; i++) {
                CHECK(domain_2_sw_sw_patch->child_ids[i] == -1);
                CHECK(domain_2_sw_se_patch->child_ids[i] == -1);
                CHECK(domain_2_sw_nw_patch->child_ids[i] == -1);
                CHECK(domain_2_sw_ne_patch->child_ids[i] == -1);
                CHECK(domain_2_se_patch->child_ids[i] == -1);
                CHECK(domain_2_nw_patch->child_ids[i] == -1);
                CHECK(domain_2_ne_patch->child_ids[i] == -1);
              }

              CHECK(domain_1_sw_patch->child_ids[0] == domain_2_sw_sw_patch->id);
              CHECK(domain_1_sw_patch->child_ids[1] == domain_2_sw_se_patch->id);
              CHECK(domain_1_sw_patch->child_ids[2] == domain_2_sw_nw_patch->id);
              CHECK(domain_1_sw_patch->child_ids[3] == domain_2_sw_ne_patch->id);

              CHECK(domain_1_se_patch->child_ids[0] == domain_2_se_patch->id);
              CHECK(domain_1_nw_patch->child_ids[0] == domain_2_nw_patch->id);
              CHECK(domain_1_ne_patch->child_ids[0] == domain_2_ne_patch->id);

              for (int i = 1; i < 4; i++) {
                CHECK(domain_1_se_patch->child_ids[i] == -1);
                CHECK(domain_1_nw_patch->child_ids[i] == -1);
                CHECK(domain_1_ne_patch->child_ids[i] == -1);
              }

              CHECK(domain_0_patch->child_ids[0] == domain_1_sw_patch->id);
              CHECK(domain_0_patch->child_ids[1] == domain_1_se_patch->id);
              CHECK(domain_0_patch->child_ids[2] == domain_1_nw_patch->id);
              CHECK(domain_0_patch->child_ids[3] == domain_1_ne_patch->id);
            }
            // SECTION("orth on parent is set correctly")
            {
              CHECK(domain_2_sw_sw_patch->orth_on_parent == Orthant<2>::sw());
              CHECK(domain_2_sw_se_patch->orth_on_parent == Orthant<2>::se());
              CHECK(domain_2_sw_nw_patch->orth_on_parent == Orthant<2>::nw());
              CHECK(domain_2_sw_ne_patch->orth_on_parent == Orthant<2>::ne());
              CHECK(domain_2_se_patch->orth_on_parent == Orthant<2>::null());
              CHECK(domain_2_nw_patch->orth_on_parent == Orthant<2>::null());
              CHECK(domain_2_ne_patch->orth_on_parent == Orthant<2>::null());

              CHECK(domain_1_sw_patch->orth_on_parent == Orthant<2>::sw());
              CHECK(domain_1_se_patch->orth_on_parent == Orthant<2>::se());
              CHECK(domain_1_nw_patch->orth_on_parent == Orthant<2>::nw());
              CHECK(domain_1_ne_patch->orth_on_parent == Orthant<2>::ne());

              CHECK(domain_0_patch->orth_on_parent == Orthant<2>::null());
            }
            // SECTION("parent ranks are set correctly")
            {
              CHECK(domain_2_sw_sw_patch->parent_rank == domain_1_sw_patch->rank);
              CHECK(domain_2_sw_se_patch->parent_rank == domain_1_sw_patch->rank);
              CHECK(domain_2_sw_nw_patch->parent_rank == domain_1_sw_patch->rank);
              CHECK(domain_2_sw_ne_patch->parent_rank == domain_1_sw_patch->rank);

              CHECK(domain_2_se_patch->parent_rank == domain_1_se_patch->rank);
              CHECK(domain_2_nw_patch->parent_rank == domain_1_nw_patch->rank);
              CHECK(domain_2_ne_patch->parent_rank == domain_1_ne_patch->rank);

              CHECK(domain_1_sw_patch->parent_rank == domain_0_patch->rank);
              CHECK(domain_1_se_patch->parent_rank == domain_0_patch->rank);
              CHECK(domain_1_nw_patch->parent_rank == domain_0_patch->rank);
              CHECK(domain_1_ne_patch->parent_rank == domain_0_patch->rank);

              CHECK(domain_0_patch->parent_rank == -1);
            }
            // SECTION("child ranks are set correctly")
            {
              for (int i = 0; i < 4; i++) {
                CHECK(domain_2_sw_sw_patch->child_ranks[i] == -1);
                CHECK(domain_2_sw_se_patch->child_ranks[i] == -1);
                CHECK(domain_2_sw_nw_patch->child_ranks[i] == -1);
                CHECK(domain_2_sw_ne_patch->child_ranks[i] == -1);
                CHECK(domain_2_se_patch->child_ranks[i] == -1);
                CHECK(domain_2_nw_patch->child_ranks[i] == -1);
                CHECK(domain_2_ne_patch->child_ranks[i] == -1);
              }

              CHECK(domain_1_sw_patch->child_ranks[0] == domain_2_sw_sw_patch->rank);
              CHECK(domain_1_sw_patch->child_ranks[1] == domain_2_sw_se_patch->rank);
              CHECK(domain_1_sw_patch->child_ranks[2] == domain_2_sw_nw_patch->rank);
              CHECK(domain_1_sw_patch->child_ranks[3] == domain_2_sw_ne_patch->rank);

              CHECK(domain_1_se_patch->child_ranks[0] == domain_2_se_patch->rank);
              CHECK(domain_1_nw_patch->child_ranks[0] == domain_2_nw_patch->rank);
              CHECK(domain_1_ne_patch->child_ranks[0] == domain_2_ne_patch->rank);

              for (int i = 1; i < 4; i++) {
                CHECK(domain_1_se_patch->child_ranks[i] == -1);
                CHECK(domain_1_nw_patch->child_ranks[i] == -1);
                CHECK(domain_1_ne_patch->child_ranks[i] == -1);
              }

              CHECK(domain_0_patch->child_ranks[0] == domain_1_sw_patch->rank);
              CHECK(domain_0_patch->child_ranks[1] == domain_1_se_patch->rank);
              CHECK(domain_0_patch->child_ranks[2] == domain_1_nw_patch->rank);
              CHECK(domain_0_patch->child_ranks[3] == domain_1_ne_patch->rank);
            }
            // SECTION("nbr_info ids are correct")
            {
              CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Side<2>::west()));
              CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_sw_se_patch->id);
              CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Side<2>::south()));
              CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_sw_nw_patch->id);

              CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_sw_sw_patch->id);
              CHECK(domain_2_sw_se_patch->getCoarseNbrInfo(Side<2>::east()).id == domain_2_se_patch->id);
              CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Side<2>::south()));
              CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_sw_ne_patch->id);

              CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Side<2>::west()));
              CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_sw_ne_patch->id);
              CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_sw_sw_patch->id);
              CHECK(domain_2_sw_nw_patch->getCoarseNbrInfo(Side<2>::north()).id == domain_2_nw_patch->id);

              CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_sw_nw_patch->id);
              CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Side<2>::east()).id == domain_2_se_patch->id);
              CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_sw_se_patch->id);
              CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Side<2>::north()).id == domain_2_nw_patch->id);

              CHECK(domain_2_se_patch->getFineNbrInfo(Side<2>::west()).ids[0] == domain_2_sw_se_patch->id);
              CHECK(domain_2_se_patch->getFineNbrInfo(Side<2>::west()).ids[1] == domain_2_sw_ne_patch->id);
              CHECK_FALSE(domain_2_se_patch->hasNbr(Side<2>::east()));
              CHECK_FALSE(domain_2_se_patch->hasNbr(Side<2>::south()));
              CHECK(domain_2_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_2_ne_patch->id);

              CHECK_FALSE(domain_2_nw_patch->hasNbr(Side<2>::west()));
              CHECK(domain_2_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_2_ne_patch->id);
              CHECK(domain_2_nw_patch->getFineNbrInfo(Side<2>::south()).ids[0] == domain_2_sw_nw_patch->id);
              CHECK(domain_2_nw_patch->getFineNbrInfo(Side<2>::south()).ids[1] == domain_2_sw_ne_patch->id);
              CHECK_FALSE(domain_2_nw_patch->hasNbr(Side<2>::north()));

              CHECK(domain_2_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_2_nw_patch->id);
              CHECK_FALSE(domain_2_ne_patch->hasNbr(Side<2>::east()));
              CHECK(domain_2_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_2_se_patch->id);
              CHECK_FALSE(domain_2_ne_patch->hasNbr(Side<2>::north()));

              CHECK_FALSE(domain_1_sw_patch->hasNbr(Side<2>::west()));
              CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_1_se_patch->id);
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Side<2>::south()));
              CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::north()).id == domain_1_nw_patch->id);

              CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::west()).id == domain_1_sw_patch->id);
              CHECK_FALSE(domain_1_se_patch->hasNbr(Side<2>::east()));
              CHECK_FALSE(domain_1_se_patch->hasNbr(Side<2>::south()));
              CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::north()).id == domain_1_ne_patch->id);

              CHECK_FALSE(domain_1_nw_patch->hasNbr(Side<2>::west()));
              CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::east()).id == domain_1_ne_patch->id);
              CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::south()).id == domain_1_sw_patch->id);
              CHECK_FALSE(domain_1_nw_patch->hasNbr(Side<2>::north()));

              CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::west()).id == domain_1_nw_patch->id);
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Side<2>::east()));
              CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::south()).id == domain_1_se_patch->id);
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Side<2>::north()));

              CHECK_FALSE(domain_0_patch->hasNbr(Side<2>::west()));
              CHECK_FALSE(domain_0_patch->hasNbr(Side<2>::east()));
              CHECK_FALSE(domain_0_patch->hasNbr(Side<2>::south()));
              CHECK_FALSE(domain_0_patch->hasNbr(Side<2>::north()));
            }
            // SECTION("nbr_info ranks are correct")
            {
              CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_sw_se_patch->rank);
              CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_sw_nw_patch->rank);

              CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_sw_sw_patch->rank);
              CHECK(domain_2_sw_se_patch->getCoarseNbrInfo(Side<2>::east()).rank == domain_2_se_patch->rank);
              CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_sw_ne_patch->rank);

              CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_sw_ne_patch->rank);
              CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_sw_sw_patch->rank);
              CHECK(domain_2_sw_nw_patch->getCoarseNbrInfo(Side<2>::north()).rank == domain_2_nw_patch->rank);

              CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_sw_nw_patch->rank);
              CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Side<2>::east()).rank == domain_2_se_patch->rank);
              CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_sw_se_patch->rank);
              CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Side<2>::north()).rank == domain_2_nw_patch->rank);

              CHECK(domain_2_se_patch->getFineNbrInfo(Side<2>::west()).ranks[0] == domain_2_sw_se_patch->rank);
              CHECK(domain_2_se_patch->getFineNbrInfo(Side<2>::west()).ranks[1] == domain_2_sw_ne_patch->rank);
              CHECK(domain_2_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_2_ne_patch->rank);

              CHECK(domain_2_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_2_ne_patch->rank);
              CHECK(domain_2_nw_patch->getFineNbrInfo(Side<2>::south()).ranks[0] == domain_2_sw_nw_patch->rank);
              CHECK(domain_2_nw_patch->getFineNbrInfo(Side<2>::south()).ranks[1] == domain_2_sw_ne_patch->rank);

              CHECK(domain_2_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_2_nw_patch->rank);
              CHECK(domain_2_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_2_se_patch->rank);

              CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_se_patch->rank);
              CHECK(domain_1_sw_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_nw_patch->rank);

              CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_sw_patch->rank);
              CHECK(domain_1_se_patch->getNormalNbrInfo(Side<2>::north()).rank == domain_1_ne_patch->rank);

              CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::east()).rank == domain_1_ne_patch->rank);
              CHECK(domain_1_nw_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_sw_patch->rank);

              CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::west()).rank == domain_1_nw_patch->rank);
              CHECK(domain_1_ne_patch->getNormalNbrInfo(Side<2>::south()).rank == domain_1_se_patch->rank);
            }
            // SECTION("nbr_info orth_on_coarse are correct")
            {
              CHECK(domain_2_sw_se_patch->getCoarseNbrInfo(Side<2>::east()).orth_on_coarse == Orthant<1>::lower());

              CHECK(domain_2_sw_nw_patch->getCoarseNbrInfo(Side<2>::north()).orth_on_coarse == Orthant<1>::lower());

              CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Side<2>::east()).orth_on_coarse == Orthant<1>::upper());
              CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Side<2>::north()).orth_on_coarse == Orthant<1>::upper());
            }
            // SECTION("corner_nbr_info ids are correct")
            {
              CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Corner<2>::sw()));
              CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Corner<2>::se()));
              CHECK_FALSE(domain_2_sw_sw_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_2_sw_sw_patch->hasNbr(Corner<2>::ne()));
              CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_2_sw_ne_patch->id);

              CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Corner<2>::sw()));
              CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_2_sw_se_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_sw_nw_patch->id);
              CHECK_FALSE(domain_2_sw_se_patch->hasNbr(Corner<2>::ne()));

              CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_2_sw_nw_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_sw_se_patch->id);
              CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Corner<2>::nw()));
              CHECK_FALSE(domain_2_sw_nw_patch->hasNbr(Corner<2>::ne()));

              CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_2_sw_sw_patch->id);
              CHECK_FALSE(domain_2_sw_ne_patch->hasNbr(Corner<2>::se()));
              CHECK_FALSE(domain_2_sw_ne_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::ne()));
              CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Corner<2>::ne()).id == domain_2_ne_patch->id);

              CHECK_FALSE(domain_2_se_patch->hasNbr(Corner<2>::sw()));
              CHECK_FALSE(domain_2_se_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_2_se_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_2_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_2_nw_patch->id);
              CHECK_FALSE(domain_2_se_patch->hasNbr(Corner<2>::ne()));

              CHECK_FALSE(domain_2_nw_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_2_nw_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_2_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_2_se_patch->id);
              CHECK_FALSE(domain_2_nw_patch->hasNbr(Corner<2>::nw()));
              CHECK_FALSE(domain_2_nw_patch->hasNbr(Corner<2>::ne()));

              CHECK(domain_2_ne_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_2_ne_patch->getFineNbrInfo(Corner<2>::sw()).ids[0] == domain_2_sw_ne_patch->id);
              CHECK_FALSE(domain_2_ne_patch->hasNbr(Corner<2>::se()));
              CHECK_FALSE(domain_2_ne_patch->hasNbr(Corner<2>::nw()));
              CHECK_FALSE(domain_2_ne_patch->hasNbr(Corner<2>::ne()));

              CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::sw()));
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::se()));
              CHECK_FALSE(domain_1_sw_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_1_sw_patch->hasNbr(Corner<2>::ne()));
              CHECK(domain_1_sw_patch->getNormalNbrInfo(Corner<2>::ne()).id == domain_1_ne_patch->id);

              CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::sw()));
              CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_1_se_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_1_se_patch->getNormalNbrInfo(Corner<2>::nw()).id == domain_1_nw_patch->id);
              CHECK_FALSE(domain_1_se_patch->hasNbr(Corner<2>::ne()));

              CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_1_nw_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_1_nw_patch->getNormalNbrInfo(Corner<2>::se()).id == domain_1_se_patch->id);
              CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::nw()));
              CHECK_FALSE(domain_1_nw_patch->hasNbr(Corner<2>::ne()));

              CHECK(domain_1_ne_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_1_ne_patch->getNormalNbrInfo(Corner<2>::sw()).id == domain_1_sw_patch->id);
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::se()));
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::nw()));
              CHECK_FALSE(domain_1_ne_patch->hasNbr(Corner<2>::ne()));

              CHECK_FALSE(domain_0_patch->hasNbr(Corner<2>::sw()));
              CHECK_FALSE(domain_0_patch->hasNbr(Corner<2>::se()));
              CHECK_FALSE(domain_0_patch->hasNbr(Corner<2>::nw()));
              CHECK_FALSE(domain_0_patch->hasNbr(Corner<2>::ne()));
            }
            // SECTION("corner_nbr_info ranks are correct")
            {
              CHECK(domain_2_sw_sw_patch->hasNbr(Corner<2>::ne()));
              CHECK(domain_2_sw_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_2_sw_ne_patch->rank);

              CHECK(domain_2_sw_se_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_2_sw_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_sw_nw_patch->rank);

              CHECK(domain_2_sw_nw_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_2_sw_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_sw_se_patch->rank);

              CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_2_sw_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_2_sw_sw_patch->rank);
              CHECK(domain_2_sw_ne_patch->hasNbr(Corner<2>::ne()));
              CHECK(domain_2_sw_ne_patch->getCoarseNbrInfo(Corner<2>::ne()).rank == domain_2_ne_patch->rank);

              CHECK(domain_2_se_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_2_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_2_nw_patch->rank);

              CHECK(domain_2_nw_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_2_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_2_se_patch->rank);

              CHECK(domain_2_ne_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_2_ne_patch->getFineNbrInfo(Corner<2>::sw()).ranks[0] == domain_2_sw_ne_patch->rank);

              CHECK(domain_1_sw_patch->hasNbr(Corner<2>::ne()));
              CHECK(domain_1_sw_patch->getNormalNbrInfo(Corner<2>::ne()).rank == domain_1_ne_patch->rank);

              CHECK(domain_1_se_patch->hasNbr(Corner<2>::nw()));
              CHECK(domain_1_se_patch->getNormalNbrInfo(Corner<2>::nw()).rank == domain_1_nw_patch->rank);

              CHECK(domain_1_nw_patch->hasNbr(Corner<2>::se()));
              CHECK(domain_1_nw_patch->getNormalNbrInfo(Corner<2>::se()).rank == domain_1_se_patch->rank);

              CHECK(domain_1_ne_patch->hasNbr(Corner<2>::sw()));
              CHECK(domain_1_ne_patch->getNormalNbrInfo(Corner<2>::sw()).rank == domain_1_sw_patch->rank);
            }
          }
        }
      }
    }
  }
}
TEST_CASE("2x1 brick", "[p4estDomGen]")
{
  p4est_connectivity_t* conn = p4est_connectivity_new_brick(2, 1, false, false);
  p4est_t* p4est = p4est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, nullptr);
  P4estDomainGenerator::BlockMapFunc bmf = [](int block_no, double unit_x, double unit_y, double& x, double& y) {
    x = unit_x;
    y = unit_y;
  };
  P4estDomainGenerator dg(p4est, { 10, 10 }, 1, bmf);
  Domain<2> domain = dg.getFinestDomain();
  CHECK(domain.getNumGlobalPatches() == 2);
  auto patch1 = domain.getPatchInfoVector()[0];
  CHECK_FALSE(patch1.hasNbr(Side<2>::west()));
  CHECK(patch1.hasNbr(Side<2>::east()));
  CHECK_FALSE(patch1.hasNbr(Side<2>::south()));
  CHECK_FALSE(patch1.hasNbr(Side<2>::north()));
  CHECK(patch1.spacings[0] == Catch::Approx(1.0 / 10));
  CHECK(patch1.spacings[1] == Catch::Approx(1.0 / 10));
  CHECK(patch1.starts[0] == Catch::Approx(0));
  CHECK(patch1.starts[1] == Catch::Approx(0));
  CHECK(patch1.ns[0] == 10);
  CHECK(patch1.ns[1] == 10);
  auto patch2 = domain.getPatchInfoVector()[1];
  CHECK(patch2.hasNbr(Side<2>::west()));
  CHECK_FALSE(patch2.hasNbr(Side<2>::east()));
  CHECK_FALSE(patch2.hasNbr(Side<2>::south()));
  CHECK_FALSE(patch2.hasNbr(Side<2>::north()));
  CHECK(patch2.spacings[0] == Catch::Approx(1.0 / 10));
  CHECK(patch2.spacings[1] == Catch::Approx(1.0 / 10));
  CHECK(patch2.starts[0] == Catch::Approx(0));
  CHECK(patch2.starts[1] == Catch::Approx(0));
  CHECK(patch2.ns[0] == 10);
  CHECK(patch2.ns[1] == 10);

  CHECK_FALSE(dg.hasCoarserDomain());
}
