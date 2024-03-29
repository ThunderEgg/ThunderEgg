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
#include "utils/DomainReader.h"
#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/DomainTools.h>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;
using namespace doctest;

constexpr auto single_mesh_file = "mesh_inputs/2d_uniform_2x2_mpi1.json";
constexpr auto refined_mesh_file = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json";
constexpr auto cross_mesh_file = "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json";
TEST_CASE("exchange uniform 2D quad BiLinearGhostFiller")
{
  for (auto nx : { 1, 2, 10, 13 }) {
    for (auto ny : { 1, 2, 10, 13 }) {
      Communicator comm(MPI_COMM_WORLD);

      double startx = 0;
      double starty = 0;
      double lengthx = 1;
      double lengthy = 1;
      int num_ghost = 1;

      auto f = [&](const std::array<double, 2> coord) -> double {
        double x = coord[0];
        double y = coord[1];
        return 1 + ((x + 0.3) * y);
      };

      /*
       *  2 | 3
       *  --+---
       *  0 | 1
       */
      vector<PatchInfo<2>> pinfos(4);

      pinfos[0].id = 0;
      pinfos[0].ns = { nx, ny };
      pinfos[0].spacings = { lengthx / nx, lengthy / ny };
      pinfos[0].starts = { startx, starty };
      pinfos[0].num_ghost_cells = num_ghost;
      pinfos[0].setNbrInfo(Side<2>::east(), new NormalNbrInfo<1>(1));
      pinfos[0].setNbrInfo(Side<2>::north(), new NormalNbrInfo<1>(2));

      pinfos[1].id = 1;
      pinfos[1].ns = { nx, ny };
      pinfos[1].spacings = { lengthx / nx, lengthy / ny };
      pinfos[1].starts = { startx + lengthx, starty };
      pinfos[1].num_ghost_cells = num_ghost;
      pinfos[1].setNbrInfo(Side<2>::west(), new NormalNbrInfo<1>(0));
      pinfos[1].setNbrInfo(Side<2>::north(), new NormalNbrInfo<1>(3));

      pinfos[2].id = 2;
      pinfos[2].ns = { nx, ny };
      pinfos[2].spacings = { lengthx / nx, lengthy / ny };
      pinfos[2].starts = { startx, starty + lengthy };
      pinfos[2].num_ghost_cells = num_ghost;
      pinfos[2].setNbrInfo(Side<2>::east(), new NormalNbrInfo<1>(3));
      pinfos[2].setNbrInfo(Side<2>::south(), new NormalNbrInfo<1>(0));

      pinfos[3].id = 3;
      pinfos[3].ns = { nx, ny };
      pinfos[3].spacings = { lengthx / nx, lengthy / ny };
      pinfos[3].starts = { startx + lengthx, starty + lengthy };
      pinfos[3].num_ghost_cells = num_ghost;
      pinfos[3].setNbrInfo(Side<2>::west(), new NormalNbrInfo<1>(2));
      pinfos[3].setNbrInfo(Side<2>::south(), new NormalNbrInfo<1>(1));

      Domain<2> d(comm, 0, { nx, ny }, num_ghost, pinfos.begin(), pinfos.end());

      Vector<2> vec(d, 1);

      DomainTools::SetValues<2>(d, vec, f);

      BiLinearGhostFiller blgf(d, GhostFillingType::Faces);
      blgf.fillGhost(vec);

      // patch 1
      {
        // check that center values weren't modified
        auto patch_1 = vec.getComponentView(0, d.getPatchInfoVector()[0].local_index);
        Loop::Nested<2>(patch_1.getStart(), patch_1.getEnd(), [&](const std::array<int, 2> coord) {
          std::array<double, 2> real_coord;
          DomainTools::GetRealCoord<2>(pinfos[0], coord, real_coord);
          CHECK_EQ(patch_1[coord], f(real_coord));
        });
        // check that west values are not modified
        {
          auto west_ghost = patch_1.getSliceOn(Side<2>::west(), { -1 });
          Loop::Nested<1>(west_ghost.getStart(), west_ghost.getEnd(), [&](const std::array<int, 1> coord) { CHECK_EQ(west_ghost[coord], 0); });
        }
        // check that east values are correct
        {
          auto east_ghost = patch_1.getSliceOn(Side<2>::east(), { -1 });
          Loop::Nested<1>(east_ghost.getStart(), east_ghost.getEnd(), [&](const std::array<int, 1> coord) {
            std::array<double, 2> real_coord;
            DomainTools::GetRealCoordGhost<2>(pinfos[0], { nx, coord[0] }, real_coord);
            CHECK_EQ(east_ghost[coord], Approx(f(real_coord)));
          });
        }
        // check that south values are not modified
        {
          auto south_ghost = patch_1.getSliceOn(Side<2>::south(), { -1 });
          Loop::Nested<1>(south_ghost.getStart(), south_ghost.getEnd(), [&](const std::array<int, 1> coord) { CHECK_EQ(south_ghost[coord], 0); });
        }
        // check that north values are correct
        {
          auto north_ghost = patch_1.getSliceOn(Side<2>::north(), { -1 });
          Loop::Nested<1>(north_ghost.getStart(), north_ghost.getEnd(), [&](const std::array<int, 1> coord) {
            std::array<double, 2> real_coord;
            DomainTools::GetRealCoordGhost<2>(pinfos[0], { coord[0], ny }, real_coord);
            CHECK_EQ(north_ghost[coord], Approx(f(real_coord)));
          });
        }
      }
      // patch 2
      {
        // check that center values weren't modified
        auto patch_2 = vec.getComponentView(0, d.getPatchInfoVector()[1].local_index);
        Loop::Nested<2>(patch_2.getStart(), patch_2.getEnd(), [&](const std::array<int, 2> coord) {
          std::array<double, 2> real_coord;
          DomainTools::GetRealCoord<2>(pinfos[1], coord, real_coord);
          CHECK_EQ(patch_2[coord], f(real_coord));
        });
        // check that west values correct
        {
          auto west_ghost = patch_2.getSliceOn(Side<2>::west(), { -1 });
          Loop::Nested<1>(west_ghost.getStart(), west_ghost.getEnd(), [&](const std::array<int, 1> coord) {
            std::array<double, 2> real_coord;
            DomainTools::GetRealCoordGhost<2>(pinfos[1], { -1, coord[0] }, real_coord);
            CHECK_EQ(west_ghost[coord], Approx(f(real_coord)));
          });
        }
        // check that east values are not modified
        {
          auto east_ghost = patch_2.getSliceOn(Side<2>::east(), { -1 });
          Loop::Nested<1>(east_ghost.getStart(), east_ghost.getEnd(), [&](const std::array<int, 1> coord) { CHECK_EQ(east_ghost[coord], 0); });
        }
        // check that south values are not modified
        {
          auto south_ghost = patch_2.getSliceOn(Side<2>::south(), { -1 });
          Loop::Nested<1>(south_ghost.getStart(), south_ghost.getEnd(), [&](const std::array<int, 1> coord) { CHECK_EQ(south_ghost[coord], 0); });
        }
        // check that north values are correct
        {
          auto north_ghost = patch_2.getSliceOn(Side<2>::north(), { -1 });
          Loop::Nested<1>(north_ghost.getStart(), north_ghost.getEnd(), [&](const std::array<int, 1> coord) {
            std::array<double, 2> real_coord;
            DomainTools::GetRealCoordGhost<2>(pinfos[1], { coord[0], ny }, real_coord);
            CHECK_EQ(north_ghost[coord], Approx(f(real_coord)));
          });
        }
      }
      // patch 3
      {
        // check that center values weren't modified
        auto patch_3 = vec.getComponentView(0, d.getPatchInfoVector()[2].local_index);
        Loop::Nested<2>(patch_3.getStart(), patch_3.getEnd(), [&](const std::array<int, 2> coord) {
          std::array<double, 2> real_coord;
          DomainTools::GetRealCoord<2>(pinfos[2], coord, real_coord);
          CHECK_EQ(patch_3[coord], f(real_coord));
        });
        // check that west values are not modified
        {
          auto west_ghost = patch_3.getSliceOn(Side<2>::west(), { -1 });
          Loop::Nested<1>(west_ghost.getStart(), west_ghost.getEnd(), [&](const std::array<int, 1> coord) { CHECK_EQ(west_ghost[coord], 0); });
        }
        // check that east values are correct
        {
          auto east_ghost = patch_3.getSliceOn(Side<2>::east(), { -1 });
          Loop::Nested<1>(east_ghost.getStart(), east_ghost.getEnd(), [&](const std::array<int, 1> coord) {
            std::array<double, 2> real_coord;
            DomainTools::GetRealCoordGhost<2>(pinfos[2], { nx, coord[0] }, real_coord);
            CHECK_EQ(east_ghost[coord], Approx(f(real_coord)));
          });
        }
        // check that south values are not modified
        {
          auto south_ghost = patch_3.getSliceOn(Side<2>::south(), { -1 });
          Loop::Nested<1>(south_ghost.getStart(), south_ghost.getEnd(), [&](const std::array<int, 1> coord) {
            std::array<double, 2> real_coord;
            DomainTools::GetRealCoordGhost<2>(pinfos[2], { coord[0], -1 }, real_coord);
            CHECK_EQ(south_ghost[coord], Approx(f(real_coord)));
          });
        }
        // check that north values are correct
        {
          auto north_ghost = patch_3.getSliceOn(Side<2>::north(), { -1 });
          Loop::Nested<1>(north_ghost.getStart(), north_ghost.getEnd(), [&](const std::array<int, 1> coord) { CHECK_EQ(north_ghost[coord], 0); });
        }
      }
      // patch 4
      {
        // check that center values weren't modified
        auto patch_4 = vec.getComponentView(0, d.getPatchInfoVector()[3].local_index);
        Loop::Nested<2>(patch_4.getStart(), patch_4.getEnd(), [&](const std::array<int, 2> coord) {
          std::array<double, 2> real_coord;
          DomainTools::GetRealCoord<2>(pinfos[3], coord, real_coord);
          CHECK_EQ(patch_4[coord], f(real_coord));
        });
        // check that west values correct
        {
          auto west_ghost = patch_4.getSliceOn(Side<2>::west(), { -1 });
          Loop::Nested<1>(west_ghost.getStart(), west_ghost.getEnd(), [&](const std::array<int, 1> coord) {
            std::array<double, 2> real_coord;
            DomainTools::GetRealCoordGhost<2>(pinfos[3], { -1, coord[0] }, real_coord);
            CHECK_EQ(west_ghost[coord], Approx(f(real_coord)));
          });
        }
        // check that east values are not modified
        {
          auto east_ghost = patch_4.getSliceOn(Side<2>::east(), { -1 });
          Loop::Nested<1>(east_ghost.getStart(), east_ghost.getEnd(), [&](const std::array<int, 1> coord) { CHECK_EQ(east_ghost[coord], 0); });
        }
        // check that south values are correct
        {
          auto south_ghost = patch_4.getSliceOn(Side<2>::south(), { -1 });
          Loop::Nested<1>(south_ghost.getStart(), south_ghost.getEnd(), [&](const std::array<int, 1> coord) {
            std::array<double, 2> real_coord;
            DomainTools::GetRealCoordGhost<2>(pinfos[3], { coord[0], -1 }, real_coord);
            CHECK_EQ(south_ghost[coord], Approx(f(real_coord)));
          });
        }
        // check that north values are not modified
        {
          auto north_ghost = patch_4.getSliceOn(Side<2>::north(), { -1 });
          Loop::Nested<1>(north_ghost.getStart(), north_ghost.getEnd(), [&](const std::array<int, 1> coord) { CHECK_EQ(north_ghost[coord], 0); });
        }
      }
    }
  }
}
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;

        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d = domain_reader.getFinerDomain();

        Vector<2> vec(d, 1);
        Vector<2> expected(d, 1);

        auto f = [&](const std::array<double, 2> coord) -> double {
          double x = coord[0];
          double y = coord[1];
          return 1 + ((x * 0.3) + y);
        };

        DomainTools::SetValues<2>(d, vec, f);
        DomainTools::SetValuesWithGhost<2>(d, expected, f);

        BiLinearGhostFiller blgf(d, GhostFillingType::Faces);
        blgf.fillGhost(vec);

        for (auto pinfo : d.getPatchInfoVector()) {
          ComponentView<double, 2> vec_ld = vec.getComponentView(0, pinfo.local_index);
          ComponentView<double, 2> expected_ld = expected.getComponentView(0, pinfo.local_index);
          Loop::Nested<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2>& coord) {
            ///
            REQUIRE_EQ(vec_ld[coord], Approx(expected_ld[coord]));
          });
          for (Side<2> s : Side<2>::getValues()) {
            View<double, 1> vec_ghost = vec_ld.getSliceOn(s, { -1 });
            View<double, 1> expected_ghost = expected_ld.getSliceOn(s, { -1 });
            if (pinfo.hasNbr(s)) {
              Loop::Nested<1>(vec_ghost.getStart(), vec_ghost.getEnd(), [&](const array<int, 1>& coord) {
                ///
                CHECK_EQ(vec_ghost[coord], Approx(expected_ghost[coord]));
              });
            }
          }
        }
      }
    }
  }
}
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller ghost already set")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;

        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d = domain_reader.getFinerDomain();

        Vector<2> vec(d, 1);
        Vector<2> expected(d, 1);

        auto f = [&](const std::array<double, 2> coord) -> double {
          double x = coord[0];
          double y = coord[0];
          return 1 + ((x * 0.3) + y);
        };

        DomainTools::SetValuesWithGhost<2>(d, vec, f);
        DomainTools::SetValuesWithGhost<2>(d, expected, f);

        BiLinearGhostFiller blgf(d, GhostFillingType::Faces);
        blgf.fillGhost(vec);

        for (auto pinfo : d.getPatchInfoVector()) {
          ComponentView<double, 2> vec_ld = vec.getComponentView(0, pinfo.local_index);
          ComponentView<double, 2> expected_ld = expected.getComponentView(0, pinfo.local_index);
          Loop::Nested<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2>& coord) {
            ///
            REQUIRE_EQ(vec_ld[coord], Approx(expected_ld[coord]));
          });
          for (Side<2> s : Side<2>::getValues()) {
            View<double, 1> vec_ghost = vec_ld.getSliceOn(s, { -1 });
            View<double, 1> expected_ghost = expected_ld.getSliceOn(s, { -1 });
            if (pinfo.hasNbr(s)) {
              Loop::Nested<1>(vec_ghost.getStart(), vec_ghost.getEnd(), [&](const array<int, 1>& coord) {
                ///
                CHECK_EQ(vec_ghost[coord], Approx(expected_ghost[coord]));
              });
            }
          }
        }
      }
    }
  }
}
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller two components")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;

        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d = domain_reader.getFinerDomain();

        Vector<2> vec(d, 2);
        Vector<2> expected(d, 2);

        auto f = [&](const std::array<double, 2> coord) -> double {
          double x = coord[0];
          double y = coord[1];
          return 1 + ((x * 0.3) + y);
        };
        auto g = [&](const std::array<double, 2> coord) -> double {
          double x = coord[0];
          double y = coord[0];
          return 99 + ((x * 7) + y * 0.1);
        };

        DomainTools::SetValues<2>(d, vec, f, g);
        DomainTools::SetValuesWithGhost<2>(d, expected, f, g);

        BiLinearGhostFiller blgf(d, GhostFillingType::Faces);
        blgf.fillGhost(vec);

        for (auto pinfo : d.getPatchInfoVector()) {
          ComponentView<double, 2> vec_ld = vec.getComponentView(0, pinfo.local_index);
          ComponentView<double, 2> expected_ld = expected.getComponentView(0, pinfo.local_index);
          ComponentView<double, 2> vec_ld2 = vec.getComponentView(1, pinfo.local_index);
          ComponentView<double, 2> expected_ld2 = expected.getComponentView(1, pinfo.local_index);
          Loop::Nested<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2>& coord) {
            ///
            REQUIRE_EQ(vec_ld[coord], Approx(expected_ld[coord]));
          });
          for (Side<2> s : Side<2>::getValues()) {
            View<double, 1> vec_ghost = vec_ld.getSliceOn(s, { -1 });
            View<double, 1> expected_ghost = expected_ld.getSliceOn(s, { -1 });
            if (pinfo.hasNbr(s)) {
              Loop::Nested<1>(vec_ghost.getStart(), vec_ghost.getEnd(), [&](const array<int, 1>& coord) {
                ///
                CHECK_EQ(vec_ghost[coord], Approx(expected_ghost[coord]));
              });
            }
          }
          Loop::Nested<2>(vec_ld2.getStart(), vec_ld2.getEnd(), [&](const array<int, 2>& coord) { REQUIRE_EQ(vec_ld2[coord], Approx(expected_ld2[coord])); });
          for (Side<2> s : Side<2>::getValues()) {
            View<double, 1> vec_ghost = vec_ld2.getSliceOn(s, { -1 });
            View<double, 1> expected_ghost = expected_ld2.getSliceOn(s, { -1 });
            if (pinfo.hasNbr(s)) {
              Loop::Nested<1>(vec_ghost.getStart(), vec_ghost.getEnd(), [&](const array<int, 1>& coord) { CHECK_EQ(vec_ghost[coord], Approx(expected_ghost[coord])); });
            }
          }
        }
      }
    }
  }
}
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller ghost already set two components")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;

        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d = domain_reader.getFinerDomain();

        Vector<2> vec(d, 2);
        Vector<2> expected(d, 2);

        auto f = [&](const std::array<double, 2> coord) -> double {
          double x = coord[0];
          double y = coord[0];
          return 1 + ((x * 0.3) + y);
        };
        auto g = [&](const std::array<double, 2> coord) -> double {
          double x = coord[0];
          double y = coord[0];
          return 99 + ((x * 7) + y * 0.1);
        };

        DomainTools::SetValuesWithGhost<2>(d, vec, f, g);
        DomainTools::SetValuesWithGhost<2>(d, expected, f, g);

        BiLinearGhostFiller blgf(d, GhostFillingType::Faces);
        blgf.fillGhost(vec);

        for (auto pinfo : d.getPatchInfoVector()) {
          ComponentView<double, 2> vec_ld = vec.getComponentView(0, pinfo.local_index);
          ComponentView<double, 2> expected_ld = expected.getComponentView(0, pinfo.local_index);
          ComponentView<double, 2> vec_ld2 = vec.getComponentView(1, pinfo.local_index);
          ComponentView<double, 2> expected_ld2 = expected.getComponentView(1, pinfo.local_index);
          Loop::Nested<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2>& coord) { REQUIRE_EQ(vec_ld[coord], Approx(expected_ld[coord])); });
          for (Side<2> s : Side<2>::getValues()) {
            View<double, 1> vec_ghost = vec_ld.getSliceOn(s, { -1 });
            View<double, 1> expected_ghost = expected_ld.getSliceOn(s, { -1 });
            if (pinfo.hasNbr(s)) {
              Loop::Nested<1>(vec_ghost.getStart(), vec_ghost.getEnd(), [&](const array<int, 1>& coord) { CHECK_EQ(vec_ghost[coord], Approx(expected_ghost[coord])); });
            }
          }
          Loop::Nested<2>(vec_ld2.getStart(), vec_ld2.getEnd(), [&](const array<int, 2>& coord) { REQUIRE_EQ(vec_ld2[coord], Approx(expected_ld2[coord])); });
          for (Side<2> s : Side<2>::getValues()) {
            View<double, 1> vec_ghost = vec_ld2.getSliceOn(s, { -1 });
            View<double, 1> expected_ghost = expected_ld2.getSliceOn(s, { -1 });
            if (pinfo.hasNbr(s)) {
              Loop::Nested<1>(vec_ghost.getStart(), vec_ghost.getEnd(), [&](const array<int, 1>& coord) { CHECK_EQ(vec_ghost[coord], Approx(expected_ghost[coord])); });
            }
          }
        }
      }
    }
  }
}
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller corners")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;

        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d = domain_reader.getFinerDomain();

        Vector<2> vec(d, 1);
        Vector<2> expected(d, 1);

        auto f = [&](const std::array<double, 2> coord) -> double {
          double x = coord[0];
          double y = coord[1];
          return 1 + ((x * 0.3) + y);
        };

        DomainTools::SetValues<2>(d, vec, f);
        DomainTools::SetValuesWithGhost<2>(d, expected, f);

        BiLinearGhostFiller blgf(d, GhostFillingType::Corners);
        blgf.fillGhost(vec);

        for (auto pinfo : d.getPatchInfoVector()) {
          ComponentView<double, 2> vec_ld = vec.getComponentView(0, pinfo.local_index);
          ComponentView<double, 2> expected_ld = expected.getComponentView(0, pinfo.local_index);
          Loop::Nested<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2>& coord) {
            ///
            REQUIRE_EQ(vec_ld[coord], Approx(expected_ld[coord]));
          });
          for (Side<2> s : Side<2>::getValues()) {
            View<double, 1> vec_ghost = vec_ld.getSliceOn(s, { -1 });
            View<double, 1> expected_ghost = expected_ld.getSliceOn(s, { -1 });
            if (pinfo.hasNbr(s)) {
              Loop::Nested<1>(vec_ghost.getStart(), vec_ghost.getEnd(), [&](const array<int, 1>& coord) {
                ///
                CHECK_EQ(vec_ghost[coord], Approx(expected_ghost[coord]));
              });
            }
          }
          for (Corner<2> c : Corner<2>::getValues()) {
            View<double, 0> vec_ghost = vec_ld.getSliceOn(c, { -1, -1 });
            View<double, 0> expected_ghost = expected_ld.getSliceOn(c, { -1, -1 });
            if (pinfo.hasNbr(c)) {
              CHECK_EQ(vec_ghost[{}], Approx(expected_ghost[{}]));
            }
          }
        }
      }
    }
  }
}
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller ghost already set corners")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;

        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d = domain_reader.getFinerDomain();

        Vector<2> vec(d, 1);
        Vector<2> expected(d, 1);

        auto f = [&](const std::array<double, 2> coord) -> double {
          double x = coord[0];
          double y = coord[0];
          return 1 + ((x * 0.3) + y);
        };

        DomainTools::SetValuesWithGhost<2>(d, vec, f);
        DomainTools::SetValuesWithGhost<2>(d, expected, f);

        BiLinearGhostFiller blgf(d, GhostFillingType::Corners);
        blgf.fillGhost(vec);

        for (auto pinfo : d.getPatchInfoVector()) {
          ComponentView<double, 2> vec_ld = vec.getComponentView(0, pinfo.local_index);
          ComponentView<double, 2> expected_ld = expected.getComponentView(0, pinfo.local_index);
          Loop::Nested<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2>& coord) {
            ///
            REQUIRE_EQ(vec_ld[coord], Approx(expected_ld[coord]));
          });
          for (Side<2> s : Side<2>::getValues()) {
            View<double, 1> vec_ghost = vec_ld.getSliceOn(s, { -1 });
            View<double, 1> expected_ghost = expected_ld.getSliceOn(s, { -1 });
            if (pinfo.hasNbr(s)) {
              Loop::Nested<1>(vec_ghost.getStart(), vec_ghost.getEnd(), [&](const array<int, 1>& coord) {
                ///
                CHECK_EQ(vec_ghost[coord], Approx(expected_ghost[coord]));
              });
            }
          }
          for (Corner<2> c : Corner<2>::getValues()) {
            View<double, 0> vec_ghost = vec_ld.getSliceOn(c, { -1, -1 });
            View<double, 0> expected_ghost = expected_ld.getSliceOn(c, { -1, -1 });
            if (pinfo.hasNbr(c)) {
              CHECK_EQ(vec_ghost[{}], Approx(expected_ghost[{}]));
            }
          }
        }
      }
    }
  }
}
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller two components corners")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;

        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d = domain_reader.getFinerDomain();

        Vector<2> vec(d, 2);
        Vector<2> expected(d, 2);

        auto f = [&](const std::array<double, 2> coord) -> double {
          double x = coord[0];
          double y = coord[1];
          return 1 + ((x * 0.3) + y);
        };
        auto g = [&](const std::array<double, 2> coord) -> double {
          double x = coord[0];
          double y = coord[0];
          return 99 + ((x * 7) + y * 0.1);
        };

        DomainTools::SetValues<2>(d, vec, f, g);
        DomainTools::SetValuesWithGhost<2>(d, expected, f, g);

        BiLinearGhostFiller blgf(d, GhostFillingType::Corners);
        blgf.fillGhost(vec);

        for (auto pinfo : d.getPatchInfoVector()) {
          ComponentView<double, 2> vec_ld = vec.getComponentView(0, pinfo.local_index);
          ComponentView<double, 2> expected_ld = expected.getComponentView(0, pinfo.local_index);
          ComponentView<double, 2> vec_ld2 = vec.getComponentView(1, pinfo.local_index);
          ComponentView<double, 2> expected_ld2 = expected.getComponentView(1, pinfo.local_index);
          Loop::Nested<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2>& coord) {
            ///
            REQUIRE_EQ(vec_ld[coord], Approx(expected_ld[coord]));
          });
          for (Side<2> s : Side<2>::getValues()) {
            View<double, 1> vec_ghost = vec_ld.getSliceOn(s, { -1 });
            View<double, 1> expected_ghost = expected_ld.getSliceOn(s, { -1 });
            if (pinfo.hasNbr(s)) {
              Loop::Nested<1>(vec_ghost.getStart(), vec_ghost.getEnd(), [&](const array<int, 1>& coord) {
                ///
                CHECK_EQ(vec_ghost[coord], Approx(expected_ghost[coord]));
              });
            }
          }
          for (Corner<2> c : Corner<2>::getValues()) {
            View<double, 0> vec_ghost = vec_ld.getSliceOn(c, { -1, -1 });
            View<double, 0> expected_ghost = expected_ld.getSliceOn(c, { -1, -1 });
            if (pinfo.hasNbr(c)) {
              CHECK_EQ(vec_ghost[{}], Approx(expected_ghost[{}]));
            }
          }
          Loop::Nested<2>(vec_ld2.getStart(), vec_ld2.getEnd(), [&](const array<int, 2>& coord) { REQUIRE_EQ(vec_ld2[coord], Approx(expected_ld2[coord])); });
          for (Side<2> s : Side<2>::getValues()) {
            View<double, 1> vec_ghost = vec_ld2.getSliceOn(s, { -1 });
            View<double, 1> expected_ghost = expected_ld2.getSliceOn(s, { -1 });
            if (pinfo.hasNbr(s)) {
              Loop::Nested<1>(vec_ghost.getStart(), vec_ghost.getEnd(), [&](const array<int, 1>& coord) { CHECK_EQ(vec_ghost[coord], Approx(expected_ghost[coord])); });
            }
          }
          for (Corner<2> c : Corner<2>::getValues()) {
            View<double, 0> vec_ghost = vec_ld2.getSliceOn(c, { -1, -1 });
            View<double, 0> expected_ghost = expected_ld2.getSliceOn(c, { -1, -1 });
            if (pinfo.hasNbr(c)) {
              CHECK_EQ(vec_ghost[{}], Approx(expected_ghost[{}]));
            }
          }
        }
      }
    }
  }
}
TEST_CASE("exchange various meshes 2D BiLinearGhostFiller ghost already set two components corners")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;

        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d = domain_reader.getFinerDomain();

        Vector<2> vec(d, 2);
        Vector<2> expected(d, 2);

        auto f = [&](const std::array<double, 2> coord) -> double {
          double x = coord[0];
          double y = coord[0];
          return 1 + ((x * 0.3) + y);
        };
        auto g = [&](const std::array<double, 2> coord) -> double {
          double x = coord[0];
          double y = coord[0];
          return 99 + ((x * 7) + y * 0.1);
        };

        DomainTools::SetValuesWithGhost<2>(d, vec, f, g);
        DomainTools::SetValuesWithGhost<2>(d, expected, f, g);

        BiLinearGhostFiller blgf(d, GhostFillingType::Corners);
        blgf.fillGhost(vec);

        for (auto pinfo : d.getPatchInfoVector()) {
          ComponentView<double, 2> vec_ld = vec.getComponentView(0, pinfo.local_index);
          ComponentView<double, 2> expected_ld = expected.getComponentView(0, pinfo.local_index);
          ComponentView<double, 2> vec_ld2 = vec.getComponentView(1, pinfo.local_index);
          ComponentView<double, 2> expected_ld2 = expected.getComponentView(1, pinfo.local_index);
          Loop::Nested<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2>& coord) { REQUIRE_EQ(vec_ld[coord], Approx(expected_ld[coord])); });
          for (Side<2> s : Side<2>::getValues()) {
            View<double, 1> vec_ghost = vec_ld.getSliceOn(s, { -1 });
            View<double, 1> expected_ghost = expected_ld.getSliceOn(s, { -1 });
            if (pinfo.hasNbr(s)) {
              Loop::Nested<1>(vec_ghost.getStart(), vec_ghost.getEnd(), [&](const array<int, 1>& coord) { CHECK_EQ(vec_ghost[coord], Approx(expected_ghost[coord])); });
            }
          }
          for (Corner<2> c : Corner<2>::getValues()) {
            View<double, 0> vec_ghost = vec_ld.getSliceOn(c, { -1, -1 });
            View<double, 0> expected_ghost = expected_ld.getSliceOn(c, { -1, -1 });
            if (pinfo.hasNbr(c)) {
              CHECK_EQ(vec_ghost[{}], Approx(expected_ghost[{}]));
            }
          }
          Loop::Nested<2>(vec_ld2.getStart(), vec_ld2.getEnd(), [&](const array<int, 2>& coord) { REQUIRE_EQ(vec_ld2[coord], Approx(expected_ld2[coord])); });
          for (Side<2> s : Side<2>::getValues()) {
            View<double, 1> vec_ghost = vec_ld2.getSliceOn(s, { -1 });
            View<double, 1> expected_ghost = expected_ld2.getSliceOn(s, { -1 });
            if (pinfo.hasNbr(s)) {
              Loop::Nested<1>(vec_ghost.getStart(), vec_ghost.getEnd(), [&](const array<int, 1>& coord) { CHECK_EQ(vec_ghost[coord], Approx(expected_ghost[coord])); });
            }
          }
          for (Corner<2> c : Corner<2>::getValues()) {
            View<double, 0> vec_ghost = vec_ld2.getSliceOn(c, { -1, -1 });
            View<double, 0> expected_ghost = expected_ld2.getSliceOn(c, { -1, -1 });
            if (pinfo.hasNbr(c)) {
              CHECK_EQ(vec_ghost[{}], Approx(expected_ghost[{}]));
            }
          }
        }
      }
    }
  }
}
