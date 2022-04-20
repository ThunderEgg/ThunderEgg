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
#include "PatchSolver_MOCKS.h"
#include "utils/DomainReader.h"
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/MPIGhostFiller.h>

#include <list>
#include <sstream>

using namespace std;
using namespace ThunderEgg;

constexpr auto single_mesh_file = "mesh_inputs/2d_uniform_2x2_mpi1.json";
constexpr auto refined_mesh_file = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json";
constexpr auto cross_mesh_file = "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json";

TEST_CASE("PatchSolver apply for various domains")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 5 }) {
      for (auto ny : { 2, 5 }) {
        for (auto u_num_components : { 1, 2, 3 }) {
          for (auto f_num_components : { 1, 2, 3 }) {
            int num_ghost = 1;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
            Domain<2> d_fine = domain_reader.getFinerDomain();

            Vector<2> u(d_fine, u_num_components);
            Vector<2> f(d_fine, f_num_components);

            MockGhostFiller<2> mgf;
            MockPatchSolver<2> mps(d_fine, mgf);

            u.setWithGhost(1);
            mps.apply(f, u);

            for (int i = 0; i < u.getNumLocalPatches(); i++) {
              for (int c = 0; c < u.getNumComponents(); c++) {
                auto ld = u.getComponentView(c, i);
                Loop::Nested<2>(ld.getStart(), ld.getEnd(), [&](const std::array<int, 2>& coord) { CHECK_EQ(ld[coord], 0); });
              }
            }
            CHECK_UNARY_FALSE(mgf.wasCalled());
            CHECK_UNARY(mps.allPatchesCalled());
          }
        }
      }
    }
  }
}
TEST_CASE("PatchSolver apply for various domains with timer")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 5 }) {
      for (auto ny : { 2, 5 }) {
        for (auto u_num_components : { 1, 2, 3 }) {
          for (auto f_num_components : { 1, 2, 3 }) {
            Communicator comm(MPI_COMM_WORLD);
            int num_ghost = 1;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
            Domain<2> d_fine = domain_reader.getFinerDomain();
            d_fine.setTimer(make_shared<Timer>(comm));

            Vector<2> u(d_fine, u_num_components);
            Vector<2> f(d_fine, f_num_components);

            MockGhostFiller<2> mgf;
            MockPatchSolver<2> mps(d_fine, mgf);

            u.setWithGhost(1);
            mps.apply(f, u);

            for (int i = 0; i < u.getNumLocalPatches(); i++) {
              for (int c = 0; c < u.getNumComponents(); c++) {
                auto ld = u.getComponentView(c, i);
                Loop::Nested<2>(ld.getStart(), ld.getEnd(), [&](const std::array<int, 2>& coord) { CHECK_EQ(ld[coord], 0); });
              }
            }
            CHECK_UNARY_FALSE(mgf.wasCalled());
            CHECK_UNARY(mps.allPatchesCalled());
            stringstream ss;
            ss << *d_fine.getTimer();
            CHECK_NE(ss.str().find("Total Patch"), string::npos);
            CHECK_NE(ss.str().find("Single Patch"), string::npos);
          }
        }
      }
    }
  }
}
TEST_CASE("PatchSolver smooth for various domains")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 5 }) {
      for (auto ny : { 2, 5 }) {
        for (auto u_num_components : { 1, 2, 3 }) {
          for (auto f_num_components : { 1, 2, 3 }) {
            Communicator comm(MPI_COMM_WORLD);
            int num_ghost = 1;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
            Domain<2> d_fine = domain_reader.getFinerDomain();
            d_fine.setTimer(make_shared<Timer>(comm));

            Vector<2> u(d_fine, u_num_components);
            Vector<2> f(d_fine, f_num_components);

            MockGhostFiller<2> mgf;
            MockPatchSolver<2> mps(d_fine, mgf);

            u.setWithGhost(1);
            mps.smooth(f, u);

            for (int i = 0; i < u.getNumLocalPatches(); i++) {
              for (int c = 0; c < u.getNumComponents(); c++) {
                auto ld = u.getComponentView(c, i);
                Loop::Nested<2>(ld.getStart(), ld.getEnd(), [&](const std::array<int, 2>& coord) { CHECK_EQ(ld[coord], 1); });
              }
            }
            CHECK_UNARY(mgf.wasCalled());
            CHECK_UNARY(mps.allPatchesCalled());
            stringstream ss;
            ss << *d_fine.getTimer();
            CHECK_NE(ss.str().find("Total Patch"), string::npos);
            CHECK_NE(ss.str().find("Single Patch"), string::npos);
          }
        }
      }
    }
  }
}
TEST_CASE("PatchSolver smooth for various domains with timer")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 5 }) {
      for (auto ny : { 2, 5 }) {
        for (auto u_num_components : { 1, 2, 3 }) {
          for (auto f_num_components : { 1, 2, 3 }) {
            int num_ghost = 1;
            DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
            Domain<2> d_fine = domain_reader.getFinerDomain();

            Vector<2> u(d_fine, u_num_components);
            Vector<2> f(d_fine, f_num_components);

            MockGhostFiller<2> mgf;
            MockPatchSolver<2> mps(d_fine, mgf);

            u.setWithGhost(1);
            mps.smooth(f, u);

            for (int i = 0; i < u.getNumLocalPatches(); i++) {
              for (int c = 0; c < u.getNumComponents(); c++) {
                auto ld = u.getComponentView(c, i);
                Loop::Nested<2>(ld.getStart(), ld.getEnd(), [&](const std::array<int, 2>& coord) { CHECK_EQ(ld[coord], 1); });
              }
            }
            CHECK_UNARY(mgf.wasCalled());
            CHECK_UNARY(mps.allPatchesCalled());
          }
        }
      }
    }
  }
}
TEST_CASE("PatchSolver getDomain")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 5 }) {
      for (auto ny : { 2, 5 }) {
        int num_ghost = 1;
        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d_fine = domain_reader.getFinerDomain();

        Vector<2> u(d_fine, 1);
        Vector<2> f(d_fine, 1);

        MockGhostFiller<2> mgf;
        MockPatchSolver<2> mps(d_fine, mgf);

        CHECK_EQ(mps.getDomain().getNumLocalPatches(), d_fine.getNumLocalPatches());
      }
    }
  }
}
TEST_CASE("PatchSolver getGhostFiller")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 5 }) {
      for (auto ny : { 2, 5 }) {
        int num_ghost = 1;
        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d_fine = domain_reader.getFinerDomain();

        Vector<2> u(d_fine, 1);
        Vector<2> f(d_fine, 1);

        MockGhostFiller<2> mgf;
        MockPatchSolver<2> mps(d_fine, mgf);

        const GhostFiller<2>& mps_mgf = mps.getGhostFiller();
        CHECK_EQ(typeid(mps_mgf), typeid(mgf));
      }
    }
  }
}
