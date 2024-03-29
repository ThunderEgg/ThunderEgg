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
#include "../utils/DomainReader.h"
#include "PatchSolver_MOCKS.h"
#include <ThunderEgg/Iterative/PatchSolver.h>

#include <list>
#include <sstream>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

constexpr auto single_mesh_file = "mesh_inputs/2d_uniform_2x2_mpi1.json";
constexpr auto refined_mesh_file = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json";
constexpr auto cross_mesh_file = "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json";

TEST_CASE("Iterative::PatchSolver passes vectors of a single patch length")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 5 }) {
      for (auto ny : { 2, 5 }) {
        for (auto num_components : { 1, 2, 3 }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();

          Vector<2> u(d_fine, num_components);
          Vector<2> f(d_fine, num_components);

          MockGhostFiller<2> mgf;
          // the patch operator is just a 0.5I operator
          MockPatchOperator<2> mpo(d_fine, mgf);
          MockSolver<2> ms([](const Operator<2>& A, Vector<2>& x, const Vector<2>& b, const Operator<2>* Mr) {
            CHECK_EQ(x.getNumLocalPatches(), 1);
            return 1;
          });

          Iterative::PatchSolver<2> bcgs_solver(ms, mpo);

          bcgs_solver.smooth(f, u);
        }
      }
    }
  }
}
TEST_CASE("Iterative::PatchSolver passes modified operator")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 5 }) {
      for (auto ny : { 2, 5 }) {
        for (auto num_components : { 1, 2, 3 }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();

          Vector<2> u(d_fine, num_components);
          Vector<2> f(d_fine, num_components);

          bool called = false;
          MockGhostFiller<2> mgf;
          // the patch operator is just a 0.5I operator
          MockPatchOperator<2> mpo(d_fine, mgf);
          MockSolver<2> ms([&](const Operator<2>& A, Vector<2>& x, const Vector<2>& b, const Operator<2>* Mr) {
            if (!called) {
              called = true;
              A.apply(b, x);
            }
            return 1;
          });

          Iterative::PatchSolver<2> bcgs_solver(ms, mpo);

          bcgs_solver.smooth(f, u);
          CHECK_EQ(mpo.getNumApplyCalls(), 1);
          CHECK_UNARY(mpo.rhsWasModified());
          CHECK_UNARY(mpo.boundaryConditionsEnforced());
          CHECK_UNARY(mpo.internalBoundaryConditionsEnforced());
        }
      }
    }
  }
}
TEST_CASE("Iterative::PatchSolver propagates BreakdownError")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 5 }) {
      for (auto ny : { 2, 5 }) {
        for (auto num_components : { 1, 2, 3 }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();

          Vector<2> u(d_fine, num_components);
          Vector<2> f(d_fine, num_components);
          f.set(1);

          MockGhostFiller<2> mgf;
          // the patch operator is just a 0.5I operator
          NonLinMockPatchOperator<2> mpo(d_fine, mgf);
          MockSolver<2> ms([](const Operator<2>& A, Vector<2>& x, const Vector<2>& b, const Operator<2>* Mr) {
            throw Iterative::BreakdownError("Blah");
            return 1;
          });

          Iterative::PatchSolver<2> bcgs_solver(ms, mpo);

          CHECK_THROWS_AS(bcgs_solver.smooth(f, u), Iterative::BreakdownError);
        }
      }
    }
  }
}
TEST_CASE("Iterative::PatchSolver does not propagate BreakdownError")
{
  for (auto mesh_file : { single_mesh_file, refined_mesh_file, cross_mesh_file }) {
    for (auto nx : { 2, 5 }) {
      for (auto ny : { 2, 5 }) {
        for (auto num_components : { 1, 2, 3 }) {
          int num_ghost = 1;
          DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();

          Vector<2> u(d_fine, num_components);
          Vector<2> f(d_fine, num_components);
          f.set(1);

          MockGhostFiller<2> mgf;
          // the patch operator is just a 0.5I operator
          NonLinMockPatchOperator<2> mpo(d_fine, mgf);
          MockSolver<2> ms([](const Operator<2>& A, Vector<2>& x, const Vector<2>& b, const Operator<2>* Mr) {
            throw Iterative::BreakdownError("Blah");
            return 1;
          });

          Iterative::PatchSolver<2> bcgs_solver(ms, mpo, true);

          CHECK_NOTHROW(bcgs_solver.smooth(f, u));
        }
      }
    }
  }
}
