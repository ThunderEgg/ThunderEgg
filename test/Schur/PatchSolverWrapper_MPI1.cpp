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
#include "PatchSolverWrapper_MOCKS.h"
#include <ThunderEgg/Schur/PatchSolverWrapper.h>

#include <limits>

using namespace std;
using namespace ThunderEgg;
#define MESHES "mesh_inputs/2d_refined_east_1x2_mpi1.json", "mesh_inputs/2d_uniform_1x2_mpi1.json"
TEST_CASE("Schur::PatchSolverWrapper<2> throws exception for non-square patches")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 5, 7 }) {
      for (auto ny : { 6, 8 }) {
        DomainReader<2> domain_reader(mesh_file, { nx, ny }, 1);
        auto domain = domain_reader.getFinerDomain();
        Schur::InterfaceDomain<2> iface_domain(domain);
        MockGhostFiller<2> ghost_filler;
        MockPatchSolver<2> solver(domain, ghost_filler);

        CHECK_THROWS_AS(Schur::PatchSolverWrapper<2>(iface_domain, solver), RuntimeError);
      }
    }
  }
}
TEST_CASE("Schur::PatchSolverWrapper<2> apply fills ghost in rhs as expected")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 7 }) {
      for (double schur_fill_value : { 1.0, 1.3, 8.0, 2.0, -1.0 }) {
        DomainReader<2> domain_reader(mesh_file, { n, n }, 1);
        auto domain = domain_reader.getFinerDomain();
        Schur::InterfaceDomain<2> iface_domain(domain);
        MockGhostFiller<2> ghost_filler;
        RHSGhostCheckingPatchSolver<2> solver(domain, ghost_filler, schur_fill_value);

        Vector<1> x = iface_domain.getNewVector();
        Vector<1> b = iface_domain.getNewVector();

        x.set(schur_fill_value);

        // checking will be done in the solver
        Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
        psw.apply(x, b);
        CHECK_UNARY(solver.wasCalled());
      }
    }
  }
}
TEST_CASE("Schur::PatchSolverWrapper<2> apply gives expected rhs value for Schur matrix")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 7 }) {
      for (auto schur_fill_value : { 1.0, 1.3, 8.0, 2.0, -1.0 }) {
        for (auto domain_fill_value : { 1.0, 1.3, 8.0, 2.0, -1.0 }) {
          DomainReader<2> domain_reader(mesh_file, { n, n }, 1);
          auto domain = domain_reader.getFinerDomain();
          Schur::InterfaceDomain<2> iface_domain(domain);
          PatchFillingGhostFiller<2> ghost_filler(domain_fill_value);
          MockPatchSolver<2> solver(domain, ghost_filler);

          Vector<1> x = iface_domain.getNewVector();
          Vector<1> b = iface_domain.getNewVector();

          x.set(schur_fill_value);

          Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
          psw.apply(x, b);
          CHECK_UNARY(solver.allPatchesCalled());
          CHECK_UNARY(ghost_filler.wasCalled());
          for (int i = 0; i < b.getNumLocalPatches(); i++) {
            auto local_data = b.getComponentView(0, i);
            Loop::Nested<1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1>& coord) { CHECK_EQ(local_data[coord], doctest::Approx(schur_fill_value - domain_fill_value)); });
          }
        }
      }
    }
  }
}
TEST_CASE("Schur::PatchSolverWrapper<2> apply gives expected rhs value for Schur matrix with rhs "
          "already set")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 7 }) {
      for (auto schur_fill_value : { 1.0, 1.3, 8.0, 2.0, -1.0 }) {
        for (auto domain_fill_value : { 1.0, 1.3, 8.0, 2.0, -1.0 }) {
          DomainReader<2> domain_reader(mesh_file, { n, n }, 1);
          auto domain = domain_reader.getFinerDomain();
          Schur::InterfaceDomain<2> iface_domain(domain);
          PatchFillingGhostFiller<2> ghost_filler(domain_fill_value);
          MockPatchSolver<2> solver(domain, ghost_filler);

          Vector<1> x = iface_domain.getNewVector();
          Vector<1> b = iface_domain.getNewVector();

          x.set(schur_fill_value);
          b.set(99);

          Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
          psw.apply(x, b);
          CHECK_UNARY(solver.allPatchesCalled());
          CHECK_UNARY(ghost_filler.wasCalled());
          for (int i = 0; i < b.getNumLocalPatches(); i++) {
            auto local_data = b.getComponentView(0, i);
            Loop::Nested<1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1>& coord) { CHECK_EQ(local_data[coord], doctest::Approx(schur_fill_value - domain_fill_value)); });
          }
        }
      }
    }
  }
}
TEST_CASE("Schur::PatchSolverWrapper<2> getSchurRHSFromDomainRHS fills ghost in rhs as expected")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 7 }) {
      for (auto domain_fill_value : { 1.0, 1.3, 8.0, 2.0, -1.0 }) {
        DomainReader<2> domain_reader(mesh_file, { n, n }, 1);
        auto domain = domain_reader.getFinerDomain();
        Schur::InterfaceDomain<2> iface_domain(domain);
        MockGhostFiller<2> ghost_filler;
        RHSGhostCheckingPatchSolver<2> solver(domain, ghost_filler, 0);

        Vector<1> schur_b = iface_domain.getNewVector();
        Vector<2> domain_b(domain, 1);

        domain_b.set(domain_fill_value);

        // checking will be done in the solver
        Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
        psw.getSchurRHSFromDomainRHS(domain_b, schur_b);
        CHECK_UNARY(solver.wasCalled());
      }
    }
  }
}
TEST_CASE("Schur::PatchSolverWrapper<2> getSchurRHSFromDomainRHS gives expected rhs value for Schur matrix")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 7 }) {
      for (auto domain_fill_value : { 1.0, 1.3, 8.0, 2.0, -1.0 }) {
        DomainReader<2> domain_reader(mesh_file, { n, n }, 1);
        auto domain = domain_reader.getFinerDomain();
        Schur::InterfaceDomain<2> iface_domain(domain);
        PatchFillingGhostFiller<2> ghost_filler(domain_fill_value);
        MockPatchSolver<2> solver(domain, ghost_filler);

        Vector<1> schur_b = iface_domain.getNewVector();
        Vector<2> domain_b(domain, 1);

        Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
        psw.getSchurRHSFromDomainRHS(domain_b, schur_b);
        CHECK_UNARY(solver.allPatchesCalled());
        CHECK_UNARY(ghost_filler.wasCalled());
        for (int i = 0; i < schur_b.getNumLocalPatches(); i++) {
          auto local_data = schur_b.getComponentView(0, i);
          Loop::Nested<1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1>& coord) { CHECK_EQ(local_data[coord], doctest::Approx(domain_fill_value)); });
        }
      }
    }
  }
}
TEST_CASE("Schur::PatchSolverWrapper<2> getSchurRHSFromDomainRHS gives expected rhs value for "
          "Schur matrix with rhs already set")
{
  for (auto mesh_file : { MESHES }) {
    for (auto n : { 5, 7 }) {
      for (auto domain_fill_value : { 1.0, 1.3, 8.0, 2.0, -1.0 }) {
        DomainReader<2> domain_reader(mesh_file, { n, n }, 1);
        auto domain = domain_reader.getFinerDomain();
        Schur::InterfaceDomain<2> iface_domain(domain);
        PatchFillingGhostFiller<2> ghost_filler(domain_fill_value);
        MockPatchSolver<2> solver(domain, ghost_filler);

        Vector<1> schur_b = iface_domain.getNewVector();
        Vector<2> domain_b(domain, 1);

        schur_b.set(99);

        Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
        psw.getSchurRHSFromDomainRHS(domain_b, schur_b);
        CHECK_UNARY(solver.allPatchesCalled());
        CHECK_UNARY(ghost_filler.wasCalled());
        for (int i = 0; i < schur_b.getNumLocalPatches(); i++) {
          auto local_data = schur_b.getComponentView(0, i);
          Loop::Nested<1>(local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1>& coord) { CHECK_EQ(local_data[coord], doctest::Approx(domain_fill_value)); });
        }
      }
    }
  }
}
