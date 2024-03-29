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
#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/LinearRestrictor.h>
#include <ThunderEgg/Poisson/DFTPatchSolver.h>
#include <ThunderEgg/Poisson/StarPatchOperator.h>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

#define MESHES "mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";

TEST_CASE("Test Poisson::DFTPatchSolver gets 2nd order convergence")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 10, 13 }) {
      for (auto ny : { 10, 13 }) {
        int num_ghost = 1;
        bitset<4> neumann;
        double errors[2];
        for (int i = 1; i <= 2; i++) {
          DomainReader<2> domain_reader(mesh_file, { i * nx, i * ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();

          auto ffun = [](const std::array<double, 2>& coord) {
            double x = coord[0];
            double y = coord[1];
            return -5 * M_PI * M_PI * sinl(M_PI * y) * cosl(2 * M_PI * x);
          };
          auto gfun = [](const std::array<double, 2>& coord) {
            double x = coord[0];
            double y = coord[1];
            return sinl(M_PI * y) * cosl(2 * M_PI * x);
          };

          Vector<2> g_vec(d_fine, 1);
          DomainTools::SetValuesWithGhost<2>(d_fine, g_vec, gfun);
          Vector<2> g_vec_expected(d_fine, 1);
          DomainTools::SetValues<2>(d_fine, g_vec_expected, gfun);

          Vector<2> f_vec(d_fine, 1);
          DomainTools::SetValues<2>(d_fine, f_vec, ffun);

          BiLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
          Poisson::StarPatchOperator<2> p_operator(d_fine, gf);
          Poisson::DFTPatchSolver<2> p_solver(p_operator, neumann);
          p_operator.addDrichletBCToRHS(f_vec, gfun);

          p_solver.smooth(f_vec, g_vec);

          Vector<2> error_vec(d_fine, 1);
          error_vec.addScaled(1.0, g_vec, -1.0, g_vec_expected);
          errors[i - 1] = error_vec.twoNorm() / g_vec_expected.twoNorm();
        }
        CHECK_GT(log(errors[0] / errors[1]) / log(2), 1.8);
      }
    }
  }
}
TEST_CASE("Test Poisson::DFTPatchSolver gets 2nd order convergence with neumann boundary")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 10, 13 }) {
      for (auto ny : { 10, 13 }) {
        int num_ghost = 1;
        bitset<4> neumann = 0xF;
        double errors[2];
        for (int i = 1; i <= 2; i++) {
          DomainReader<2> domain_reader(mesh_file, { i * nx, i * ny }, num_ghost);
          Domain<2> d_fine = domain_reader.getFinerDomain();

          auto ffun = [](const std::array<double, 2>& coord) {
            double x = coord[0];
            double y = coord[1];
            return -5 * M_PI * M_PI * sinl(M_PI * y) * cosl(2 * M_PI * x);
          };
          auto gfun = [](const std::array<double, 2>& coord) {
            double x = coord[0];
            double y = coord[1];
            return sinl(M_PI * y) * cosl(2 * M_PI * x);
          };
          auto gfun_x = [](const std::array<double, 2>& coord) {
            double x = coord[0];
            double y = coord[1];
            return -2 * M_PI * sin(M_PI * y) * sin(2 * M_PI * x);
          };
          auto gfun_y = [](const std::array<double, 2>& coord) {
            double x = coord[0];
            double y = coord[1];
            return M_PI * cos(M_PI * y) * cos(2 * M_PI * x);
          };

          Vector<2> g_vec(d_fine, 1);
          DomainTools::SetValuesWithGhost<2>(d_fine, g_vec, gfun);
          Vector<2> g_vec_expected(d_fine, 1);
          DomainTools::SetValues<2>(d_fine, g_vec_expected, gfun);

          Vector<2> f_vec(d_fine, 1);
          DomainTools::SetValues<2>(d_fine, f_vec, ffun);

          BiLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
          Poisson::StarPatchOperator<2> p_operator(d_fine, gf, true);
          Poisson::DFTPatchSolver<2> p_solver(p_operator, neumann);
          p_operator.addNeumannBCToRHS(f_vec, gfun, { gfun_x, gfun_y });

          p_solver.smooth(f_vec, g_vec);

          Vector<2> error_vec(d_fine, 1);
          g_vec.shift(-DomainTools::Integrate<2>(d_fine, g_vec) / d_fine.volume());
          g_vec_expected.shift(-DomainTools::Integrate<2>(d_fine, g_vec_expected) / d_fine.volume());
          error_vec.addScaled(1.0, g_vec, -1.0, g_vec_expected);
          errors[i - 1] = error_vec.twoNorm() / g_vec_expected.twoNorm();
        }
        CHECK_GT(log(errors[0] / errors[1]) / log(2), 1.8);
      }
    }
  }
}
TEST_CASE("Test Poisson::DFTPatchSolver gets 2nd order convergence with neumann boundary single patch")
{
  for (auto nx : { 10, 13 }) {
    for (auto ny : { 10, 13 }) {
      auto mesh_file = "mesh_inputs/2d_uniform_2x2_mpi1.json";
      int num_ghost = 1;
      bitset<4> neumann = 0xF;
      double errors[2];
      for (int i = 1; i <= 2; i++) {
        DomainReader<2> domain_reader(mesh_file, { i * nx, i * ny }, num_ghost);
        Domain<2> d_fine = domain_reader.getCoarserDomain();

        auto ffun = [](const std::array<double, 2>& coord) {
          double x = coord[0];
          double y = coord[1];
          return -5 * M_PI * M_PI * sinl(M_PI * y) * cosl(2 * M_PI * x);
        };
        auto gfun = [](const std::array<double, 2>& coord) {
          double x = coord[0];
          double y = coord[1];
          return sinl(M_PI * y) * cosl(2 * M_PI * x);
        };
        auto gfun_x = [](const std::array<double, 2>& coord) {
          double x = coord[0];
          double y = coord[1];
          return -2 * M_PI * sin(M_PI * y) * sin(2 * M_PI * x);
        };
        auto gfun_y = [](const std::array<double, 2>& coord) {
          double x = coord[0];
          double y = coord[1];
          return M_PI * cos(M_PI * y) * cos(2 * M_PI * x);
        };

        Vector<2> g_vec(d_fine, 1);
        DomainTools::SetValuesWithGhost<2>(d_fine, g_vec, gfun);
        Vector<2> g_vec_expected(d_fine, 1);
        DomainTools::SetValues<2>(d_fine, g_vec_expected, gfun);

        Vector<2> f_vec(d_fine, 1);
        DomainTools::SetValues<2>(d_fine, f_vec, ffun);

        BiLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
        Poisson::StarPatchOperator<2> p_operator(d_fine, gf, true);
        Poisson::DFTPatchSolver<2> p_solver(p_operator, neumann);
        p_operator.addNeumannBCToRHS(f_vec, gfun, { gfun_x, gfun_y });

        p_solver.smooth(f_vec, g_vec);

        Vector<2> error_vec(d_fine, 1);
        g_vec.shift(-DomainTools::Integrate<2>(d_fine, g_vec) / d_fine.volume());
        g_vec_expected.shift(-DomainTools::Integrate<2>(d_fine, g_vec_expected) / d_fine.volume());
        error_vec.addScaled(1.0, g_vec, -1.0, g_vec_expected);
        errors[i - 1] = error_vec.twoNorm() / g_vec_expected.twoNorm();
      }
      CHECK_GT(log(errors[0] / errors[1]) / log(2), 1.8);
    }
  }
}
