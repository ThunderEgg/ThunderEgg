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

#include "../utils/DomainReader.h"
#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/Iterative/CG.h>
#include <ThunderEgg/Poisson/StarPatchOperator.h>

#include <sstream>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::Iterative;

TEST_CASE("CG default max iterations")
{
  CG<2> bcgs;
  CHECK_EQ(bcgs.getMaxIterations(), 1000);
}
TEST_CASE("CG set max iterations")
{
  for (int iterations : { 1, 2, 3 }) {
    CG<2> bcgs;
    bcgs.setMaxIterations(iterations);
    CHECK_EQ(bcgs.getMaxIterations(), iterations);
  }
}
TEST_CASE("CG default tolerance")
{
  CG<2> bcgs;
  CHECK_EQ(bcgs.getTolerance(), 1e-12);
}
TEST_CASE("CG set tolerance")
{
  for (double tolerance : { 1.2, 2.3, 3.4 }) {
    CG<2> bcgs;
    bcgs.setTolerance(tolerance);
    CHECK_EQ(bcgs.getTolerance(), tolerance);
  }
}
TEST_CASE("CG default timer")
{
  CG<2> bcgs;
  CHECK_EQ(bcgs.getTimer(), nullptr);
}
TEST_CASE("CG set timer")
{
  Communicator comm(MPI_COMM_WORLD);
  CG<2> bcgs;
  auto timer = make_shared<Timer>(comm);
  bcgs.setTimer(timer);
  CHECK_EQ(bcgs.getTimer(), timer);
}
TEST_CASE("CG clone")
{
  for (int iterations : { 1, 2, 3 }) {
    for (double tolerance : { 1.2, 2.3, 3.4 }) {
      CG<2> bcgs;
      bcgs.setMaxIterations(iterations);

      bcgs.setTolerance(tolerance);

      Communicator comm(MPI_COMM_WORLD);
      auto timer = make_shared<Timer>(comm);
      bcgs.setTimer(timer);

      unique_ptr<CG<2>> clone(bcgs.clone());
      CHECK_EQ(bcgs.getTimer(), clone->getTimer());
      CHECK_EQ(bcgs.getMaxIterations(), clone->getMaxIterations());
      CHECK_EQ(bcgs.getTolerance(), clone->getTolerance());
    }
  }
}
TEST_CASE("CG solves poisson problem within given tolerance")
{
  for (double tolerance : { 1e-9, 1e-7, 1e-5 }) {
    string mesh_file = "mesh_inputs/2d_uniform_2x2_mpi1.json";
    DomainReader<2> domain_reader(mesh_file, { 32, 32 }, 1);
    Domain<2> domain = domain_reader.getCoarserDomain();

    auto ffun = [](const std::array<double, 2>& coord) {
      double x = coord[0];
      double y = coord[1];
      return -5 * M_PI * M_PI * sin(M_PI * y) * cos(2 * M_PI * x);
    };
    auto gfun = [](const std::array<double, 2>& coord) {
      double x = coord[0];
      double y = coord[1];
      return sin(M_PI * y) * cos(2 * M_PI * x);
    };

    Vector<2> f_vec(domain, 1);
    DomainTools::SetValues<2>(domain, f_vec, ffun);
    Vector<2> residual(domain, 1);

    Vector<2> g_vec(domain, 1);

    BiLinearGhostFiller gf(domain, GhostFillingType::Faces);

    Poisson::StarPatchOperator<2> p_operator(domain, gf);
    p_operator.addDrichletBCToRHS(f_vec, gfun);

    CG<2> solver;
    solver.setMaxIterations(1000);
    solver.setTolerance(tolerance);
    solver.solve(p_operator, g_vec, f_vec);

    p_operator.apply(g_vec, residual);
    residual.addScaled(-1, f_vec);
    CHECK_LE(residual.dot(residual) / f_vec.dot(f_vec), tolerance);
  }
}
TEST_CASE("CG handles zero rhs vector")
{
  for (double tolerance : { 1e-9, 1e-7, 1e-5 }) {
    string mesh_file = "mesh_inputs/2d_uniform_2x2_mpi1.json";
    DomainReader<2> domain_reader(mesh_file, { 32, 32 }, 1);
    Domain<2> domain = domain_reader.getCoarserDomain();

    Vector<2> f_vec(domain, 1);

    Vector<2> g_vec(domain, 1);

    BiLinearGhostFiller gf(domain, GhostFillingType::Faces);

    Poisson::StarPatchOperator<2> p_operator(domain, gf);

    CG<2> solver;
    solver.setMaxIterations(1000);
    solver.setTolerance(tolerance);
    solver.solve(p_operator, g_vec, f_vec);

    CHECK_EQ(g_vec.infNorm(), 0);
  }
}
TEST_CASE("CG outputs iteration count and residual to output")
{
  string mesh_file = "mesh_inputs/2d_uniform_2x2_mpi1.json";
  DomainReader<2> domain_reader(mesh_file, { 32, 32 }, 1);
  Domain<2> domain = domain_reader.getCoarserDomain();

  auto ffun = [](const std::array<double, 2>& coord) {
    double x = coord[0];
    double y = coord[1];
    return -5 * M_PI * M_PI * sin(M_PI * y) * cos(2 * M_PI * x);
  };
  auto gfun = [](const std::array<double, 2>& coord) {
    double x = coord[0];
    double y = coord[1];
    return sin(M_PI * y) * cos(2 * M_PI * x);
  };

  Vector<2> f_vec(domain, 1);
  DomainTools::SetValues<2>(domain, f_vec, ffun);
  Vector<2> residual(domain, 1);

  Vector<2> g_vec(domain, 1);

  BiLinearGhostFiller gf(domain, GhostFillingType::Faces);

  Poisson::StarPatchOperator<2> p_operator(domain, gf);
  p_operator.addDrichletBCToRHS(f_vec, gfun);

  double tolerance = 1e-7;

  std::stringstream ss;

  CG<2> solver;
  solver.setMaxIterations(1000);
  solver.setTolerance(tolerance);
  solver.solve(p_operator, g_vec, f_vec, nullptr, true, ss);

  int prev_iteration;
  double resid;
  ss >> prev_iteration >> resid;
  while (prev_iteration < 18) {
    int iteration;
    ss >> iteration >> resid;
    CHECK_EQ(iteration, prev_iteration + 1);
    prev_iteration = iteration;
  }
}
TEST_CASE("CG giving a good initial guess reduces the iterations")
{
  string mesh_file = "mesh_inputs/2d_uniform_2x2_mpi1.json";
  DomainReader<2> domain_reader(mesh_file, { 32, 32 }, 1);
  Domain<2> domain = domain_reader.getCoarserDomain();

  auto ffun = [](const std::array<double, 2>& coord) {
    double x = coord[0];
    double y = coord[1];
    return -5 * M_PI * M_PI * sin(M_PI * y) * cos(2 * M_PI * x);
  };
  auto gfun = [](const std::array<double, 2>& coord) {
    double x = coord[0];
    double y = coord[1];
    return sin(M_PI * y) * cos(2 * M_PI * x);
  };

  Vector<2> f_vec(domain, 1);
  DomainTools::SetValues<2>(domain, f_vec, ffun);
  Vector<2> residual(domain, 1);

  Vector<2> g_vec(domain, 1);

  BiLinearGhostFiller gf(domain, GhostFillingType::Faces);

  Poisson::StarPatchOperator<2> p_operator(domain, gf);
  p_operator.addDrichletBCToRHS(f_vec, gfun);

  double tolerance = 1e-5;

  CG<2> solver;
  solver.setMaxIterations(1000);
  solver.setTolerance(tolerance);
  solver.solve(p_operator, g_vec, f_vec);

  int iterations_with_solved_guess = solver.solve(p_operator, g_vec, f_vec);

  CHECK_EQ(iterations_with_solved_guess, 0);
}
