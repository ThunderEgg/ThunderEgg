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
#include <ThunderEgg/BiQuadraticGhostFiller.h>
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/LinearRestrictor.h>
#include <ThunderEgg/Poisson/StarPatchOperator.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

#define MESHES                                                                                     \
  "mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_4x4_mpi1.json",                  \
    "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json",                                             \
    "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";

TEST_CASE("Test Poisson::StarPatchOperator add ghost to RHS", "[Poisson::StarPatchOperator]")
{
  auto nx = GENERATE(2, 10);
  auto ny = GENERATE(2, 10);
  int num_ghost = 1;
  DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
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

  Vector<2> f_vec(d_fine, 1);
  DomainTools::SetValuesWithGhost<2>(d_fine, f_vec, ffun);

  Vector<2> g_vec(d_fine, 1);
  DomainTools::SetValuesWithGhost<2>(d_fine, g_vec, gfun);

  Vector<2> g_zero(d_fine, 1);
  DomainTools::SetValuesWithGhost<2>(d_fine, g_zero, gfun);
  g_zero.set(0);

  BiLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
  Poisson::StarPatchOperator<2> p_operator(d_fine, gf);
  p_operator.addDrichletBCToRHS(f_vec, gfun);

  Vector<2> f_expected = f_vec;
  for (auto pinfo : d_fine.getPatchInfoVector()) {
    auto u = g_vec.getComponentView(0, pinfo.local_index);
    auto f = f_expected.getComponentView(0, pinfo.local_index);
    for (Side<2> s : Side<2>::getValues()) {
      if (pinfo.hasNbr(s)) {
        double h2 = std::pow(pinfo.spacings[s.getAxisIndex()], 2);
        auto f_slice = f.getSliceOn(s, { 0 });
        auto u_inner = u.getSliceOn(s, { 0 });
        auto u_ghost = u.getSliceOn(s, { -1 });
        Loop::Nested<1>(f_slice.getStart(), f_slice.getEnd(), [&](std::array<int, 1> coord) {
          f_slice[coord] += -(u_inner[coord] + u_ghost[coord]) / (h2);
        });
      }
    }
  }

  for (auto pinfo : d_fine.getPatchInfoVector()) {
    auto gs = g_vec.getPatchView(pinfo.local_index);
    auto fs = f_vec.getPatchView(pinfo.local_index);
    p_operator.modifyRHSForInternalBoundaryConditions(pinfo, gs, fs);
  }

  for (auto pinfo : d_fine.getPatchInfoVector()) {
    INFO("Patch: " << pinfo.id);
    INFO("x:     " << pinfo.starts[0]);
    INFO("y:     " << pinfo.starts[1]);
    INFO("nx:    " << pinfo.ns[0]);
    INFO("ny:    " << pinfo.ns[1]);
    ComponentView<double, 2> vec_ld = f_vec.getComponentView(0, pinfo.local_index);
    ComponentView<double, 2> expected_ld = f_expected.getComponentView(0, pinfo.local_index);
    Loop::Nested<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2>& coord) {
      INFO("xi:    " << coord[0]);
      INFO("yi:    " << coord[1]);
      REQUIRE(vec_ld[coord] == Catch::Approx(expected_ld[coord]));
    });
  }
}
TEST_CASE("Test Poisson::StarPatchOperator apply on linear lhs constant coeff",
          "[Poisson::StarPatchOperator]")
{
  auto mesh_file = GENERATE(as<std::string>{}, MESHES);
  INFO("MESH FILE " << mesh_file);
  auto nx = GENERATE(2, 10);
  auto ny = GENERATE(2, 10);
  int num_ghost = 1;
  DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
  Domain<2> d_fine = domain_reader.getFinerDomain();

  auto gfun = [](const std::array<double, 2>& coord) {
    double x = coord[0];
    double y = coord[1];
    return 1 + 5 * x + 13 * y;
  };
  auto ffun = [](const std::array<double, 2>& coord) { return 0; };

  Vector<2> f_vec(d_fine, 1);
  Vector<2> f_vec_expected(d_fine, 1);
  DomainTools::SetValuesWithGhost<2>(d_fine, f_vec_expected, ffun);

  Vector<2> g_vec(d_fine, 1);
  DomainTools::SetValuesWithGhost<2>(d_fine, g_vec, gfun);

  BiLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
  Poisson::StarPatchOperator<2> p_operator(d_fine, gf);
  p_operator.addDrichletBCToRHS(f_vec_expected, gfun);

  p_operator.apply(g_vec, f_vec);

  for (auto pinfo : d_fine.getPatchInfoVector()) {
    INFO("Patch: " << pinfo.id);
    INFO("x:     " << pinfo.starts[0]);
    INFO("y:     " << pinfo.starts[1]);
    INFO("nx:    " << pinfo.ns[0]);
    INFO("ny:    " << pinfo.ns[1]);
    ComponentView<double, 2> vec_ld = f_vec.getComponentView(0, pinfo.local_index);
    ComponentView<double, 2> expected_vec_ld =
      f_vec_expected.getComponentView(0, pinfo.local_index);
    Loop::Nested<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2>& coord) {
      INFO("xi:    " << coord[0]);
      INFO("yi:    " << coord[1]);
      CHECK(vec_ld[coord] + 1 == Catch::Approx(1 + expected_vec_ld[coord]));
    });
  }
}
TEST_CASE("Test Poisson::StarPatchOperator apply on linear lhs constant coeff with Neumann BC",
          "[Poisson::StarPatchOperator]")
{
  auto mesh_file = GENERATE(as<std::string>{}, MESHES);
  INFO("MESH FILE " << mesh_file);
  auto nx = GENERATE(2, 10);
  auto ny = GENERATE(2, 10);
  int num_ghost = 1;
  DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
  Domain<2> d_fine = domain_reader.getFinerDomain();

  auto gfun = [](const std::array<double, 2>& coord) {
    double x = coord[0];
    double y = coord[1];
    return 1 + 5 * x + 13 * y;
  };
  auto ffun = [](const std::array<double, 2>& coord) { return 0; };
  auto gfun_x = [](const std::array<double, 2>& coord) { return 5; };
  auto gfun_y = [](const std::array<double, 2>& coord) { return 13; };

  Vector<2> f_vec(d_fine, 1);
  Vector<2> f_vec_expected(d_fine, 1);
  DomainTools::SetValuesWithGhost<2>(d_fine, f_vec_expected, ffun);

  Vector<2> g_vec(d_fine, 1);
  DomainTools::SetValuesWithGhost<2>(d_fine, g_vec, gfun);

  BiLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
  Poisson::StarPatchOperator<2> p_operator(d_fine, gf, true);
  p_operator.addNeumannBCToRHS(f_vec_expected, gfun, { gfun_x, gfun_y });

  p_operator.apply(g_vec, f_vec);

  for (auto pinfo : d_fine.getPatchInfoVector()) {
    INFO("Patch: " << pinfo.id);
    INFO("x:     " << pinfo.starts[0]);
    INFO("y:     " << pinfo.starts[1]);
    INFO("nx:    " << pinfo.ns[0]);
    INFO("ny:    " << pinfo.ns[1]);
    ComponentView<double, 2> vec_ld = f_vec.getComponentView(0, pinfo.local_index);
    ComponentView<double, 2> exptected_vec_ld =
      f_vec_expected.getComponentView(0, pinfo.local_index);
    Loop::Nested<2>(vec_ld.getStart(), vec_ld.getEnd(), [&](const array<int, 2>& coord) {
      INFO("xi:    " << coord[0]);
      INFO("yi:    " << coord[1]);
      CHECK(vec_ld[coord] + 1 == Catch::Approx(1 + exptected_vec_ld[coord]));
    });
  }
}
TEST_CASE("Test Poisson::StarPatchOperator gets 2nd order convergence",
          "[Poisson::StarPatchOperator]")
{
  auto mesh_file = GENERATE(as<std::string>{}, MESHES);
  INFO("MESH FILE " << mesh_file);
  int ns[2] = { 32, 64 };
  int num_ghost = 1;
  double errors[2];
  for (int i = 0; i < 2; i++) {
    DomainReader<2> domain_reader(mesh_file, { ns[i], ns[i] }, num_ghost);
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
    auto hfun = [](const std::array<double, 2>& coord) { return 1; };

    Vector<2> f_vec(d_fine, 1);
    Vector<2> f_vec_expected(d_fine, 1);
    DomainTools::SetValues<2>(d_fine, f_vec_expected, ffun);

    Vector<2> g_vec(d_fine, 1);
    DomainTools::SetValues<2>(d_fine, g_vec, gfun);

    Vector<2> h_vec(d_fine, 1);
    DomainTools::SetValuesWithGhost<2>(d_fine, h_vec, hfun);

    BiLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
    Poisson::StarPatchOperator<2> p_operator(d_fine, gf);
    p_operator.addDrichletBCToRHS(f_vec_expected, gfun);

    p_operator.apply(g_vec, f_vec);

    Vector<2> error_vec(d_fine, 1);
    error_vec.addScaled(1.0, f_vec, -1.0, f_vec_expected);
    errors[i] = error_vec.twoNorm() / f_vec_expected.twoNorm();
  }
  INFO("Errors: " << errors[0] << ", " << errors[1]);
  CHECK(log(errors[0] / errors[1]) / log(2) > 1.8);
}
TEST_CASE("Test Poisson::StarPatchOperator gets 2nd order convergence with neumann BC",
          "[Poisson::StarPatchOperator]")
{
  auto mesh_file = GENERATE(as<std::string>{}, MESHES);
  INFO("MESH FILE " << mesh_file);
  int ns[2] = { 32, 64 };
  int num_ghost = 1;
  double errors[2];
  for (int i = 0; i < 2; i++) {
    DomainReader<2> domain_reader(mesh_file, { ns[i], ns[i] }, num_ghost);
    Domain<2> d_fine = domain_reader.getFinerDomain();

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

    Vector<2> f_vec(d_fine, 1);
    Vector<2> f_vec_expected(d_fine, 1);
    DomainTools::SetValues<2>(d_fine, f_vec_expected, ffun);

    Vector<2> g_vec(d_fine, 1);
    DomainTools::SetValues<2>(d_fine, g_vec, gfun);

    BiQuadraticGhostFiller gf(d_fine, GhostFillingType::Faces);
    Poisson::StarPatchOperator<2> p_operator(d_fine, gf, true);
    p_operator.addNeumannBCToRHS(f_vec_expected, gfun, { gfun_x, gfun_y });

    p_operator.apply(g_vec, f_vec);

    Vector<2> error_vec(d_fine, 1);
    error_vec.addScaled(1.0, f_vec, -1.0, f_vec_expected);
    errors[i] = error_vec.twoNorm() / f_vec_expected.twoNorm();
  }
  INFO("Errors: " << errors[0] << ", " << errors[1]);
  CHECK(log(errors[0] / errors[1]) / log(2) > 1.8);
}
TEST_CASE("Test Poisson::StarPatchOperator constructor throws exception with no ghost cells",
          "[Poisson::StarPatchOperator]")
{
  auto mesh_file = GENERATE(as<std::string>{}, MESHES);
  INFO("MESH FILE " << mesh_file);
  int n = 32;
  int num_ghost = 0;

  DomainReader<2> domain_reader(mesh_file, { n, n }, num_ghost);
  Domain<2> d_fine = domain_reader.getFinerDomain();

  BiLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
  CHECK_THROWS_AS(Poisson::StarPatchOperator<2>(d_fine, gf), ThunderEgg::RuntimeError);
}