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
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/LinearRestrictor.h>
#include <ThunderEgg/PETSc/MatWrapper.h>
#include <ThunderEgg/Poisson/MatrixHelper.h>
#include <ThunderEgg/Poisson/StarPatchOperator.h>
#include <ThunderEgg/TriLinearGhostFiller.h>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

#define MESHES "mesh_inputs/3d_uniform_2x2x2_mpi1.json", "mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", "mesh_inputs/3d_mid_refine_4x4x4_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";

TEST_CASE("Poisson::MatrixHelper gives equivalent operator to Poisson::StarPatchOperator")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 8, 10 }) {
      for (auto ny : { 8, 10 }) {
        for (auto nz : { 8, 10 }) {
          int num_ghost = 1;
          bitset<6> neumann;
          DomainReader<3> domain_reader(mesh_file, { nx, ny, nz }, num_ghost);
          Domain<3> d_fine = domain_reader.getFinerDomain();

          auto gfun = [](const std::array<double, 3>& coord) {
            double x = coord[0];
            double y = coord[1];
            double z = coord[2];
            return sin(M_PI * y) * cos(2 * M_PI * x) * cos(M_PI * z);
          };

          Vector<3> f_vec(d_fine, 1);
          Vector<3> f_vec_expected(d_fine, 1);

          Vector<3> g_vec(d_fine, 1);
          DomainTools::SetValues<3>(d_fine, g_vec, gfun);

          TriLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
          Poisson::StarPatchOperator<3> p_operator(d_fine, gf);
          p_operator.apply(g_vec, f_vec_expected);

          // generate matrix with matrix_helper
          Poisson::MatrixHelper mh(d_fine, neumann);
          Mat A = mh.formCRSMatrix();
          PETSc::MatWrapper<3> m_operator(A);
          m_operator.apply(g_vec, f_vec);

          REQUIRE_GT(f_vec.infNorm(), 0);

          for (auto pinfo : d_fine.getPatchInfoVector()) {
            ComponentView<double, 3> f_vec_ld = f_vec.getComponentView(0, pinfo.local_index);
            ComponentView<double, 3> f_vec_expected_ld = f_vec_expected.getComponentView(0, pinfo.local_index);
            Loop::Nested<3>(f_vec_ld.getStart(), f_vec_ld.getEnd(), [&](const array<int, 3>& coord) { CHECK_EQ(f_vec_ld[coord], doctest::Approx(f_vec_expected_ld[coord])); });
          }
          MatDestroy(&A);
        }
      }
    }
  }
}
TEST_CASE("Poisson::MatrixHelper gives equivalent operator to Poisson::StarPatchOperator with Neumann BC")
{
  for (auto mesh_file : { MESHES }) {
    for (auto nx : { 8, 10 }) {
      for (auto ny : { 8, 10 }) {
        for (auto nz : { 8, 10 }) {
          int num_ghost = 1;
          bitset<6> neumann = 0xFF;
          DomainReader<3> domain_reader(mesh_file, { nx, ny, nz }, num_ghost);
          Domain<3> d_fine = domain_reader.getFinerDomain();

          auto gfun = [](const std::array<double, 3>& coord) {
            double x = coord[0];
            double y = coord[1];
            double z = coord[2];
            return sin(M_PI * y) * cos(2 * M_PI * x) * cos(M_PI * z);
          };

          Vector<3> f_vec(d_fine, 1);
          Vector<3> f_vec_expected(d_fine, 1);

          Vector<3> g_vec(d_fine, 1);
          DomainTools::SetValues<3>(d_fine, g_vec, gfun);

          TriLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
          Poisson::StarPatchOperator<3> p_operator(d_fine, gf, true);
          p_operator.apply(g_vec, f_vec_expected);

          // generate matrix with matrix_helper
          Poisson::MatrixHelper mh(d_fine, neumann);
          Mat A = mh.formCRSMatrix();
          PETSc::MatWrapper<3> m_operator(A);
          m_operator.apply(g_vec, f_vec);

          REQUIRE_GT(f_vec.infNorm(), 0);

          for (auto pinfo : d_fine.getPatchInfoVector()) {
            ComponentView<double, 3> f_vec_ld = f_vec.getComponentView(0, pinfo.local_index);
            ComponentView<double, 3> f_vec_expected_ld = f_vec_expected.getComponentView(0, pinfo.local_index);
            Loop::Nested<3>(f_vec_ld.getStart(), f_vec_ld.getEnd(), [&](const array<int, 3>& coord) { CHECK_EQ(f_vec_ld[coord], doctest::Approx(f_vec_expected_ld[coord])); });
          }
          MatDestroy(&A);
        }
      }
    }
  }
}
TEST_CASE("Poisson::MatrixHelper constructor throws error with odd number of cells")
{
  for (auto mesh_file : { MESHES }) {
    for (auto axis : { 0, 1, 2 }) {
      int n_even = 10;
      int n_odd = 11;
      int num_ghost = 1;
      bitset<6> neumann;

      array<int, 3> ns;
      ns.fill(n_even);
      ns[axis] = n_odd;
      DomainReader<3> domain_reader(mesh_file, ns, num_ghost);
      Domain<3> d = domain_reader.getFinerDomain();

      CHECK_THROWS_AS(Poisson::MatrixHelper(d, neumann), RuntimeError);
    }
  }
}
