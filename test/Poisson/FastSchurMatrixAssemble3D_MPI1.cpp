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
#include <ThunderEgg/Poisson/FFTWPatchSolver.h>
#include <ThunderEgg/Poisson/FastSchurMatrixAssemble3D.h>
#include <ThunderEgg/Poisson/StarPatchOperator.h>
#include <ThunderEgg/Schur/PatchSolverWrapper.h>
#include <ThunderEgg/TriLinearGhostFiller.h>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

#define MESHES "mesh_inputs/3d_uniform_2x2x2_mpi1.json", "mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json", "mesh_inputs/3d_mid_refine_4x4x4_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";

TEST_CASE("Poisson::FastSchurMatrixAssemble3D throws exception for non-square patches")
{
  for (auto mesh_file : { MESHES }) {
    for (int nx : { 2, 8 }) {
      for (int ny : { 4, 10 }) {
        for (int nz : { 6, 12 }) {
          int num_ghost = 1;
          bitset<6> neumann;
          DomainReader<3> domain_reader(mesh_file, { nx, ny, nz }, num_ghost);
          Domain<3> d_fine = domain_reader.getFinerDomain();
          Schur::InterfaceDomain<3> iface_domain(d_fine);

          TriLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
          Poisson::StarPatchOperator<3> p_operator(d_fine, gf);
          Poisson::FFTWPatchSolver<3> p_solver(p_operator, neumann);

          CHECK_THROWS_AS(Poisson::FastSchurMatrixAssemble3D(iface_domain, p_solver), RuntimeError);
        }
      }
    }
  }
}
namespace ThunderEgg {
namespace {
template<int D>
class MockGhostFiller : public GhostFiller<D>
{
public:
  MockGhostFiller<D>* clone() const override { return new MockGhostFiller<D>(*this); }
  void fillGhost(const Vector<D>& u) const override {}
};
} // namespace
} // namespace ThunderEgg
TEST_CASE("Poisson::FastSchurMatrixAssemble3D throws with unsupported ghost filler")
{
  for (auto mesh_file : { MESHES }) {
    int num_ghost = 1;
    bitset<6> neumann;
    DomainReader<3> domain_reader(mesh_file, { 10, 10, 10 }, num_ghost);
    Domain<3> d_fine = domain_reader.getFinerDomain();
    Schur::InterfaceDomain<3> iface_domain(d_fine);

    MockGhostFiller<3> gf;
    Poisson::StarPatchOperator<3> p_operator(d_fine, gf);
    Poisson::FFTWPatchSolver<3> p_solver(p_operator, neumann);

    CHECK_THROWS_AS(Poisson::FastSchurMatrixAssemble3D(iface_domain, p_solver), RuntimeError);
  }
}
TEST_CASE("Poisson::FastSchurMatrixAssemble3D gives equivalent operator to "
          "Poisson::StarPatchOperator trilinear ghost filler")
{
  for (auto mesh_file : { MESHES }) {
    int n = 4;
    int num_ghost = 1;
    bitset<6> neumann;
    DomainReader<3> domain_reader(mesh_file, { n, n, n }, num_ghost);
    Domain<3> d_fine = domain_reader.getFinerDomain();
    Schur::InterfaceDomain<3> iface_domain(d_fine);

    Vector<2> f_vec = iface_domain.getNewVector();
    Vector<2> g_vec = iface_domain.getNewVector();
    Vector<2> f_vec_expected = iface_domain.getNewVector();

    int index = 0;
    for (auto iface_info : iface_domain.getInterfaces()) {
      View<double, 2> view = g_vec.getComponentView(0, iface_info->local_index);
      for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
          double x = (index + 0.5) / g_vec.getNumLocalCells();
          view(i, j) = sin(M_PI * x);
          index++;
        }
      }
    }

    TriLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
    Poisson::StarPatchOperator<3> p_operator(d_fine, gf);
    Poisson::FFTWPatchSolver<3> p_solver(p_operator, neumann);
    Schur::PatchSolverWrapper<3> p_solver_wrapper(iface_domain, p_solver);
    p_solver_wrapper.apply(g_vec, f_vec_expected);

    // generate matrix with matrix_helper
    Mat A = Poisson::FastSchurMatrixAssemble3D(iface_domain, p_solver);

    /*
    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "matrix.p", FILE_MODE_WRITE, &viewer);
    MatView(A, viewer);
    PetscViewerDestroy(&viewer);
    exit(0);
    */

    auto m_operator = make_shared<PETSc::MatWrapper<2>>(A);
    m_operator->apply(g_vec, f_vec);

    CHECK_EQ(f_vec.infNorm(), doctest::Approx(f_vec_expected.infNorm()));
    CHECK_EQ(f_vec.twoNorm(), doctest::Approx(f_vec_expected.twoNorm()));
    REQUIRE_GT(f_vec.infNorm(), 0);

    for (auto iface : iface_domain.getInterfaces()) {
      ComponentView<double, 2> f_vec_ld = f_vec.getComponentView(0, iface->local_index);
      ComponentView<double, 2> f_vec_expected_ld = f_vec_expected.getComponentView(0, iface->local_index);
      Loop::Nested<2>(f_vec_ld.getStart(), f_vec_ld.getEnd(), [&](const array<int, 2>& coord) { CHECK_EQ(f_vec_ld[coord], doctest::Approx(f_vec_expected_ld[coord])); });
    }
    MatDestroy(&A);
  }
}
TEST_CASE("Poisson::FastSchurMatrixAssemble3D gives equivalent operator to "
          "Poisson::StarPatchOperator with Neumann BC trilinear ghost filler")
{
  for (auto mesh_file : { MESHES }) {
    int n = 4;
    int num_ghost = 1;
    bitset<6> neumann = 0xFF;
    DomainReader<3> domain_reader(mesh_file, { n, n, n }, num_ghost);
    Domain<3> d_fine = domain_reader.getFinerDomain();
    Schur::InterfaceDomain<3> iface_domain(d_fine);

    Vector<2> f_vec = iface_domain.getNewVector();
    Vector<2> f_vec_expected = iface_domain.getNewVector();

    Vector<2> g_vec = iface_domain.getNewVector();

    int index = 0;
    for (auto iface_info : iface_domain.getInterfaces()) {
      View<double, 2> view = g_vec.getComponentView(0, iface_info->local_index);
      for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
          double x = (index + 0.5) / g_vec.getNumLocalCells();
          view(i, j) = sin(M_PI * x);
          index++;
        }
      }
    }

    TriLinearGhostFiller gf(d_fine, GhostFillingType::Faces);
    Poisson::StarPatchOperator<3> p_operator(d_fine, gf, true);
    Poisson::FFTWPatchSolver<3> p_solver(p_operator, neumann);
    Schur::PatchSolverWrapper<3> p_solver_wrapper(iface_domain, p_solver);
    p_solver_wrapper.apply(g_vec, f_vec_expected);

    // generate matrix with matrix_helper
    Mat A = Poisson::FastSchurMatrixAssemble3D(iface_domain, p_solver);
    auto m_operator = make_shared<PETSc::MatWrapper<2>>(A);
    m_operator->apply(g_vec, f_vec);

    CHECK_EQ(f_vec.infNorm(), doctest::Approx(f_vec_expected.infNorm()));
    CHECK_EQ(f_vec.twoNorm(), doctest::Approx(f_vec_expected.twoNorm()));
    REQUIRE_GT(f_vec.infNorm(), 0);

    for (int i = 0; i < f_vec.getNumLocalPatches(); i++) {
      ComponentView<double, 2> f_vec_ld = f_vec.getComponentView(0, i);
      ComponentView<double, 2> f_vec_expected_ld = f_vec_expected.getComponentView(0, i);
      Loop::Nested<2>(f_vec_ld.getStart(), f_vec_ld.getEnd(), [&](const array<int, 2>& coord) { CHECK_EQ(f_vec_ld[coord], doctest::Approx(f_vec_expected_ld[coord])); });
    }
    MatDestroy(&A);
  }
}
