/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
 *  top-level directory.
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
#include "catch.hpp"
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/LinearRestrictor.h>
#include <ThunderEgg/PETSc/MatWrapper.h>
#include <ThunderEgg/PETSc/VecWrapper.h>
#include <ThunderEgg/Poisson/FFTWPatchSolver.h>
#include <ThunderEgg/Poisson/FastSchurMatrixAssemble3D.h>
#include <ThunderEgg/Poisson/StarPatchOperator.h>
#include <ThunderEgg/Schur/PatchSolverWrapper.h>
#include <ThunderEgg/Schur/VecWrapperGenerator.h>
#include <ThunderEgg/TriLinearGhostFiller.h>
using namespace std;
using namespace ThunderEgg;
#define MESHES                                                                                     \
	"mesh_inputs/3d_uniform_2x2x2_mpi1.json", "mesh_inputs/3d_refined_bnw_2x2x2_mpi1.json",        \
	"mesh_inputs/3d_mid_refine_4x4x4_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";
TEST_CASE("Poisson::FastSchurMatrixAssemble3D throws exception for non-square patches",
          "[Poisson::FastSchurMatrixAssemble3D]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int                   nx        = GENERATE(2, 8);
	int                   ny        = GENERATE(4, 10);
	int                   nz        = GENERATE(6, 12);
	int                   num_ghost = 1;
	DomainReader<3>       domain_reader(mesh_file, {nx, ny, nz}, num_ghost);
	shared_ptr<Domain<3>> d_fine       = domain_reader.getFinerDomain();
	auto                  iface_domain = make_shared<Schur::InterfaceDomain<3>>(d_fine);

	auto gf         = make_shared<TriLinearGhostFiller>(d_fine);
	auto p_operator = make_shared<Poisson::StarPatchOperator<3>>(d_fine, gf);
	auto p_solver   = make_shared<Poisson::FFTWPatchSolver<3>>(p_operator);

	CHECK_THROWS_AS(Poisson::FastSchurMatrixAssemble3D(iface_domain, p_solver), RuntimeError);
}
namespace ThunderEgg
{
namespace
{
template <int D> class MockGhostFiller : public GhostFiller<D>
{
	public:
	void fillGhost(std::shared_ptr<const Vector<D>> u) const override {}
};
} // namespace
} // namespace ThunderEgg
TEST_CASE("Poisson::FastSchurMatrixAssemble3D throws with unsupported ghost filler",
          "[Poisson::FastSchurMatrixAssemble3D]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int                   num_ghost = 1;
	DomainReader<3>       domain_reader(mesh_file, {10, 10, 10}, num_ghost);
	shared_ptr<Domain<3>> d_fine       = domain_reader.getFinerDomain();
	auto                  iface_domain = make_shared<Schur::InterfaceDomain<3>>(d_fine);

	auto gf         = make_shared<MockGhostFiller<3>>();
	auto p_operator = make_shared<Poisson::StarPatchOperator<3>>(d_fine, gf);
	auto p_solver   = make_shared<Poisson::FFTWPatchSolver<3>>(p_operator);

	CHECK_THROWS_AS(Poisson::FastSchurMatrixAssemble3D(iface_domain, p_solver), RuntimeError);
}
TEST_CASE(
"Poisson::FastSchurMatrixAssemble3D gives equivalent operator to Poisson::StarPatchOperator trilinear ghost filler",
"[Poisson::FastSchurMatrixAssemble3D]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int                   n         = 4;
	int                   num_ghost = 1;
	DomainReader<3>       domain_reader(mesh_file, {n, n, n}, num_ghost);
	shared_ptr<Domain<3>> d_fine       = domain_reader.getFinerDomain();
	auto                  iface_domain = make_shared<Schur::InterfaceDomain<3>>(d_fine);

	Schur::VecWrapperGenerator<2> vg(iface_domain);
	auto                          f_vec          = vg.getNewVecWrapper();
	auto                          f_vec_expected = vg.getNewVecWrapper();

	auto    g_vec = vg.getNewVecWrapper();
	double *g_view;
	VecGetArray(g_vec->getVec(), &g_view);
	for (int i = 0; i < g_vec->getNumLocalCells(); i++) {
		double x  = (i + 0.5) / g_vec->getNumLocalCells();
		g_view[i] = sin(M_PI * x);
	}
	VecRestoreArray(g_vec->getVec(), &g_view);

	auto gf               = make_shared<TriLinearGhostFiller>(d_fine);
	auto p_operator       = make_shared<Poisson::StarPatchOperator<3>>(d_fine, gf);
	auto p_solver         = make_shared<Poisson::FFTWPatchSolver<3>>(p_operator);
	auto p_solver_wrapper = make_shared<Schur::PatchSolverWrapper<3>>(iface_domain, p_solver);
	p_solver_wrapper->apply(g_vec, f_vec_expected);

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

	CHECK(f_vec->infNorm() == Approx(f_vec_expected->infNorm()));
	CHECK(f_vec->twoNorm() == Approx(f_vec_expected->twoNorm()));
	REQUIRE(f_vec->infNorm() > 0);

	for (auto iface : iface_domain->getInterfaces()) {
		INFO("ID: " << iface->id);
		INFO("LOCAL_INDEX: " << iface->local_index);
		INFO("type: " << iface->patches.size());
		LocalData<2> f_vec_ld          = f_vec->getLocalData(0, iface->local_index);
		LocalData<2> f_vec_expected_ld = f_vec_expected->getLocalData(0, iface->local_index);
		nested_loop<2>(f_vec_ld.getStart(), f_vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			INFO("xi:    " << coord[0]);
			CHECK(f_vec_ld[coord] == Approx(f_vec_expected_ld[coord]));
		});
	}
	MatDestroy(&A);
}
TEST_CASE(
"Poisson::FastSchurMatrixAssemble3D gives equivalent operator to Poisson::StarPatchOperator with Neumann BC trilinear ghost filler",
"[Poisson::FastSchurMatrixAssemble3D]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int                   n         = 4;
	int                   num_ghost = 1;
	DomainReader<3>       domain_reader(mesh_file, {n, n, n}, num_ghost, true);
	shared_ptr<Domain<3>> d_fine       = domain_reader.getFinerDomain();
	auto                  iface_domain = make_shared<Schur::InterfaceDomain<3>>(d_fine);

	Schur::VecWrapperGenerator<2> vg(iface_domain);
	auto                          f_vec          = vg.getNewVecWrapper();
	auto                          f_vec_expected = vg.getNewVecWrapper();

	auto    g_vec = vg.getNewVecWrapper();
	double *g_view;
	Vec     g = g_vec->getVec();
	VecGetArray(g, &g_view);
	for (int i = 0; i < g_vec->getNumLocalCells(); i++) {
		double x  = (i + 0.5) / g_vec->getNumLocalCells();
		g_view[i] = sin(M_PI * x);
	}
	VecRestoreArray(g, &g_view);

	auto gf               = make_shared<TriLinearGhostFiller>(d_fine);
	auto p_operator       = make_shared<Poisson::StarPatchOperator<3>>(d_fine, gf, true);
	auto p_solver         = make_shared<Poisson::FFTWPatchSolver<3>>(p_operator);
	auto p_solver_wrapper = make_shared<Schur::PatchSolverWrapper<3>>(iface_domain, p_solver);
	p_solver_wrapper->apply(g_vec, f_vec_expected);

	// generate matrix with matrix_helper
	Mat  A          = Poisson::FastSchurMatrixAssemble3D(iface_domain, p_solver);
	auto m_operator = make_shared<PETSc::MatWrapper<2>>(A);
	m_operator->apply(g_vec, f_vec);

	CHECK(f_vec->infNorm() == Approx(f_vec_expected->infNorm()));
	CHECK(f_vec->twoNorm() == Approx(f_vec_expected->twoNorm()));
	REQUIRE(f_vec->infNorm() > 0);

	for (int i = 0; i < f_vec->getNumLocalPatches(); i++) {
		LocalData<2> f_vec_ld          = f_vec->getLocalData(0, i);
		LocalData<2> f_vec_expected_ld = f_vec_expected->getLocalData(0, i);
		nested_loop<2>(f_vec_ld.getStart(), f_vec_ld.getEnd(), [&](const array<int, 2> &coord) {
			INFO("xi:    " << coord[0]);
			CHECK(f_vec_ld[coord] == Approx(f_vec_expected_ld[coord]));
		});
	}
	MatDestroy(&A);
}