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
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/PETSc/MatShellCreator.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

#define MESHES \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";

// mock operator
class HalfIdentity : public Operator<2>
{
	public:
	HalfIdentity *clone() const override
	{
		return new HalfIdentity(*this);
	}
	void apply(const Vector<2> &x, Vector<2> &b) const override
	{
		b.copy(x);
		b.scale(0.5);
	}
};
TEST_CASE("PETSc::MatShellCreator works with 0.5I", "[PETSc::MatShellCreator]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int             n         = 32;
	int             num_ghost = 0;
	DomainReader<2> domain_reader(mesh_file, {n, n}, num_ghost);
	Domain<2>       d_fine           = domain_reader.getFinerDomain();
	auto            vector_allocator = [&] {
        return Vector<2>(d_fine, 1);
	};

	Vec x;
	VecCreateMPI(MPI_COMM_WORLD, d_fine.getNumLocalCells() * 1, PETSC_DETERMINE, &x);
	double *x_view;
	VecGetArray(x, &x_view);
	for (int i = 0; i < d_fine.getNumLocalCells() * 1; i++) {
		x_view[i] = i;
	}
	VecRestoreArray(x, &x_view);
	Vec b;
	VecCreateMPI(MPI_COMM_WORLD, d_fine.getNumLocalCells() * 1, PETSC_DETERMINE, &b);

	// create an Identity matrix
	HalfIdentity TE_A;
	Mat          A = PETSc::MatShellCreator<2>::GetNewMatShell(TE_A, vector_allocator);

	MatMult(A, x, b);

	double *b_view;
	VecGetArray(x, &x_view);
	VecGetArray(b, &b_view);
	for (int i = 0; i < d_fine.getNumLocalCells() * 1; i++) {
		CHECK(x_view[i] * 0.5 == b_view[i]);
	}
	VecRestoreArray(x, &x_view);
	VecRestoreArray(b, &b_view);

	MatDestroy(&A);
	VecDestroy(&x);
	VecDestroy(&b);
}
TEST_CASE("PETSc::MatShellCreator works with 0.5I two components", "[PETSc::MatShellCreator]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH FILE " << mesh_file);
	int             n         = 32;
	int             num_ghost = 0;
	DomainReader<2> domain_reader(mesh_file, {n, n}, num_ghost);
	Domain<2>       d_fine           = domain_reader.getFinerDomain();
	auto            vector_allocator = [&] {
        return Vector<2>(d_fine, 2);
	};

	Vec x;
	VecCreateMPI(MPI_COMM_WORLD, d_fine.getNumLocalCells() * 2, PETSC_DETERMINE, &x);
	double *x_view;
	VecGetArray(x, &x_view);
	for (int i = 0; i < d_fine.getNumLocalCells() * 2; i++) {
		x_view[i] = i;
	}
	VecRestoreArray(x, &x_view);
	Vec b;
	VecCreateMPI(MPI_COMM_WORLD, d_fine.getNumLocalCells() * 2, PETSC_DETERMINE, &b);

	// create an Identity matrix
	HalfIdentity TE_A;
	Mat          A = PETSc::MatShellCreator<2>::GetNewMatShell(TE_A, vector_allocator);

	MatMult(A, x, b);

	double *b_view;
	VecGetArray(x, &x_view);
	VecGetArray(b, &b_view);
	for (int i = 0; i < d_fine.getNumLocalCells() * 2; i++) {
		CHECK(x_view[i] * 0.5 == b_view[i]);
	}
	VecRestoreArray(x, &x_view);
	VecRestoreArray(b, &b_view);

	MatDestroy(&A);
	VecDestroy(&x);
	VecDestroy(&b);
}