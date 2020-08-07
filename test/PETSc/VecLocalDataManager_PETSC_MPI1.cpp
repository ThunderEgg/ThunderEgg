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
#include <ThunderEgg/PETSc/VecLocalDataManager.h>
using namespace std;
using namespace ThunderEgg;
#define MESHES                                                                                     \
	"mesh_inputs/2d_uniform_2x2_mpi1.json", "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json"
const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";
TEST_CASE("PETSc::VecLocalDataManager getVecView", "[PETSc::VecLocalDataManager]")
{
	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, 100, PETSC_DETERMINE, &vec);
	PETSc::VecLocalDataManager vldm(vec, false);

	double *view;
	VecGetArray(vec, &view);

	CHECK(vldm.getVecView() == view);

	VecRestoreArray(vec, &view);
	VecDestroy(&vec);
}
TEST_CASE("PETSc::VecLocalDataManager getVecView read only", "[PETSc::VecLocalDataManager]")
{
	Vec vec;
	VecCreateMPI(MPI_COMM_WORLD, 100, PETSC_DETERMINE, &vec);
	PETSc::VecLocalDataManager vldm(vec, true);

	double *view;
	VecGetArray(vec, &view);

	CHECK(vldm.getVecView() == view);

	VecRestoreArray(vec, &view);
	VecDestroy(&vec);
}