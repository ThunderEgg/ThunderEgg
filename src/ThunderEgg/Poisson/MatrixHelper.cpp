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

#include <ThunderEgg/Poisson/MatrixHelper.h>
#include <ThunderEgg/Poisson/StencilHelper.h>
using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::Poisson;
MatrixHelper::MatrixHelper(std::shared_ptr<Domain<3>> domain) : domain(domain)
{
	for (int n : domain->getNs()) {
		if (n % 2 != 0) {
			throw RuntimeError("TriLinearGhostFiller only supports even number patch sizes");
		}
	}
}
/**
 * @brief Add the the center coefficients to the Matrix
 *
 * @param A the matrix
 * @param pinfo the patch we are processing
 */
static void addCenterCoefficients(Mat A, std::shared_ptr<const PatchInfo<3>> pinfo)
{
	int    nx    = pinfo->ns[0];
	int    ny    = pinfo->ns[0];
	int    nz    = pinfo->ns[0];
	double h_x   = pinfo->spacings[0];
	double h_y   = pinfo->spacings[1];
	double h_z   = pinfo->spacings[2];
	int    start = nx * ny * nz * pinfo->global_index;
	// center coeffs
	double coeff = -2.0 / (h_x * h_x) - 2.0 / (h_y * h_y) - 2.0 / (h_z * h_z);
	for (int z_i = 0; z_i < nz; z_i++) {
		for (int y_i = 0; y_i < ny; y_i++) {
			for (int x_i = 0; x_i < nx; x_i++) {
				int row = start + x_i + nx * y_i + nx * ny * z_i;
				MatSetValues(A, 1, &row, 1, &row, &coeff, ADD_VALUES);
			}
		}
	}
}
/**
 * @brief Add the the west coefficients to the Matrix
 * will not add ghost coefficients
 *
 * @param A the matrix
 * @param pinfo the patch we are processing
 */
static void addWestCoefficients(Mat A, std::shared_ptr<const PatchInfo<3>> pinfo)
{
	int    nx    = pinfo->ns[0];
	int    ny    = pinfo->ns[0];
	int    nz    = pinfo->ns[0];
	double h_x   = pinfo->spacings[0];
	double h_y   = pinfo->spacings[1];
	double h_z   = pinfo->spacings[2];
	int    start = nx * ny * nz * pinfo->global_index;
	// west coeffs
	double coeff = 1.0 / (h_x * h_x);
	for (int z_i = 0; z_i < nz; z_i++) {
		for (int y_i = 0; y_i < ny; y_i++) {
			for (int x_i = 1; x_i < nx; x_i++) {
				int row = start + x_i + nx * y_i + nx * ny * z_i;
				int col = start + x_i - 1 + nx * y_i + nx * ny * z_i;
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
	}
}
/**
 * @brief Add the the east coefficients to the Matrix
 * will not add ghost coefficients
 *
 * @param A the matrix
 * @param pinfo the patch we are processing
 */
static void addEastCoefficients(Mat A, std::shared_ptr<const PatchInfo<3>> pinfo)
{
	int    nx    = pinfo->ns[0];
	int    ny    = pinfo->ns[0];
	int    nz    = pinfo->ns[0];
	double h_x   = pinfo->spacings[0];
	int    start = nx * ny * nz * pinfo->global_index;
	// east coeffs
	double coeff = 1.0 / (h_x * h_x);
	for (int z_i = 0; z_i < nz; z_i++) {
		for (int y_i = 0; y_i < ny; y_i++) {
			for (int x_i = 0; x_i < nx - 1; x_i++) {
				int row = start + x_i + nx * y_i + nx * ny * z_i;
				int col = start + x_i + 1 + nx * y_i + nx * ny * z_i;
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
	}
}
/**
 * @brief Add the the north coefficients to the Matrix
 * will not add ghost coefficients
 *
 * @param A the matrix
 * @param pinfo the patch we are processing
 */
static void addNorthCoefficients(Mat A, std::shared_ptr<const PatchInfo<3>> pinfo)
{
	int    nx    = pinfo->ns[0];
	int    ny    = pinfo->ns[0];
	int    nz    = pinfo->ns[0];
	double h_y   = pinfo->spacings[1];
	int    start = nx * ny * nz * pinfo->global_index;
	// north coeffs
	double coeff = 1.0 / (h_y * h_y);
	for (int z_i = 0; z_i < nz; z_i++) {
		for (int y_i = 0; y_i < ny - 1; y_i++) {
			for (int x_i = 0; x_i < nx; x_i++) {
				int row = start + x_i + nx * y_i + nx * ny * z_i;
				int col = start + x_i + nx * (y_i + 1) + nx * ny * z_i;
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
	}
}
/**
 * @brief Add the the south coefficients to the Matrix
 * will not add ghost coefficients
 *
 * @param A the matrix
 * @param pinfo the patch we are processing
 */
static void addSouthCoefficients(Mat A, std::shared_ptr<const PatchInfo<3>> pinfo)
{
	int    nx    = pinfo->ns[0];
	int    ny    = pinfo->ns[0];
	int    nz    = pinfo->ns[0];
	double h_y   = pinfo->spacings[1];
	int    start = nx * ny * nz * pinfo->global_index;
	// south coeffs
	double coeff = 1.0 / (h_y * h_y);
	for (int z_i = 0; z_i < nz; z_i++) {
		for (int y_i = 1; y_i < ny; y_i++) {
			for (int x_i = 0; x_i < nx; x_i++) {
				int row = start + x_i + nx * y_i + nx * ny * z_i;
				int col = start + x_i + nx * (y_i - 1) + nx * ny * z_i;
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
	}
}
/**
 * @brief Add the the top coefficients to the Matrix
 * will not add ghost coefficients
 *
 * @param A the matrix
 * @param pinfo the patch we are processing
 */
static void addTopCoefficients(Mat A, std::shared_ptr<const PatchInfo<3>> pinfo)
{
	int    nx    = pinfo->ns[0];
	int    ny    = pinfo->ns[0];
	int    nz    = pinfo->ns[0];
	double h_z   = pinfo->spacings[2];
	int    start = nx * ny * nz * pinfo->global_index;
	// top coeffs
	double coeff = 1.0 / (h_z * h_z);
	for (int z_i = 0; z_i < nz - 1; z_i++) {
		for (int y_i = 0; y_i < ny; y_i++) {
			for (int x_i = 0; x_i < nx; x_i++) {
				int row = start + x_i + nx * y_i + nx * ny * z_i;
				int col = start + x_i + nx * y_i + nx * ny * (z_i + 1);
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
	}
}
/**
 * @brief Add the the bottom coefficients to the Matrix
 * will not add ghost coefficients
 *
 * @param A the matrix
 * @param pinfo the patch we are processing
 */
static void addBottomCoefficients(Mat A, std::shared_ptr<const PatchInfo<3>> pinfo)
{
	int    nx    = pinfo->ns[0];
	int    ny    = pinfo->ns[0];
	int    nz    = pinfo->ns[0];
	double h_z   = pinfo->spacings[2];
	int    start = nx * ny * nz * pinfo->global_index;
	// top coeffs
	double coeff = 1.0 / (h_z * h_z);
	for (int z_i = 1; z_i < nz; z_i++) {
		for (int y_i = 0; y_i < ny; y_i++) {
			for (int x_i = 0; x_i < nx; x_i++) {
				int row = start + x_i + nx * y_i + nx * ny * z_i;
				int col = start + x_i + nx * y_i + nx * ny * (z_i - 1);
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
	}
}
Mat MatrixHelper::formCRSMatrix()
{
	Mat A;
	MatCreate(MPI_COMM_WORLD, &A);
	int nx          = domain->getNs()[0];
	int ny          = domain->getNs()[1];
	int nz          = domain->getNs()[2];
	int local_size  = domain->getNumLocalPatches() * nx * ny * nz;
	int global_size = domain->getNumGlobalPatches() * nx * ny * nz;
	MatSetSizes(A, local_size, local_size, global_size, global_size);
	MatSetType(A, MATMPIAIJ);
	MatMPIAIJSetPreallocation(A, 19, nullptr, 19, nullptr);

	for (auto pinfo : domain->getPatchInfoVector()) {
		addCenterCoefficients(A, pinfo);
		addWestCoefficients(A, pinfo);
		addEastCoefficients(A, pinfo);
		addNorthCoefficients(A, pinfo);
		addSouthCoefficients(A, pinfo);
		addTopCoefficients(A, pinfo);
		addBottomCoefficients(A, pinfo);
		// boundaries
		for (Side<3> s : Side<3>::getValues()) {
			unique_ptr<StencilHelper> sh = getStencilHelper(pinfo, s);
			for (int yi = 0; yi < sh->ny; yi++) {
				for (int xi = 0; xi < sh->nx; xi++) {
					int           row    = sh->row(xi, yi);
					int           size   = sh->size(xi, yi);
					const double *coeffs = sh->coeffs(xi, yi);
					const int *   cols   = sh->cols(xi, yi);
					MatSetValues(A, 1, &row, size, cols, coeffs, ADD_VALUES);
				}
			}
		}
	}
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	return A;
}
