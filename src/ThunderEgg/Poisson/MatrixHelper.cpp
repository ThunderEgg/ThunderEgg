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
using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::Poisson;
MatrixHelper::MatrixHelper(const Domain<3> &domain, std::bitset<6> neumann) : domain(domain), neumann(neumann)
{
	for (int n : domain.getNs()) {
		if (n % 2 != 0) {
			throw RuntimeError("TriLinearGhostFiller only supports even number patch sizes");
		}
	}
}
namespace
{
/**
 * @brief Helps with adding stencil values at boundary patches
 */
class StencilHelper
{
	public:
	/**
	 * @brief the number of x cells on the patch boundary
	 */
	int nx;
	/**
	 * @brief the number of y cells on the patch boundary
	 */
	int ny;
	/**
	 * @brief Destroy the StencilHelper object
	 */
	virtual ~StencilHelper() {}
	/**
	 * @brief get the matrix row for a given coordinate on the boundary
	 *
	 * @param xi the x coordinate
	 * @param yi the y coordinate
	 * @return int the row of the matrix
	 */
	virtual int row(int xi, int yi) = 0;
	/**
	 * @brief get the number of coefficients for a given coordinate on the boundary
	 *
	 * @param xi the x coordinate
	 * @param yi the y coordinate
	 * @return int the number of coefficients
	 */
	virtual int size(int xi, int yi) = 0;
	/**
	 * @brief get the coefficients for a given coordinate on the boundary
	 *
	 * @param xi the x coordinate
	 * @param yi the y coordinate
	 * @return int the coefficients
	 */
	virtual const double *coeffs(int xi, int yi) = 0;
	/**
	 * @brief get the column indexes for the coefficients for a given coordinate on the boundary
	 *
	 * @param xi the x coordinate
	 * @param yi the y coordinate
	 * @return int the indexes
	 */
	virtual const int *cols(int xi, int yi) = 0;
};
/**
 * @brief Modifies stencil value for Dirichlet boundary conditions
 */
class DirichletSH : public StencilHelper
{
	private:
	double coeff;
	int    col;
	int    start;
	int    stridex;
	int    stridey;

	public:
	/**
	 * @brief Construct a new DirichletSH object
	 *
	 * @param pinfo the domain
	 * @param s the side of the boundary we are processing
	 */
	DirichletSH(const PatchInfo<3> &pinfo, Side<3> s)
	{
		double h   = 0;
		int    idx = pinfo.global_index * pinfo.ns[0] * pinfo.ns[1] * pinfo.ns[2];
		if (s == Side<3>::west()) {
			h       = pinfo.spacings[0];
			start   = idx;
			stridex = pinfo.ns[0];
			stridey = pinfo.ns[0] * pinfo.ns[1];
			nx      = pinfo.ns[1];
			ny      = pinfo.ns[2];
		} else if (s == Side<3>::east()) {
			h       = pinfo.spacings[0];
			start   = idx + (pinfo.ns[0] - 1);
			stridex = pinfo.ns[0];
			stridey = pinfo.ns[0] * pinfo.ns[1];
			nx      = pinfo.ns[1];
			ny      = pinfo.ns[2];
		} else if (s == Side<3>::south()) {
			h       = pinfo.spacings[1];
			start   = idx;
			stridex = 1;
			stridey = pinfo.ns[0] * pinfo.ns[1];
			nx      = pinfo.ns[0];
			ny      = pinfo.ns[2];
		} else if (s == Side<3>::north()) {
			h       = pinfo.spacings[1];
			start   = idx + (pinfo.ns[1] - 1) * pinfo.ns[0];
			stridex = 1;
			stridey = pinfo.ns[0] * pinfo.ns[1];
			nx      = pinfo.ns[0];
			ny      = pinfo.ns[2];
		} else if (s == Side<3>::bottom()) {
			h       = pinfo.spacings[2];
			start   = idx;
			stridex = 1;
			stridey = pinfo.ns[0];
			nx      = pinfo.ns[0];
			ny      = pinfo.ns[1];
		} else if (s == Side<3>::top()) {
			h       = pinfo.spacings[2];
			start   = idx + (pinfo.ns[2] - 1) * pinfo.ns[1] * pinfo.ns[0];
			stridex = 1;
			stridey = pinfo.ns[0];
			nx      = pinfo.ns[0];
			ny      = pinfo.ns[1];
		}
		coeff = -1.0 / (h * h);
	}
	int row(int xi, int yi)
	{
		return start + stridex * xi + stridey * yi;
	}
	int size(int xi, int yi)
	{
		return 1;
	}
	const double *coeffs(int xi, int yi)
	{
		return &coeff;
	}
	const int *cols(int xi, int yi)
	{
		col = start + stridex * xi + stridey * yi;
		return &col;
	}
};
/**
 * @brief Modifies stencil value for Neumann boundary conditions
 */
class NeumannSH : public StencilHelper
{
	private:
	double coeff;
	int    col;
	int    start;
	int    stridex;
	int    stridey;

	public:
	/**
	 * @brief Construct a new NeumannSH object
	 *
	 * @param pinfo the domain
	 * @param s the side of the boundary we are processing
	 */
	NeumannSH(const PatchInfo<3> &pinfo, Side<3> s)
	{
		double h   = 0;
		int    idx = pinfo.global_index * pinfo.ns[0] * pinfo.ns[1] * pinfo.ns[2];
		if (s == Side<3>::west()) {
			h       = pinfo.spacings[0];
			start   = idx;
			stridex = pinfo.ns[0];
			stridey = pinfo.ns[0] * pinfo.ns[1];
			nx      = pinfo.ns[1];
			ny      = pinfo.ns[2];
		} else if (s == Side<3>::east()) {
			h       = pinfo.spacings[0];
			start   = idx + (pinfo.ns[0] - 1);
			stridex = pinfo.ns[0];
			stridey = pinfo.ns[0] * pinfo.ns[1];
			nx      = pinfo.ns[1];
			ny      = pinfo.ns[2];
		} else if (s == Side<3>::south()) {
			h       = pinfo.spacings[1];
			start   = idx;
			stridex = 1;
			stridey = pinfo.ns[0] * pinfo.ns[1];
			nx      = pinfo.ns[0];
			ny      = pinfo.ns[2];
		} else if (s == Side<3>::north()) {
			h       = pinfo.spacings[1];
			start   = idx + (pinfo.ns[1] - 1) * pinfo.ns[0];
			stridex = 1;
			stridey = pinfo.ns[0] * pinfo.ns[1];
			nx      = pinfo.ns[0];
			ny      = pinfo.ns[2];
		} else if (s == Side<3>::bottom()) {
			h       = pinfo.spacings[2];
			start   = idx;
			stridex = 1;
			stridey = pinfo.ns[0];
			nx      = pinfo.ns[0];
			ny      = pinfo.ns[1];
		} else if (s == Side<3>::top()) {
			h       = pinfo.spacings[2];
			start   = idx + (pinfo.ns[2] - 1) * pinfo.ns[1] * pinfo.ns[0];
			stridex = 1;
			stridey = pinfo.ns[0];
			nx      = pinfo.ns[0];
			ny      = pinfo.ns[1];
		}
		coeff = 1.0 / (h * h);
	}
	int row(int xi, int yi)
	{
		return start + stridex * xi + stridey * yi;
	}
	int size(int xi, int yi)
	{
		return 1;
	}
	const double *coeffs(int xi, int yi)
	{
		return &coeff;
	}
	const int *cols(int xi, int yi)
	{
		col = start + stridex * xi + stridey * yi;
		return &col;
	}
};
/**
 * @brief Adds stencil for ghost values that are on a normal neighbor
 */
class NormalSH : public StencilHelper
{
	private:
	double coeff;
	int    col;
	int    start;
	int    nbr_start;
	int    stridex;
	int    stridey;

	public:
	/**
	 * @brief Construct a new NormalSH object
	 *
	 * @param pinfo the domain
	 * @param s the side of the boundary we are processing
	 */
	NormalSH(const PatchInfo<3> &pinfo, Side<3> s)
	{
		const NormalNbrInfo<2> &nbr_info = pinfo.getNormalNbrInfo(s);
		double                  h        = 0;
		int                     idx      = pinfo.global_index * pinfo.ns[0] * pinfo.ns[1] * pinfo.ns[2];
		int                     nbr_idx  = nbr_info.global_index * pinfo.ns[0] * pinfo.ns[1] * pinfo.ns[2];
		if (s == Side<3>::west()) {
			h         = pinfo.spacings[0];
			start     = idx;
			nbr_start = nbr_idx + (pinfo.ns[0] - 1);
			stridex   = pinfo.ns[0];
			stridey   = pinfo.ns[0] * pinfo.ns[1];
			nx        = pinfo.ns[1];
			ny        = pinfo.ns[2];
		} else if (s == Side<3>::east()) {
			h         = pinfo.spacings[0];
			start     = idx + (pinfo.ns[0] - 1);
			nbr_start = nbr_idx;
			stridex   = pinfo.ns[0];
			stridey   = pinfo.ns[0] * pinfo.ns[1];
			nx        = pinfo.ns[1];
			ny        = pinfo.ns[2];
		} else if (s == Side<3>::south()) {
			h         = pinfo.spacings[1];
			start     = idx;
			nbr_start = nbr_idx + (pinfo.ns[1] - 1) * pinfo.ns[0];
			stridex   = 1;
			stridey   = pinfo.ns[0] * pinfo.ns[1];
			nx        = pinfo.ns[0];
			ny        = pinfo.ns[2];
		} else if (s == Side<3>::north()) {
			h         = pinfo.spacings[1];
			start     = idx + (pinfo.ns[1] - 1) * pinfo.ns[0];
			nbr_start = nbr_idx;
			stridex   = 1;
			stridey   = pinfo.ns[0] * pinfo.ns[1];
			nx        = pinfo.ns[0];
			ny        = pinfo.ns[2];
		} else if (s == Side<3>::bottom()) {
			h         = pinfo.spacings[2];
			start     = idx;
			nbr_start = nbr_idx + (pinfo.ns[2] - 1) * pinfo.ns[1] * pinfo.ns[0];
			stridex   = 1;
			stridey   = pinfo.ns[0];
			nx        = pinfo.ns[0];
			ny        = pinfo.ns[1];
		} else if (s == Side<3>::top()) {
			h         = pinfo.spacings[2];
			start     = idx + (pinfo.ns[2] - 1) * pinfo.ns[1] * pinfo.ns[0];
			nbr_start = nbr_idx;
			stridex   = 1;
			stridey   = pinfo.ns[0];
			nx        = pinfo.ns[0];
			ny        = pinfo.ns[1];
		}
		coeff = 1.0 / (h * h);
	}
	int row(int xi, int yi)
	{
		return start + stridex * xi + stridey * yi;
	}
	int size(int xi, int yi)
	{
		return 1;
	}
	const double *coeffs(int xi, int yi)
	{
		return &coeff;
	}
	const int *cols(int xi, int yi)
	{
		col = nbr_start + stridex * xi + stridey * yi;
		return &col;
	}
};
/**
 * @brief Adds stencil for ghost values that are on a coarse neighbor
 */
class CoarseSH : public StencilHelper
{
	private:
	std::valarray<double> coeff = {{5.0 / 6, -1.0 / 6, -1.0 / 6, -1.0 / 6, 4.0 / 6}};
	int                   colz[5];
	int                   quad;
	int                   start;
	int                   nbr_start;
	int                   stridex;
	int                   stridey;

	public:
	/**
	 * @brief Construct a new CoarseSH object
	 *
	 * @param pinfo the domain
	 * @param s the side of the boundary we are processing
	 */
	CoarseSH(const PatchInfo<3> &pinfo, Side<3> s)
	{
		const CoarseNbrInfo<2> &nbr_info = pinfo.getCoarseNbrInfo(s);
		double                  h        = 0;
		int                     idx      = pinfo.global_index * pinfo.ns[0] * pinfo.ns[1] * pinfo.ns[2];
		int                     nbr_idx  = nbr_info.global_index * pinfo.ns[0] * pinfo.ns[1] * pinfo.ns[2];
		quad                             = (int) nbr_info.orth_on_coarse.getIndex();
		if (s == Side<3>::west()) {
			h         = pinfo.spacings[0];
			start     = idx;
			nbr_start = nbr_idx + (pinfo.ns[0] - 1);
			stridex   = pinfo.ns[0];
			stridey   = pinfo.ns[0] * pinfo.ns[1];
			nx        = pinfo.ns[1];
			ny        = pinfo.ns[2];
		} else if (s == Side<3>::east()) {
			h         = pinfo.spacings[0];
			start     = idx + (pinfo.ns[0] - 1);
			nbr_start = nbr_idx;
			stridex   = pinfo.ns[0];
			stridey   = pinfo.ns[0] * pinfo.ns[1];
			nx        = pinfo.ns[1];
			ny        = pinfo.ns[2];
		} else if (s == Side<3>::south()) {
			h         = pinfo.spacings[1];
			start     = idx;
			nbr_start = nbr_idx + (pinfo.ns[1] - 1) * pinfo.ns[0];
			stridex   = 1;
			stridey   = pinfo.ns[0] * pinfo.ns[1];
			nx        = pinfo.ns[0];
			ny        = pinfo.ns[2];
		} else if (s == Side<3>::north()) {
			h         = pinfo.spacings[1];
			start     = idx + (pinfo.ns[1] - 1) * pinfo.ns[0];
			nbr_start = nbr_idx;
			stridex   = 1;
			stridey   = pinfo.ns[0] * pinfo.ns[1];
			nx        = pinfo.ns[0];
			ny        = pinfo.ns[2];
		} else if (s == Side<3>::bottom()) {
			h         = pinfo.spacings[2];
			start     = idx;
			nbr_start = nbr_idx + (pinfo.ns[2] - 1) * pinfo.ns[1] * pinfo.ns[0];
			stridex   = 1;
			stridey   = pinfo.ns[0];
			nx        = pinfo.ns[0];
			ny        = pinfo.ns[1];
		} else if (s == Side<3>::top()) {
			h         = pinfo.spacings[2];
			start     = idx + (pinfo.ns[2] - 1) * pinfo.ns[1] * pinfo.ns[0];
			nbr_start = nbr_idx;
			stridex   = 1;
			stridey   = pinfo.ns[0];
			nx        = pinfo.ns[0];
			ny        = pinfo.ns[1];
		}
		coeff /= (h * h);
	}
	int row(int xi, int yi)
	{
		return start + stridex * xi + stridey * yi;
	}
	int size(int xi, int yi)
	{
		return (int) coeff.size();
	}
	const double *coeffs(int xi, int yi)
	{
		return &coeff[0];
	}
	const int *cols(int xi, int yi)
	{
		switch (quad) {
			case 0:
				colz[4] = nbr_start + stridex * (xi / 2) + stridey * (yi / 2);
				break;
			case 1:
				colz[4] = nbr_start + stridex * ((xi + nx) / 2) + stridey * (yi / 2);
				break;
			case 2:
				colz[4] = nbr_start + stridex * (xi / 2) + stridey * ((yi + ny) / 2);
				break;
			case 3:
				colz[4] = nbr_start + stridex * ((xi + nx) / 2) + stridey * ((yi + ny) / 2);
				break;
			default:
				break;
		}
		int nxi;
		if (xi % 2 == 0) {
			nxi = xi + 1;
		} else {
			nxi = xi - 1;
		}
		int nyi = 0;
		if (yi % 2 == 0) {
			nyi = yi + 1;
		} else {
			nyi = yi - 1;
		}
		colz[0] = start + stridex * xi + stridey * yi;
		colz[1] = start + stridex * nxi + stridey * yi;
		colz[2] = start + stridex * xi + stridey * nyi;
		colz[3] = start + stridex * nxi + stridey * nyi;
		return colz;
	}
};
class FineSH : public StencilHelper
{
	private:
	std::valarray<double> coeff = {{-1.0 / 3, 1.0 / 3, 1.0 / 3, 1.0 / 3, 1.0 / 3}};
	int                   colz[5];
	int                   start;
	int                   nbr_start[4];
	int                   stridex;
	int                   stridey;

	public:
	/**
	 * @brief Construct a new FineSH object
	 *
	 * @param pinfo the domain
	 * @param s the side of the boundary we are processing
	 */
	FineSH(const PatchInfo<3> &pinfo, Side<3> s)
	{
		const FineNbrInfo<2> &nbr_info = pinfo.getFineNbrInfo(s);
		double                h        = 0;
		int                   idx      = pinfo.global_index * pinfo.ns[0] * pinfo.ns[1] * pinfo.ns[2];
		int                   nbr_idx[4];
		for (int i = 0; i < 4; i++) {
			nbr_idx[i] = nbr_info.global_indexes[i] * pinfo.ns[0] * pinfo.ns[1] * pinfo.ns[2];
		}
		if (s == Side<3>::west()) {
			h     = pinfo.spacings[0];
			start = idx;
			for (int i = 0; i < 4; i++) {
				nbr_start[i] = nbr_idx[i] + (pinfo.ns[0] - 1);
			}
			stridex = pinfo.ns[0];
			stridey = pinfo.ns[0] * pinfo.ns[1];
			nx      = pinfo.ns[1];
			ny      = pinfo.ns[2];
		} else if (s == Side<3>::east()) {
			h     = pinfo.spacings[0];
			start = idx + (pinfo.ns[0] - 1);
			for (int i = 0; i < 4; i++) {
				nbr_start[i] = nbr_idx[i];
			}
			stridex = pinfo.ns[0];
			stridey = pinfo.ns[0] * pinfo.ns[1];
			nx      = pinfo.ns[1];
			ny      = pinfo.ns[2];
		} else if (s == Side<3>::south()) {
			h     = pinfo.spacings[1];
			start = idx;
			for (int i = 0; i < 4; i++) {
				nbr_start[i] = nbr_idx[i] + (pinfo.ns[1] - 1) * pinfo.ns[0];
			}
			stridex = 1;
			stridey = pinfo.ns[0] * pinfo.ns[1];
			nx      = pinfo.ns[0];
			ny      = pinfo.ns[2];
		} else if (s == Side<3>::north()) {
			h     = pinfo.spacings[1];
			start = idx + (pinfo.ns[1] - 1) * pinfo.ns[0];
			for (int i = 0; i < 4; i++) {
				nbr_start[i] = nbr_idx[i];
			}
			stridex = 1;
			stridey = pinfo.ns[0] * pinfo.ns[1];
			nx      = pinfo.ns[0];
			ny      = pinfo.ns[2];
		} else if (s == Side<3>::bottom()) {
			h     = pinfo.spacings[2];
			start = idx;
			for (int i = 0; i < 4; i++) {
				nbr_start[i] = nbr_idx[i] + (pinfo.ns[2] - 1) * pinfo.ns[1] * pinfo.ns[0];
			}
			stridex = 1;
			stridey = pinfo.ns[0];
			nx      = pinfo.ns[0];
			ny      = pinfo.ns[1];
		} else if (s == Side<3>::top()) {
			h     = pinfo.spacings[2];
			start = idx + (pinfo.ns[2] - 1) * pinfo.ns[1] * pinfo.ns[0];
			for (int i = 0; i < 4; i++) {
				nbr_start[i] = nbr_idx[i];
			}
			stridex = 1;
			stridex = 1;
			stridey = pinfo.ns[0];
			nx      = pinfo.ns[0];
			ny      = pinfo.ns[1];
		}
		coeff /= (h * h);
	}
	int row(int xi, int yi)
	{
		return start + stridex * xi + stridey * yi;
	}
	int size(int xi, int yi)
	{
		return (int) coeff.size();
	}
	const double *coeffs(int xi, int yi)
	{
		return &coeff[0];
	}
	const int *cols(int xi, int yi)
	{
		colz[0]  = start + stridex * xi + stridey * yi;
		int quad = (xi >= nx / 2) | ((yi >= ny / 2) << 1);
		int nxi  = xi % (nx / 2) * 2;
		int nyi  = yi % (ny / 2) * 2;
		colz[1]  = nbr_start[quad] + stridex * nxi + stridey * nyi;
		colz[2]  = nbr_start[quad] + stridex * (nxi + 1) + stridey * nyi;
		colz[3]  = nbr_start[quad] + stridex * nxi + stridey * (nyi + 1);
		colz[4]  = nbr_start[quad] + stridex * (nxi + 1) + stridey * (nyi + 1);
		return colz;
	}
};
/**
 * @brief Get the a StencilHelper object for a particular side
 *
 * @param pinfo the patch we are processing
 * @param s the side of the patch we are processing
 * @param neumann neumann boundary conditions
 * @return std::unique_ptr<StencilHelper>
 */
std::unique_ptr<StencilHelper> getStencilHelper(const PatchInfo<3> &pinfo, Side<3> s, std::bitset<6> neumann)
{
	StencilHelper *retval = nullptr;
	if (pinfo.hasNbr(s)) {
		switch (pinfo.getNbrType(s)) {
			case NbrType::Normal:
				retval = new NormalSH(pinfo, s);
				break;
			case NbrType::Fine:
				retval = new FineSH(pinfo, s);
				break;
			case NbrType::Coarse:
				retval = new CoarseSH(pinfo, s);
				break;
			default:
				throw RuntimeError("Unsupported NbrType");
				break;
		}
	} else {
		if (!pinfo.hasNbr(s) && neumann[s.getIndex()]) {
			retval = new NeumannSH(pinfo, s);
		} else {
			retval = new DirichletSH(pinfo, s);
		}
	}
	return std::unique_ptr<StencilHelper>(retval);
}
} // namespace
/**
 * @brief Add the the center coefficients to the Matrix
 *
 * @param A the matrix
 * @param pinfo the patch we are processing
 */
static void addCenterCoefficients(Mat A, const PatchInfo<3> &pinfo)
{
	int    nx    = pinfo.ns[0];
	int    ny    = pinfo.ns[1];
	int    nz    = pinfo.ns[2];
	double h_x   = pinfo.spacings[0];
	double h_y   = pinfo.spacings[1];
	double h_z   = pinfo.spacings[2];
	int    start = nx * ny * nz * pinfo.global_index;
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
static void addWestCoefficients(Mat A, const PatchInfo<3> &pinfo)
{
	int    nx    = pinfo.ns[0];
	int    ny    = pinfo.ns[1];
	int    nz    = pinfo.ns[2];
	double h_x   = pinfo.spacings[0];
	int    start = nx * ny * nz * pinfo.global_index;
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
static void addEastCoefficients(Mat A, const PatchInfo<3> &pinfo)
{
	int    nx    = pinfo.ns[0];
	int    ny    = pinfo.ns[1];
	int    nz    = pinfo.ns[2];
	double h_x   = pinfo.spacings[0];
	int    start = nx * ny * nz * pinfo.global_index;
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
static void addNorthCoefficients(Mat A, const PatchInfo<3> &pinfo)
{
	int    nx    = pinfo.ns[0];
	int    ny    = pinfo.ns[1];
	int    nz    = pinfo.ns[2];
	double h_y   = pinfo.spacings[1];
	int    start = nx * ny * nz * pinfo.global_index;
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
static void addSouthCoefficients(Mat A, const PatchInfo<3> &pinfo)
{
	int    nx    = pinfo.ns[0];
	int    ny    = pinfo.ns[1];
	int    nz    = pinfo.ns[2];
	double h_y   = pinfo.spacings[1];
	int    start = nx * ny * nz * pinfo.global_index;
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
static void addTopCoefficients(Mat A, const PatchInfo<3> &pinfo)
{
	int    nx    = pinfo.ns[0];
	int    ny    = pinfo.ns[1];
	int    nz    = pinfo.ns[2];
	double h_z   = pinfo.spacings[2];
	int    start = nx * ny * nz * pinfo.global_index;
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
static void addBottomCoefficients(Mat A, const PatchInfo<3> &pinfo)
{
	int    nx    = pinfo.ns[0];
	int    ny    = pinfo.ns[1];
	int    nz    = pinfo.ns[2];
	double h_z   = pinfo.spacings[2];
	int    start = nx * ny * nz * pinfo.global_index;
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
	int nx          = domain.getNs()[0];
	int ny          = domain.getNs()[1];
	int nz          = domain.getNs()[2];
	int local_size  = domain.getNumLocalPatches() * nx * ny * nz;
	int global_size = domain.getNumGlobalPatches() * nx * ny * nz;
	MatSetSizes(A, local_size, local_size, global_size, global_size);
	MatSetType(A, MATMPIAIJ);
	MatMPIAIJSetPreallocation(A, 19, nullptr, 19, nullptr);

	for (auto pinfo : domain.getPatchInfoVector()) {
		addCenterCoefficients(A, pinfo);
		addWestCoefficients(A, pinfo);
		addEastCoefficients(A, pinfo);
		addNorthCoefficients(A, pinfo);
		addSouthCoefficients(A, pinfo);
		addTopCoefficients(A, pinfo);
		addBottomCoefficients(A, pinfo);
		// boundaries
		for (Side<3> s : Side<3>::getValues()) {
			unique_ptr<StencilHelper> sh = getStencilHelper(pinfo, s, neumann);
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
