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

#ifndef THUNDEREGG_POISSON_DFTPATCHSOLVER_H
#define THUNDEREGG_POISSON_DFTPATCHSOLVER_H
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/PatchOperator.h>
#include <ThunderEgg/PatchSolver.h>
#include <ThunderEgg/ValVector.h>
#include <bitset>
#include <map>
#include <valarray>

extern "C" void dgemv_(char &, int &, int &, double &, double *, int &, double *, int &, double &,
                       double *, int &);

namespace ThunderEgg
{
namespace Poisson
{
/**
 * @brief This patch solver uses DFT transforms to solve for the Poisson equation
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> class DFTPatchSolver : public PatchSolver<D>
{
	private:
	/**
	 * @brief Enum of DFT types
	 */
	enum DftType { DCT_II, DCT_III, DCT_IV, DST_II, DST_III, DST_IV };
	/**
	 * @brief Comparator used in the maps, patches with the same spacings and boundary conditions
	 * will be equal
	 */
	struct CompareByBoundaryAndSpacings {
		bool operator()(const std::shared_ptr<const PatchInfo<D>> &a,
		                const std::shared_ptr<const PatchInfo<D>> &b) const
		{
			return std::forward_as_tuple(a->neumann.to_ulong(), a->spacings[0])
			       < std::forward_as_tuple(b->neumann.to_ulong(), b->spacings[0]);
		}
	};
	/**
	 * @brief The patch opertar that we are solving for
	 */
	std::shared_ptr<const PatchOperator<D>> op;
	/**
	 * @brief Map of patchinfo to DFT plan
	 */
	std::map<std::shared_ptr<const PatchInfo<D>>,
	         std::array<std::shared_ptr<std::valarray<double>>, D>, CompareByBoundaryAndSpacings>
	plan1;
	/**
	 * @brief Map of patchinfo to inverse DFT plan
	 */
	std::map<std::shared_ptr<const PatchInfo<D>>,
	         std::array<std::shared_ptr<std::valarray<double>>, D>, CompareByBoundaryAndSpacings>
	plan2;
	/**
	 * @brief Temporary copy for the modified right hand side
	 */
	std::shared_ptr<ValVector<D>> f_copy;
	/**
	 * @brief Temporary work vector
	 */
	std::shared_ptr<ValVector<D>> tmp;
	/*
	 * @brief Temporary work vector
	 */
	std::shared_ptr<ValVector<D>> local_tmp;
	/**
	 * @brief Map of PatchInfo object to it's respective eigenvalue array.
	 */
	std::map<std::shared_ptr<const PatchInfo<D>>, std::valarray<double>,
	         CompareByBoundaryAndSpacings>
	eigen_vals;
	/**
	 * @brief map of DFT transforms for each type and size
	 */
	std::map<std::tuple<DftType, int>, std::shared_ptr<std::valarray<double>>> transform_matrixes;

	/**
	 * @brief get arrays of coefficients necessary for each transform.
	 *
	 * @param types the tranforms for each axis
	 * @return std::array<std::shared_ptr<std::valarray<double>>, D> the coefficients for each
	 * transform
	 */
	std::array<std::shared_ptr<std::valarray<double>>, D> plan(std::array<DftType, D> types)
	{
		std::array<std::shared_ptr<std::valarray<double>>, D> retval;
		for (size_t axis = 0; axis < D; axis++) {
			retval[axis] = getTransformArray(types[axis], this->domain->getNs()[axis]);
		}
		return retval;
	}

	/**
	 * @brief Get a dft Transform matrix for a certain type and size
	 *
	 * @param type the DFT type
	 * @param n the size of the matrix
	 * @return std::shared_ptr<std::valarray<double>> an nxn DFT matrix
	 */
	std::shared_ptr<std::valarray<double>> getTransformArray(DftType type, int n)
	{
		std::shared_ptr<std::valarray<double>> matrix_ptr;

		if (transform_matrixes.count(std::make_tuple(type, n)) == 0) {
			matrix_ptr = std::make_shared<std::valarray<double>>(n * n);
			transform_matrixes[std::make_tuple(type, n)] = matrix_ptr;
			std::valarray<double> &matrix                = *matrix_ptr;

			switch (type) {
				case DftType::DCT_II:
					for (int j = 0; j < n; j++) {
						for (int i = 0; i < n; i++) {
							matrix[i * n + j] = cos(M_PI / n * (i * (j + 0.5)));
						}
					}
					break;
				case DftType::DCT_III:
					for (int i = 0; i < n; i++) {
						matrix[i * n] = 0.5;
					}
					for (int j = 1; j < n; j++) {
						for (int i = 0; i < n; i++) {
							matrix[i * n + j] = cos(M_PI / n * ((i + 0.5) * j));
						}
					}
					break;
				case DftType::DCT_IV:
					for (int j = 0; j < n; j++) {
						for (int i = 0; i < n; i++) {
							matrix[i * n + j] = cos(M_PI / n * ((i + 0.5) * (j + 0.5)));
						}
					}
					break;
				case DftType::DST_II:
					for (int i = 0; i < n; i++) {
						for (int j = 0; j < n; j++) {
							matrix[i * n + j] = sin(M_PI / n * ((i + 1) * (j + 0.5)));
						}
					}
					break;
				case DftType::DST_III:
					for (int i = 0; i < n; i += 2) {
						matrix[i * n + n - 1] = 0.5;
					}
					for (int i = 1; i < n; i += 2) {
						matrix[i * n + n - 1] = -0.5;
					}
					for (int i = 0; i < n; i++) {
						for (int j = 0; j < n - 1; j++) {
							matrix[i * n + j] = sin(M_PI / n * ((i + 0.5) * (j + 1)));
						}
					}
					break;
				case DftType::DST_IV:
					for (int j = 0; j < n; j++) {
						for (int i = 0; i < n; i++) {
							matrix[i * n + j] = sin(M_PI / n * ((i + 0.5) * (j + 0.5)));
						}
					}
			}
		} else {
			matrix_ptr = transform_matrixes.at(std::make_tuple(type, n));
		}
		return matrix_ptr;
	}
	/**
	 * @brief Execute a given DFT plan
	 *
	 * @param plan the plan (the matrixes for each axis)
	 * @param in the input values
	 * @param out the resulting values after the transform
	 * @param inverse weather we a re calculate
	 */
	void executePlan(const std::array<std::shared_ptr<std::valarray<double>>, D> &plan,
	                 LocalData<D> in, LocalData<D> out) const
	{
		LocalData<D> prev_result = in;

		for (size_t axis = 0; axis < D; axis++) {
			int                n     = this->domain->getNs()[axis];
			std::array<int, D> start = in.getStart();
			std::array<int, D> end   = in.getEnd();
			end[axis]                = 0;

			std::valarray<double> &matrix = *plan[axis];

			LocalData<D> new_result;
			if (D % 2) {
				if (axis % 2) {
					new_result = in;
				} else {
					new_result = out;
				}
			} else {
				if (axis == D - 1) {
					new_result = out;
				} else if (axis == D - 2) {
					new_result = local_tmp->getLocalData(0, 0);
				} else if (axis % 2) {
					new_result = in;
				} else {
					new_result = out;
				}
			}

			int pstride = prev_result.getStrides()[axis];
			int nstride = new_result.getStrides()[axis];

			char   T    = 'T';
			double one  = 1;
			double zero = 0;
			nested_loop<D>(start, end, [&](std::array<int, D> coord) {
				dgemv_(T, n, n, one, &matrix[0], n, &prev_result[coord], pstride, zero,
				       &new_result[coord], nstride);
			});

			prev_result = new_result;
		}
	}
	/**
	 * @brief Get the dft transform types for a patch
	 *
	 * @param pinfo the patch
	 * @return std::array<DftType, D> an array of tranforms for each axis
	 */
	std::array<DftType, D>
	getTransformsForPatch(std::shared_ptr<const ThunderEgg::PatchInfo<D>> pinfo)
	{
		// get transform types for each axis
		std::array<DftType, D> transforms;
		for (size_t axis = 0; axis < D; axis++) {
			if (pinfo->isNeumann(Side<D>::LowerSideOnAxis(axis))
			    && pinfo->isNeumann(Side<D>::HigherSideOnAxis(axis))) {
				transforms[axis] = DftType::DCT_II;
			} else if (pinfo->isNeumann(Side<D>::LowerSideOnAxis(axis))) {
				transforms[axis] = DftType::DCT_IV;
			} else if (pinfo->isNeumann(Side<D>::HigherSideOnAxis(axis))) {
				transforms[axis] = DftType::DST_IV;
			} else {
				transforms[axis] = DftType::DST_II;
			}
		}
		return transforms;
	}
	/**
	 * @brief Get the inverse dft transform types for a patch
	 *
	 * @param pinfo the patch
	 * @return std::array<DftType, D> an array of tranforms for each axis
	 */
	std::array<DftType, D>
	getInverseTransformsForPatch(std::shared_ptr<const ThunderEgg::PatchInfo<D>> pinfo)
	{
		// get transform types for each axis
		std::array<DftType, D> transforms_inv;
		for (size_t axis = 0; axis < D; axis++) {
			if (pinfo->isNeumann(Side<D>::LowerSideOnAxis(axis))
			    && pinfo->isNeumann(Side<D>::HigherSideOnAxis(axis))) {
				transforms_inv[axis] = DftType::DCT_III;
			} else if (pinfo->isNeumann(Side<D>::LowerSideOnAxis(axis))) {
				transforms_inv[axis] = DftType::DCT_IV;
			} else if (pinfo->isNeumann(Side<D>::HigherSideOnAxis(axis))) {
				transforms_inv[axis] = DftType::DST_IV;
			} else {
				transforms_inv[axis] = DftType::DST_III;
			}
		}
		return transforms_inv;
	}
	/**
	 * @brief Get an array of eigenvalues for a patch
	 *
	 * @param pinfo the patch
	 * @return std::valarray<double> the eigen values
	 */
	std::valarray<double> getEigenValues(std::shared_ptr<const PatchInfo<D>> pinfo)
	{
		std::valarray<double> retval(this->domain->getNumCellsInPatch());

		std::valarray<size_t> all_strides(D);
		size_t                curr_stride = 1;
		for (size_t i = 0; i < D; i++) {
			all_strides[i] = curr_stride;
			curr_stride *= pinfo->ns[i];
		}

		for (size_t axis = 0; axis < D; axis++) {
			int    n = pinfo->ns[axis];
			double h = pinfo->spacings[axis];

			std::valarray<size_t> sizes(D - 1);
			std::valarray<size_t> strides(D - 1);
			size_t                ones_size = 1;
			for (size_t i = 0; i < axis; i++) {
				strides[i] = all_strides[i];
				sizes[i]   = pinfo->ns[i];
				ones_size *= pinfo->ns[i];
			}
			int input_stride = (int) all_strides[axis];
			for (size_t i = axis + 1; i < D; i++) {
				strides[i - 1] = all_strides[i];
				sizes[i - 1]   = pinfo->ns[i];
				ones_size *= pinfo->ns[i];
			}

			std::valarray<double> ones(ones_size);
			ones = 1;

			if (pinfo->isNeumann(Side<D>::LowerSideOnAxis(axis))
			    && pinfo->isNeumann(Side<D>::HigherSideOnAxis(axis))) {
				for (int xi = 0; xi < n; xi++) {
					retval[std::gslice(xi * input_stride, sizes, strides)]
					-= 4 / (h * h) * pow(sin(xi * M_PI / (2 * n)), 2) * ones;
				}
			} else if (pinfo->isNeumann(Side<D>::LowerSideOnAxis(axis))
			           || pinfo->isNeumann(Side<D>::HigherSideOnAxis(axis))) {
				for (int xi = 0; xi < n; xi++) {
					retval[std::gslice(xi * input_stride, sizes, strides)]
					-= 4 / (h * h) * pow(sin((xi + 0.5) * M_PI / (2 * n)), 2) * ones;
				}
			} else {
				for (int xi = 0; xi < n; xi++) {
					retval[std::gslice(xi * input_stride, sizes, strides)]
					-= 4 / (h * h) * pow(sin((xi + 1) * M_PI / (2 * n)), 2) * ones;
				}
			}
		}
		return retval;
	}
	/**
	 * @brief add a patch to the solver
	 *
	 * This will calculate the necessary coefficients needed for the patch
	 *
	 * @param pinfo the patch
	 */
	void addPatch(std::shared_ptr<const PatchInfo<D>> pinfo)
	{
		if (plan1.count(pinfo) == 0) {
			plan1[pinfo]      = plan(getTransformsForPatch(pinfo));
			plan2[pinfo]      = plan(getInverseTransformsForPatch(pinfo));
			eigen_vals[pinfo] = getEigenValues(pinfo);
		}
	}

	public:
	/**
	 * @brief Construct a new DFTPatchSolver object
	 *
	 * @param op_in the Poisson PatchOperator that cooresponds to this DFTPatchSolver
	 */
	explicit DFTPatchSolver(std::shared_ptr<const PatchOperator<D>> op_in)
	: PatchSolver<D>(op_in->getDomain(), op_in->getGhostFiller()), op(op_in)
	{
		f_copy = std::make_shared<ValVector<D>>(MPI_COMM_SELF, this->domain->getNs(), 0, 1, 1);
		tmp    = std::make_shared<ValVector<D>>(MPI_COMM_SELF, this->domain->getNs(), 0, 1, 1);
		if (!(D % 2)) {
			local_tmp
			= std::make_shared<ValVector<D>>(MPI_COMM_SELF, this->domain->getNs(), 0, 1, 1);
		}
		// process patches
		for (auto pinfo : this->domain->getPatchInfoVector()) {
			addPatch(pinfo);
		}
	}
	void solveSinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo, LocalData<D> u,
	                      const LocalData<D> f) const override
	{
		LocalData<D> f_copy_ld = f_copy->getLocalData(0, 0);
		LocalData<D> tmp_ld    = tmp->getLocalData(0, 0);

		nested_loop<D>(f_copy_ld.getStart(), f_copy_ld.getEnd(),
		               [&](std::array<int, D> coord) { f_copy_ld[coord] = f[coord]; });

		const std::vector<LocalData<D>> us = {u};
		std::vector<LocalData<D>>       fs = {f_copy_ld};
		op->addGhostToRHS(pinfo, us, fs);

		executePlan(plan1.at(pinfo), f_copy_ld, tmp_ld);

		tmp->getValArray() /= eigen_vals.at(pinfo);

		if (pinfo->neumann.all()) {
			tmp->getValArray()[0] = 0;
		}

		executePlan(plan2.at(pinfo), tmp_ld, u);

		double scale = 1;
		for (size_t axis = 0; axis < D; axis++) {
			scale *= 2.0 / this->domain->getNs()[axis];
		}
		nested_loop<D>(u.getStart(), u.getEnd(),
		               [&](std::array<int, D> coord) { u[coord] *= scale; });
	}
};

extern template class DFTPatchSolver<2>;
extern template class DFTPatchSolver<3>;
} // namespace Poisson
} // namespace ThunderEgg
#endif
