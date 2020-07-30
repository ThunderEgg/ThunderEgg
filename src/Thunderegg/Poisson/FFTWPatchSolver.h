/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
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

#ifndef THUNDEREGG_POISSON_SCHUR_FFTWPATCHSOLVER_H
#define THUNDEREGG_POISSON_SCHUR_FFTWPATCHSOLVER_H
#include <Thunderegg/PatchOperator.h>
#include <Thunderegg/PatchSolver.h>
#include <Thunderegg/ValVector.h>
#include <bitset>
#include <fftw3.h>
#include <map>
namespace Thunderegg
{
namespace Poisson
{
/**
 * @brief This patch solver uses FFT transforms to solve for the Poisson equation
 *
 * @tparam D the number of Cartesian dimensions
 */
template <size_t D> class FFTWPatchSolver : public PatchSolver<D>
{
	private:
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
	std::map<std::shared_ptr<const PatchInfo<D>>, fftw_plan, CompareByBoundaryAndSpacings> plan1;
	/**
	 * @brief Map of patchinfo to inverse DFT plan
	 */
	std::map<std::shared_ptr<const PatchInfo<D>>, fftw_plan, CompareByBoundaryAndSpacings> plan2;
	/**
	 * @brief Temporary copy for the modified right hand side
	 */
	std::shared_ptr<ValVector<D>> f_copy;
	/**
	 * @brief Temporary work vector
	 */
	std::shared_ptr<ValVector<D>> tmp;
	/*
	 * @brief Temporary work vector for solution
	 */
	std::shared_ptr<ValVector<D>> sol;
	/**
	 * @brief Map of PatchInfo object to it's respective eigenvalue array.
	 */
	std::map<std::shared_ptr<const PatchInfo<D>>, std::valarray<double>,
	         CompareByBoundaryAndSpacings>
	eigen_vals;
	/**
	 * @brief Get the fft transform types for a patch
	 *
	 * @param pinfo the patch
	 * @return std::array<fftw_r2r_kind, D> an array of tranforms for each axis, the order of
	 * dimensions is reversed because FFTW uses row-major format
	 */
	std::array<fftw_r2r_kind, D>
	getTransformsForPatch(std::shared_ptr<const Thunderegg::PatchInfo<D>> pinfo)
	{
		// get transform types for each axis
		std::array<fftw_r2r_kind, D> transforms;
		for (size_t axis = 0; axis < D; axis++) {
			if (pinfo->isNeumann(Side<D>::LowerSideOnAxis(axis))
			    && pinfo->isNeumann(Side<D>::HigherSideOnAxis(axis))) {
				transforms[D - 1 - axis] = FFTW_REDFT10;
			} else if (pinfo->isNeumann(Side<D>::LowerSideOnAxis(axis))) {
				transforms[D - 1 - axis] = FFTW_REDFT11;
			} else if (pinfo->isNeumann(Side<D>::HigherSideOnAxis(axis))) {
				transforms[D - 1 - axis] = FFTW_RODFT11;
			} else {
				transforms[D - 1 - axis] = FFTW_RODFT10;
			}
		}
		return transforms;
	}
	/**
	 * @brief Get the inverse fft transform types for a patch
	 *
	 * @param pinfo the patch
	 * @return std::array<fftw_r2r_kind, D> an array of tranforms for each axis, the order of
	 * dimensions is reversed because FFTW uses row-major format
	 */
	std::array<fftw_r2r_kind, D>
	getInverseTransformsForPatch(std::shared_ptr<const Thunderegg::PatchInfo<D>> pinfo)
	{
		// get transform types for each axis
		std::array<fftw_r2r_kind, D> transforms_inv;
		for (size_t axis = 0; axis < D; axis++) {
			if (pinfo->isNeumann(Side<D>::LowerSideOnAxis(axis))
			    && pinfo->isNeumann(Side<D>::HigherSideOnAxis(axis))) {
				transforms_inv[D - 1 - axis] = FFTW_REDFT01;
			} else if (pinfo->isNeumann(Side<D>::LowerSideOnAxis(axis))) {
				transforms_inv[D - 1 - axis] = FFTW_REDFT11;
			} else if (pinfo->isNeumann(Side<D>::HigherSideOnAxis(axis))) {
				transforms_inv[D - 1 - axis] = FFTW_RODFT11;
			} else {
				transforms_inv[D - 1 - axis] = FFTW_RODFT01;
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
			// revers ns because FFTW is row major
			std::array<int, D> ns_reversed;
			for (size_t i = 0; i < D; i++) {
				ns_reversed[D - 1 - i] = pinfo->ns[i];
			}
			std::array<fftw_r2r_kind, D> transforms     = getTransformsForPatch(pinfo);
			std::array<fftw_r2r_kind, D> transforms_inv = getInverseTransformsForPatch(pinfo);

			plan1[pinfo] = fftw_plan_r2r(D, ns_reversed.data(), &f_copy->vec[0], &tmp->vec[0],
			                             transforms.data(), FFTW_MEASURE | FFTW_DESTROY_INPUT);
			plan2[pinfo] = fftw_plan_r2r(D, ns_reversed.data(), &tmp->vec[0], &sol->vec[0],
			                             transforms_inv.data(), FFTW_MEASURE | FFTW_DESTROY_INPUT);

			eigen_vals[pinfo] = getEigenValues(pinfo);
		}
	}

	public:
	/**
	 * @brief Construct a new FftwPatchSolver object
	 *
	 * @param op_in the Poisson PatchOperator that cooresponds to this DftPatchSolver
	 */
	explicit FFTWPatchSolver(std::shared_ptr<const PatchOperator<D>> op_in)
	: PatchSolver<D>(op_in->getDomain(), op_in->getGhostFiller()), op(op_in)
	{
		f_copy = std::make_shared<ValVector<D>>(MPI_COMM_SELF, this->domain->getNs(), 0, 1);
		tmp    = std::make_shared<ValVector<D>>(MPI_COMM_SELF, this->domain->getNs(), 0, 1);
		sol    = std::make_shared<ValVector<D>>(MPI_COMM_SELF, this->domain->getNs(), 0, 1);
		// process patches
		for (auto pinfo : this->domain->getPatchInfoVector()) {
			addPatch(pinfo);
		}
	}
	void solveSinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo, LocalData<D> u,
	                      const LocalData<D> f) const override
	{
		LocalData<D> f_copy_ld = f_copy->getLocalData(0);

		nested_loop<D>(f_copy_ld.getStart(), f_copy_ld.getEnd(),
		               [&](std::array<int, D> coord) { f_copy_ld[coord] = f[coord]; });

		op->addGhostToRHS(pinfo, u, f_copy_ld);

		fftw_execute(plan1.at(pinfo));

		tmp->vec /= eigen_vals.at(pinfo);

		if (pinfo->neumann.all()) {
			tmp->vec[0] = 0;
		}

		fftw_execute(plan2.at(pinfo));

		LocalData<D> sol_ld = sol->getLocalData(0);

		double scale = 1;
		for (size_t axis = 0; axis < D; axis++) {
			scale *= 2.0 * this->domain->getNs()[axis];
		}
		nested_loop<D>(u.getStart(), u.getEnd(),
		               [&](std::array<int, D> coord) { u[coord] = sol_ld[coord] / scale; });
	}
};
extern template class FFTWPatchSolver<2>;
extern template class FFTWPatchSolver<3>;
} // namespace Poisson
} // namespace Thunderegg
#endif