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

#ifndef THUNDEREGG_POISSON_SCHUR_FFTWPATCHSOLVER_H
#define THUNDEREGG_POISSON_SCHUR_FFTWPATCHSOLVER_H
#include <ThunderEgg/PatchOperator.h>
#include <ThunderEgg/PatchSolver.h>
#include <ThunderEgg/ValVector.h>
#include <bitset>
#include <fftw3.h>
#include <map>
namespace ThunderEgg
{
namespace Poisson
{
/**
 * @brief This patch solver uses FFT transforms to solve for the Poisson equation
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> class FFTWPatchSolver : public PatchSolver<D>
{
	private:
	/**
	 * @brief Comparator used in the maps, patches with the same spacings and boundary conditions
	 * will be equal
	 */
	std::function<bool(const std::shared_ptr<const PatchInfo<D>> &, const std::shared_ptr<const PatchInfo<D>> &)> CompareByBoundaryAndSpacings
	= [&](const std::shared_ptr<const PatchInfo<D>> &a, const std::shared_ptr<const PatchInfo<D>> &b) -> bool {
		std::bitset<Side<D>::num_sides> a_neumann;
		std::bitset<Side<D>::num_sides> b_neumann;
		for (Side<D> s : Side<D>::getValues()) {
			a_neumann[s.getIndex()] = patchIsNeumannOnSide(a, s);
			b_neumann[s.getIndex()] = patchIsNeumannOnSide(b, s);
		}
		return std::forward_as_tuple(a_neumann.to_ulong(), a->spacings[0]) < std::forward_as_tuple(b_neumann.to_ulong(), b->spacings[0]);
	};
	/**
	 * @brief The patch opertar that we are solving for
	 */
	std::shared_ptr<const PatchOperator<D>> op;
	/**
	 * @brief Map of patchinfo to DFT plan
	 */
	std::map<std::shared_ptr<const PatchInfo<D>>, fftw_plan, decltype(CompareByBoundaryAndSpacings)> plan1;
	/**
	 * @brief Map of patchinfo to inverse DFT plan
	 */
	std::map<std::shared_ptr<const PatchInfo<D>>, fftw_plan, decltype(CompareByBoundaryAndSpacings)> plan2;
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
	std::map<std::shared_ptr<const PatchInfo<D>>, std::valarray<double>, decltype(CompareByBoundaryAndSpacings)> eigen_vals;
	/**
	 * @brief Neumann boundary conditions for domain
	 */
	std::bitset<Side<D>::num_sides> neumann;

	/**
	 * @brief Return if a patch has a neumann boundary condition on a particular side
	 *
	 * @param pinfo the patch
	 * @param s the side
	 * @return true if neumann
	 * @return false if not neumann
	 */
	bool patchIsNeumannOnSide(const std::shared_ptr<const PatchInfo<D>> &pinfo, Side<D> s)
	{
		return !pinfo->hasNbr(s) && neumann[s.getIndex()];
	}
	/**
	 * @brief Get the fft transform types for a patch
	 *
	 * @param pinfo the patch
	 * @return std::array<fftw_r2r_kind, D> an array of tranforms for each axis, the order of
	 * dimensions is reversed because FFTW uses row-major format
	 */
	std::array<fftw_r2r_kind, D> getTransformsForPatch(std::shared_ptr<const ThunderEgg::PatchInfo<D>> pinfo)
	{
		// get transform types for each axis
		std::array<fftw_r2r_kind, D> transforms;
		for (size_t axis = 0; axis < D; axis++) {
			if (patchIsNeumannOnSide(pinfo, Side<D>::LowerSideOnAxis(axis)) && patchIsNeumannOnSide(pinfo, Side<D>::HigherSideOnAxis(axis))) {
				transforms[D - 1 - axis] = FFTW_REDFT10;
			} else if (patchIsNeumannOnSide(pinfo, Side<D>::LowerSideOnAxis(axis))) {
				transforms[D - 1 - axis] = FFTW_REDFT11;
			} else if (patchIsNeumannOnSide(pinfo, Side<D>::HigherSideOnAxis(axis))) {
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
	std::array<fftw_r2r_kind, D> getInverseTransformsForPatch(std::shared_ptr<const ThunderEgg::PatchInfo<D>> pinfo)
	{
		// get transform types for each axis
		std::array<fftw_r2r_kind, D> transforms_inv;
		for (size_t axis = 0; axis < D; axis++) {
			if (patchIsNeumannOnSide(pinfo, Side<D>::LowerSideOnAxis(axis)) && patchIsNeumannOnSide(pinfo, Side<D>::HigherSideOnAxis(axis))) {
				transforms_inv[D - 1 - axis] = FFTW_REDFT01;
			} else if (patchIsNeumannOnSide(pinfo, Side<D>::LowerSideOnAxis(axis))) {
				transforms_inv[D - 1 - axis] = FFTW_REDFT11;
			} else if (patchIsNeumannOnSide(pinfo, Side<D>::HigherSideOnAxis(axis))) {
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

			if (patchIsNeumannOnSide(pinfo, Side<D>::LowerSideOnAxis(axis)) && patchIsNeumannOnSide(pinfo, Side<D>::HigherSideOnAxis(axis))) {
				for (int xi = 0; xi < n; xi++) {
					retval[std::gslice(xi * input_stride, sizes, strides)] -= 4 / (h * h) * pow(sin(xi * M_PI / (2 * n)), 2) * ones;
				}
			} else if (patchIsNeumannOnSide(pinfo, Side<D>::LowerSideOnAxis(axis)) || patchIsNeumannOnSide(pinfo, Side<D>::HigherSideOnAxis(axis))) {
				for (int xi = 0; xi < n; xi++) {
					retval[std::gslice(xi * input_stride, sizes, strides)] -= 4 / (h * h) * pow(sin((xi + 0.5) * M_PI / (2 * n)), 2) * ones;
				}
			} else {
				for (int xi = 0; xi < n; xi++) {
					retval[std::gslice(xi * input_stride, sizes, strides)] -= 4 / (h * h) * pow(sin((xi + 1) * M_PI / (2 * n)), 2) * ones;
				}
			}
		}
		return retval;
	}

	public:
	/**
	 * @brief Construct a new FftwPatchSolver object
	 *
	 * @param op_in the Poisson PatchOperator that cooresponds to this DftPatchSolver
	 * @param neumann true if domain has neumann boundary conditions on a side
	 */
	FFTWPatchSolver(std::shared_ptr<const PatchOperator<D>> op_in, std::bitset<Side<D>::num_sides> neumann)
	: PatchSolver<D>(op_in->getDomain(), op_in->getGhostFiller()),
	  op(op_in),
	  plan1(CompareByBoundaryAndSpacings),
	  plan2(CompareByBoundaryAndSpacings),
	  eigen_vals(CompareByBoundaryAndSpacings),
	  neumann(neumann)
	{
		f_copy = std::make_shared<ValVector<D>>(MPI_COMM_SELF, this->domain->getNs(), 0, 1, 1);
		tmp    = std::make_shared<ValVector<D>>(MPI_COMM_SELF, this->domain->getNs(), 0, 1, 1);
		sol    = std::make_shared<ValVector<D>>(MPI_COMM_SELF, this->domain->getNs(), 0, 1, 1);
		// process patches
		for (auto pinfo : this->domain->getPatchInfoVector()) {
			addPatch(pinfo);
		}
	}
	void
	solveSinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo, const std::vector<LocalData<D>> &fs, std::vector<LocalData<D>> &us) const override
	{
		LocalData<D> f_copy_ld = f_copy->getLocalData(0, 0);

		nested_loop<D>(f_copy_ld.getStart(), f_copy_ld.getEnd(), [&](std::array<int, D> coord) { f_copy_ld[coord] = fs[0][coord]; });

		std::vector<LocalData<D>> f_copy_lds = {f_copy_ld};
		op->addGhostToRHS(pinfo, us, f_copy_lds);

		fftw_execute(plan1.at(pinfo));

		tmp->getValArray() /= eigen_vals.at(pinfo);

		if (neumann.all() && pinfo->nbr_info == std::array<std::unique_ptr<NbrInfo<D>>, Side<D>::num_sides>()) {
			tmp->getValArray()[0] = 0;
		}

		fftw_execute(plan2.at(pinfo));

		LocalData<D> sol_ld = sol->getLocalData(0, 0);

		double scale = 1;
		for (size_t axis = 0; axis < D; axis++) {
			scale *= 2.0 * this->domain->getNs()[axis];
		}
		nested_loop<D>(us[0].getStart(), us[0].getEnd(), [&](std::array<int, D> coord) { us[0][coord] = sol_ld[coord] / scale; });
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

			plan1[pinfo] = fftw_plan_r2r(
			D, ns_reversed.data(), &f_copy->getValArray()[0], &tmp->getValArray()[0], transforms.data(), FFTW_MEASURE | FFTW_DESTROY_INPUT);
			plan2[pinfo] = fftw_plan_r2r(
			D, ns_reversed.data(), &tmp->getValArray()[0], &sol->getValArray()[0], transforms_inv.data(), FFTW_MEASURE | FFTW_DESTROY_INPUT);

			eigen_vals[pinfo] = getEigenValues(pinfo);
		}
	}
	/**
	 * @brief Get the neumann boundary conditions for this operator
	 *
	 * @return std::bitset<Side<D>::num_sides> the boundary conditions
	 */
	std::bitset<Side<D>::num_sides> getNeumann() const
	{
		return neumann;
	}
};
extern template class FFTWPatchSolver<2>;
extern template class FFTWPatchSolver<3>;
} // namespace Poisson
} // namespace ThunderEgg
#endif