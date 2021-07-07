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
#include <ThunderEgg/PatchArray.h>
#include <ThunderEgg/PatchOperator.h>
#include <ThunderEgg/PatchSolver.h>
#include <ThunderEgg/Vector.h>
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
	using CompareFunction = std::function<bool(const PatchInfo<D> &, const PatchInfo<D> &)>;

	/**
	 * @brief The patch opertar that we are solving for
	 */
	std::shared_ptr<const PatchOperator<D>> op;
	/**
	 * @brief Map of patchinfo to DFT plan
	 */
	std::map<const PatchInfo<D>, std::shared_ptr<fftw_plan>, CompareFunction> plan1;
	/**
	 * @brief Map of patchinfo to inverse DFT plan
	 */
	std::map<const PatchInfo<D>, std::shared_ptr<fftw_plan>, CompareFunction> plan2;
	/**
	 * @brief Map of PatchInfo object to it's respective eigenvalue array.
	 */
	std::map<const PatchInfo<D>, PatchArray<D>, CompareFunction> eigen_vals;
	/**
	 * @brief Neumann boundary conditions for domain
	 */
	std::bitset<Side<D>::number_of> neumann;

	/**
	 * @brief Return if a patch has a neumann boundary condition on a particular side
	 *
	 * @param pinfo the patch
	 * @param s the side
	 * @return true if neumann
	 * @return false if not neumann
	 */
	bool patchIsNeumannOnSide(const PatchInfo<D> &pinfo, Side<D> s)
	{
		return !pinfo.hasNbr(s) && neumann[s.getIndex()];
	}
	/**
	 * @brief Get the fft transform types for a patch
	 *
	 * @param pinfo the patch
	 * @return std::array<fftw_r2r_kind, D> an array of tranforms for each axis, the order of
	 * dimensions is reversed because FFTW uses row-major format
	 */
	std::array<fftw_r2r_kind, D> getTransformsForPatch(const PatchInfo<D> &pinfo)
	{
		// get transform types for each axis
		std::array<fftw_r2r_kind, D> transforms;
		for (size_t axis = 0; axis < D; axis++) {
			if (patchIsNeumannOnSide(pinfo, LowerSideOnAxis<D>(axis)) && patchIsNeumannOnSide(pinfo, HigherSideOnAxis<D>(axis))) {
				transforms[D - 1 - axis] = FFTW_REDFT10;
			} else if (patchIsNeumannOnSide(pinfo, LowerSideOnAxis<D>(axis))) {
				transforms[D - 1 - axis] = FFTW_REDFT11;
			} else if (patchIsNeumannOnSide(pinfo, HigherSideOnAxis<D>(axis))) {
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
	std::array<fftw_r2r_kind, D> getInverseTransformsForPatch(const PatchInfo<D> &pinfo)
	{
		// get transform types for each axis
		std::array<fftw_r2r_kind, D> transforms_inv;
		for (size_t axis = 0; axis < D; axis++) {
			if (patchIsNeumannOnSide(pinfo, LowerSideOnAxis<D>(axis)) && patchIsNeumannOnSide(pinfo, HigherSideOnAxis<D>(axis))) {
				transforms_inv[D - 1 - axis] = FFTW_REDFT01;
			} else if (patchIsNeumannOnSide(pinfo, LowerSideOnAxis<D>(axis))) {
				transforms_inv[D - 1 - axis] = FFTW_REDFT11;
			} else if (patchIsNeumannOnSide(pinfo, HigherSideOnAxis<D>(axis))) {
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
	 * @return PatchArray<D> the eigen values
	 */
	PatchArray<D> getEigenValues(const PatchInfo<D> &pinfo)
	{
		PatchArray<D> retval(this->getDomain().getNs(), 1, 0);

		std::valarray<size_t> all_strides(D);
		size_t                curr_stride = 1;
		for (size_t i = 0; i < D; i++) {
			all_strides[i] = curr_stride;
			curr_stride *= pinfo.ns[i];
		}

		for (size_t axis = 0; axis < D; axis++) {
			int    n = pinfo.ns[axis];
			double h = pinfo.spacings[axis];

			if (patchIsNeumannOnSide(pinfo, LowerSideOnAxis<D>(axis)) && patchIsNeumannOnSide(pinfo, HigherSideOnAxis<D>(axis))) {
				for (int xi = 0; xi < n; xi++) {
					double          val   = 4 / (h * h) * pow(sin(xi * M_PI / (2 * n)), 2);
					View<double, D> slice = retval.getSliceOn(Side<D>(2 * axis), {xi});
					loop_over_interior_indexes<D>(slice, [&](const std::array<int, D> &coord) { slice[coord] -= val; });
				}
			} else if (patchIsNeumannOnSide(pinfo, LowerSideOnAxis<D>(axis)) || patchIsNeumannOnSide(pinfo, HigherSideOnAxis<D>(axis))) {
				for (int xi = 0; xi < n; xi++) {
					double          val   = 4 / (h * h) * pow(sin((xi + 0.5) * M_PI / (2 * n)), 2);
					View<double, D> slice = retval.getSliceOn(Side<D>(2 * axis), {xi});
					loop_over_interior_indexes<D>(slice, [&](const std::array<int, D> &coord) { slice[coord] -= val; });
				}
			} else {
				for (int xi = 0; xi < n; xi++) {
					double          val   = 4 / (h * h) * pow(sin((xi + 1) * M_PI / (2 * n)), 2);
					View<double, D> slice = retval.getSliceOn(Side<D>(2 * axis), {xi});
					loop_over_interior_indexes<D>(slice, [&](const std::array<int, D> &coord) { slice[coord] -= val; });
				}
			}
		}
		return retval;
	}

	public:
	/**
	 * @brief Construct a new FftwPatchSolver object
	 *
	 * @param op the Poisson PatchOperator that cooresponds to this DftPatchSolver
	 * @param neumann true if domain has neumann boundary conditions on a side
	 */
	FFTWPatchSolver(const PatchOperator<D> &op, std::bitset<Side<D>::number_of> neumann)
	: PatchSolver<D>(op.getDomain(), op.getGhostFiller()),
	  op(op.clone()),
	  neumann(neumann)
	{
		CompareFunction compare = [&](const PatchInfo<D> &a, const PatchInfo<D> &b) {
			std::bitset<Side<D>::number_of> a_neumann;
			std::bitset<Side<D>::number_of> b_neumann;
			for (Side<D> s : Side<D>::getValues()) {
				a_neumann[s.getIndex()] = patchIsNeumannOnSide(a, s);
				b_neumann[s.getIndex()] = patchIsNeumannOnSide(b, s);
			}
			return std::forward_as_tuple(a_neumann.to_ulong(), a.spacings[0]) < std::forward_as_tuple(b_neumann.to_ulong(), b.spacings[0]);
		};

		plan1      = std::map<const PatchInfo<D>, std::shared_ptr<fftw_plan>, CompareFunction>(compare);
		plan2      = std::map<const PatchInfo<D>, std::shared_ptr<fftw_plan>, CompareFunction>(compare);
		eigen_vals = std::map<const PatchInfo<D>, PatchArray<D>, CompareFunction>(compare);

		// process patches
		for (auto pinfo : this->getDomain().getPatchInfoVector()) {
			addPatch(pinfo);
		}
	}
	/**
	 * @brief Clone this patch solver
	 *
	 * @return FFTWPatchSolver<D>* a newly allocated copy of this patch solver
	 */
	FFTWPatchSolver<D> *clone() const override
	{
		return new FFTWPatchSolver<D>(*this);
	}
	void solveSinglePatch(const PatchInfo<D> &pinfo, const PatchView<const double, D> &f_view, const PatchView<double, D> &u_view) const override
	{
		PatchArray<D> f_copy(pinfo.ns, 1, 0);
		PatchArray<D> tmp(pinfo.ns, 1, 0);
		PatchArray<D> sol(pinfo.ns, 1, 0);

		loop_over_interior_indexes<D + 1>(f_copy, [&](std::array<int, D + 1> coord) { f_copy[coord] = f_view[coord]; });

		op->modifyRHSForZeroDirichletAtInternalBoundaries(pinfo, u_view, f_copy.getView());

		fftw_execute_r2r(*plan1.at(pinfo), &f_copy[f_copy.getStart()], &tmp[tmp.getStart()]);

		const PatchArray<D> &eigen_vals_view = eigen_vals.at(pinfo);
		loop_over_interior_indexes<D + 1>(tmp, [&](std::array<int, D + 1> coord) { tmp[coord] /= eigen_vals_view[coord]; });

		if (neumann.all() && !pinfo.hasNbr()) {
			tmp[tmp.getStart()] = 0;
		}

		fftw_execute_r2r(*plan2.at(pinfo), &tmp[tmp.getStart()], &sol[sol.getStart()]);

		double scale = 1;
		for (size_t axis = 0; axis < D; axis++) {
			scale *= 2.0 * this->getDomain().getNs()[axis];
		}
		loop_over_interior_indexes<D + 1>(u_view, [&](std::array<int, D + 1> coord) { u_view[coord] = sol[coord] / scale; });
	}
	/**
	 * @brief add a patch to the solver
	 *
	 * This will calculate the necessary coefficients needed for the patch
	 *
	 * @param pinfo the patch
	 */
	void addPatch(const PatchInfo<D> &pinfo)
	{
		if (plan1.count(pinfo) == 0) {
			// revers ns because FFTW is row major
			std::array<int, D> ns_reversed;
			for (size_t i = 0; i < D; i++) {
				ns_reversed[D - 1 - i] = pinfo.ns[i];
			}
			std::array<fftw_r2r_kind, D> transforms     = getTransformsForPatch(pinfo);
			std::array<fftw_r2r_kind, D> transforms_inv = getInverseTransformsForPatch(pinfo);

			PatchArray<D> f_copy(pinfo.ns, 1, 0);
			PatchArray<D> tmp(pinfo.ns, 1, 0);
			PatchArray<D> sol(pinfo.ns, 1, 0);

			fftw_plan *fftw_plan1 = new fftw_plan();

			*fftw_plan1 = fftw_plan_r2r(D,
			                            ns_reversed.data(),
			                            &f_copy[f_copy.getStart()],
			                            &tmp[tmp.getStart()],
			                            transforms.data(),
			                            FFTW_MEASURE | FFTW_DESTROY_INPUT | FFTW_UNALIGNED);

			plan1[pinfo] = std::shared_ptr<fftw_plan>(fftw_plan1, [](fftw_plan *plan) {
				fftw_destroy_plan(*plan);
				delete plan;
			});

			fftw_plan *fftw_plan2 = new fftw_plan();

			*fftw_plan2 = fftw_plan_r2r(D,
			                            ns_reversed.data(),
			                            &tmp[tmp.getStart()],
			                            &sol[sol.getStart()],
			                            transforms_inv.data(),
			                            FFTW_MEASURE | FFTW_DESTROY_INPUT | FFTW_UNALIGNED);

			plan2[pinfo] = std::shared_ptr<fftw_plan>(fftw_plan2, [](fftw_plan *plan) {
				fftw_destroy_plan(*plan);
				delete plan;
			});

			eigen_vals.emplace(pinfo, getEigenValues(pinfo));
		}
	}
	/**
	 * @brief Get the neumann boundary conditions for this operator
	 *
	 * @return std::bitset<Side<D>::number_of> the boundary conditions
	 */
	std::bitset<Side<D>::number_of> getNeumann() const
	{
		return neumann;
	}
};
extern template class FFTWPatchSolver<2>;
extern template class FFTWPatchSolver<3>;
} // namespace Poisson
} // namespace ThunderEgg
#endif