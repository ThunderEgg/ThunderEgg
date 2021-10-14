/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2018-2021 Scott Aiton
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

#ifndef THUNDEREGG_GMG_VCYCLE_H
#define THUNDEREGG_GMG_VCYCLE_H
/**
 * @file
 *
 * @brief VCycle class
 */
#include <ThunderEgg/GMG/Cycle.h>
#include <ThunderEgg/GMG/CycleOpts.h>
namespace ThunderEgg::GMG
{
/**
 * @brief Implementation of a V-cycle
 */
template <int D> class VCycle : public Cycle<D>
{
	private:
	int num_pre_sweeps    = 1;
	int num_post_sweeps   = 1;
	int num_coarse_sweeps = 1;

	protected:
	/*
	 * @brief Implements V-cycle. Pre-smooth, visit coarser level and then post-smooth.
	 */
	void visit(const Level<D> &level, const Vector<D> &f, Vector<D> &u) const override
	{
		if (level.coarsest()) {
			for (int i = 0; i < num_coarse_sweeps; i++) {
				level.getSmoother().smooth(f, u);
			}
		} else {
			for (int i = 0; i < num_pre_sweeps; i++) {
				level.getSmoother().smooth(f, u);
			}

			Vector<D> coarser_f = this->restrict(level, f, u);

			const Level<D> &coarser_level = level.getCoarser();
			Vector<D>       coarser_u     = coarser_f.getZeroClone();

			this->visit(coarser_level, coarser_f, coarser_u);

			coarser_level.getInterpolator().interpolate(coarser_u, u);

			for (int i = 0; i < num_post_sweeps; i++) {
				level.getSmoother().smooth(f, u);
			}
		}
	}

	public:
	/**
	 * @brief Create new V-cycle
	 *
	 * @param finest_level a pointer to the finest level
	 */
	VCycle(const Level<D> &finest_level, const CycleOpts &opts) : Cycle<D>(finest_level)
	{
		num_pre_sweeps    = opts.pre_sweeps;
		num_post_sweeps   = opts.post_sweeps;
		num_coarse_sweeps = opts.coarse_sweeps;
	}
	/**
	 * @brief Get a clone of this VCycle
	 *
	 * @return VCycle<D>* a newly allocated copy
	 */
	VCycle<D> *clone() const override
	{
		return new VCycle(*this);
	}
};
extern template class VCycle<2>;
extern template class VCycle<3>;
} // namespace ThunderEgg::GMG
#endif