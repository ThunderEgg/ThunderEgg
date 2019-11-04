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

#include <Thunderegg/GMG/InterLevelComm.h>
#include <Thunderegg/GMG/Level.h>
#include <Thunderegg/GMG/VCycle.h>
#include <Thunderegg/GMG/WCycle.h>
namespace Thunderegg
{
namespace GMG
{
/**
 * @brief Factory for producing GMG cycles.
 *
 * User has to provide functions that return Restrictors, Interpolators, Smoothers, and Operators
 * for the GMG Levels.
 *
 * @tparam D the number of cartesian dimensions
 */
template <size_t D> class CycleFactory
{
	public:
	/**
	 * @brief Function that returns a Restrictor for a given level
	 */
	using RestrictorGenerator
	= std::function<std::shared_ptr<const Restrictor<D>>(std::shared_ptr<const Level<D>> level)>;
	/**
	 * @brief Function that returns an Interpolator for a given level
	 */
	using InterpolatorGenerator
	= std::function<std::shared_ptr<const Interpolator<D>>(std::shared_ptr<const Level<D>> level)>;
	/**
	 * @brief Function that returns a Smoother for a given level
	 */
	using SmootherGenerator
	= std::function<std::shared_ptr<const Smoother<D>>(std::shared_ptr<const Level<D>> level)>;
	/**
	 * @brief Function that returns an Operator for a given level
	 */
	using OperatorGenerator
	= std::function<std::shared_ptr<const Operator<D>>(std::shared_ptr<const Level<D>> level)>;
	static std::shared_ptr<Cycle<2>>
	/**
	 * @brief Generate a Cycle object
	 * 
	 * @param opts Cycle options
	 * @param domain_gen DomainGenerator for a GMG level
	 * @param restrictor_gen returns a Restrictor for a given level (restricts to coarser level)
	 * @param interpolator_gen returns an Interpolator for a given level (interpolates to finer level)
	 * @param smoother_gen returns a Smoother for a given level
	 * @param operator_gen returns an operator for a given level
	 */
	getCycle(const CycleOpts &opts, std::shared_ptr<DomainGenerator<D>> domain_gen,
	         RestrictorGenerator restrictor_gen, InterpolatorGenerator interpolator_gen,
	         SmootherGenerator smoother_gen, OperatorGenerator operator_gen)
	{
		int size;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		// finest level
		std::shared_ptr<Level<D>> finest_level;
		{
			std::shared_ptr<Domain<D>>             domain = domain_gen->getFinestDomain();
			std::shared_ptr<Schur::SchurHelper<D>> sh(new Schur::SchurHelper<2>(domain));
			std::shared_ptr<VectorGenerator<D>>    vg(new DomainVG<2>(domain));
			finest_level.reset(new Level<D>(domain, vg));
			finest_level->setOperator(operator_gen(finest_level));
			finest_level->setSmoother(smoother_gen(finest_level));
		}
		std::shared_ptr<Level<D>> finer_level = finest_level;
		// other levels
		int curr_level = 1;
		while (domain_gen->hasCoarserDomain()
		       && (opts.max_levels <= 0 || curr_level < opts.max_levels)) {
			// create new level
			std::shared_ptr<Domain<D>> domain = domain_gen->getCoarserDomain();
			if ((domain->getNumGlobalPatches() + 0.0) / size < opts.patches_per_proc) { break; }
			std::shared_ptr<VectorGenerator<D>> vg(new DomainVG<2>(domain));
			std::shared_ptr<Level<D>>           coarser_level(new Level<2>(domain, vg));

			// link levels
			coarser_level->setFiner(finer_level);
			finer_level->setCoarser(coarser_level);

			// set restrictor and interpolator operators
			finer_level->setRestrictor(restrictor_gen(finer_level));
			coarser_level->setInterpolator(interpolator_gen(coarser_level));

			coarser_level->setOperator(operator_gen(coarser_level));
			coarser_level->setSmoother(smoother_gen(coarser_level));

			curr_level++;
			finer_level = coarser_level;
		}
		std::shared_ptr<Cycle<D>> cycle;
		if (opts.cycle_type == "V") {
			cycle.reset(new VCycle<D>(finest_level, opts));
		} else if (opts.cycle_type == "W") {
			cycle.reset(new WCycle<D>(finest_level, opts));
		} else {
			throw 3;
		}
		return cycle;
	}
};
} // namespace GMG
} // namespace Thunderegg
