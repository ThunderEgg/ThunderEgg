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

#ifndef THUNDEREGG_GMG_CYCLEBUILDER_H
#define THUNDEREGG_GMG_CYCLEBUILDER_H
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/GMG/VCycle.h>
#include <ThunderEgg/GMG/WCycle.h>
#include <ThunderEgg/RuntimeError.h>
namespace ThunderEgg
{
namespace GMG
{
/**
 * @brief Builder for building GMG cycles.
 *
 * User will provide Operator,Smoother,Restrictor, and Interpolator objects for each level.
 *
 * addFinestLevel has to be called first, then addIntermediateLevel (if there are any), and finally
 * addCoarsestLevel. Then getCycle can be called to get the completed Cycle. If these are called in
 * the wrong order, an exception will be thrown.
 *
 * @tparam D the number of cartesian dimensions
 */
template <size_t D> class CycleBuilder
{
	public:
	/**
	 * @brief Construct a new CycleBuilder object
	 *
	 * @param opts the options to use for the GMG cycle
	 */
	CycleBuilder(const CycleOpts &opts) {}
	/**
	 * @brief Add the finest level to the Cycle
	 *
	 * @param op the Operator for the level
	 * @param smoother the Smoother for the level
	 * @param restrictor the Restrictor that restricts from this level to the coarser level
	 * @param vg the VectorGenerator for the level
	 */
	void addFinestLevel(std::shared_ptr<Operator<D>> op, std::shared_ptr<Smoother<D>> smoother,
	                    std::shared_ptr<Restrictor<D>>      restrictor,
	                    std::shared_ptr<VectorGenerator<D>> vg)
	{
		if (op == nullptr) {
			throw RuntimeError("Operator is nullptr");
		}
		if (smoother == nullptr) {
			throw RuntimeError("Smoother is nullptr");
		}
		if (restrictor == nullptr) {
			throw RuntimeError("Restrictor is nullptr");
		}
		if (vg == nullptr) {
			throw RuntimeError("VectorGenerator is nullptr");
		}
	}
	/**
	 * @brief Add the next intermediate level to the Cycle
	 *
	 * @param op the Operator for the level
	 * @param smoother the Smoother for the level
	 * @param restrictor the Restrictor that restricts from this level to the coarser level
	 * @param interpolator the Interpolator that restricts from this level to the finer level
	 * @param vg the VectorGenerator for the level
	 */
	void addIntermediateLevel(std::shared_ptr<Operator<D>>        op,
	                          std::shared_ptr<Smoother<D>>        smoother,
	                          std::shared_ptr<Restrictor<D>>      restrictor,
	                          std::shared_ptr<Interpolator<D>>    interpolator,
	                          std::shared_ptr<VectorGenerator<D>> vg)
	{
		if (op == nullptr) {
			throw RuntimeError("Operator is nullptr");
		}
		if (smoother == nullptr) {
			throw RuntimeError("Smoother is nullptr");
		}
		if (restrictor == nullptr) {
			throw RuntimeError("Restrictor is nullptr");
		}
		if (interpolator == nullptr) {
			throw RuntimeError("Interpolator is nullptr");
		}
		if (vg == nullptr) {
			throw RuntimeError("VectorGenerator is nullptr");
		}
	}
	/**
	 * @brief Add the next intermediate level to the Cycle
	 *
	 * @param op the Operator for the level
	 * @param smoother the Smoother for the level
	 * @param interpolator the Interpolator that restricts from this level to the finer level
	 * @param vg the VectorGenerator for the level
	 */
	void addCoarsestLevel(std::shared_ptr<Operator<D>> op, std::shared_ptr<Smoother<D>> smoother,
	                      std::shared_ptr<Interpolator<D>>    interpolator,
	                      std::shared_ptr<VectorGenerator<D>> vg)
	{
		if (op == nullptr) {
			throw RuntimeError("Operator is nullptr");
		}
		if (smoother == nullptr) {
			throw RuntimeError("Smoother is nullptr");
		}
		if (interpolator == nullptr) {
			throw RuntimeError("Interpolator is nullptr");
		}
		if (vg == nullptr) {
			throw RuntimeError("VectorGenerator is nullptr");
		}
	}
	/**
	 * @brief Get the completed Cycle object
	 *
	 * @return std::shared_ptr<Cycle<D>> the completed Cycle, will throw an exception if it is
	 * incomplete
	 */
	std::shared_ptr<Cycle<D>> getCycle() const {}
};
} // namespace GMG
} // namespace ThunderEgg

#endif