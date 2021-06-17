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

#ifndef THUNDEREGG_GMG_CYCLE_H
#define THUNDEREGG_GMG_CYCLE_H

#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/Vector.h>
#include <list>

namespace ThunderEgg
{
namespace GMG
{
/**
 * @brief Base class for cycles. Includes functions for preparing vectors for finer and coarser
 * levels, and a function to run an iteration of smoothing on a level. Derived cycle classes
 * need to implement the visit function.
 */
template <int D> class Cycle : public Operator<D>
{
	private:
	/**
	 * @brief pointer to the finest level
	 */
	std::shared_ptr<Level<D>> finest_level;

	protected:
	/**
	 * @brief Prepare vectors for coarser level.
	 *
	 * @param level the current level
	 */
	void restrict(const Level<D> &level, const Vector<D> &f, const Vector<D> &u, Vector<D> &coarser_f) const
	{
		// calculate residual
		std::shared_ptr<Vector<D>> r = level.getVectorGenerator()->getNewVector();
		level.getOperator()->apply(u, *r);
		r->scaleThenAdd(-1, f);
		// create vectors for coarser levels
		level.getRestrictor()->restrict(*r, coarser_f);
	}

	/**
	 * @brief Virtual visit function that needs to be implemented in derived classes.
	 *
	 * @param level the level currently begin visited.
	 */
	virtual void visit(const Level<D> &level, const Vector<D> &f, Vector<D> &u) const = 0;

	public:
	/**
	 * @brief Create new cycle object.
	 *
	 * @param finest_level pointer to the finest level object.
	 */
	Cycle(std::shared_ptr<Level<D>> finest_level)
	{
		this->finest_level = finest_level;
	}
	/**
	 * @brief Run one iteration of the cycle.
	 *
	 * @param f the RHS vector.
	 * @param u the current solution vector. Output will be updated solution vector.
	 */
	void apply(const Vector<D> &f, Vector<D> &u) const
	{
		u.setWithGhost(0);
		visit(*finest_level, f, u);
	}
	/**
	 * @brief Get the finest Level object
	 *
	 * @return std::shared_ptr<const Level<D>> the Level
	 */
	std::shared_ptr<const Level<D>> getFinestLevel() const
	{
		return finest_level;
	}
};
extern template class Cycle<2>;
extern template class Cycle<3>;
} // namespace GMG
} // namespace ThunderEgg
#endif