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

#ifndef THUNDEREGG_VALVECTORGENERATOR_H
#define THUNDEREGG_VALVECTORGENERATOR_H
#include <ThunderEgg/ValVector.h>
#include <ThunderEgg/VectorGenerator.h>
#include <valarray>
namespace ThunderEgg
{
/**
 * @brief Generates new ValVector objects for a given Domain
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> class ValVectorGenerator : public VectorGenerator<D>
{
	private:
	/**
	 * @brief the Domain to generate ValVector objects for
	 */
	std::shared_ptr<const Domain<D>> domain;
	/**
	 * @brief The number of components in each cell
	 */
	int num_components;

	public:
	/**
	 * @brief Construct a new ValVectorGenerator object
	 *
	 * @param domain the Domain to generate ValVector objects for
	 * @param num_components the number of components for each cell
	 */
	explicit ValVectorGenerator(std::shared_ptr<const Domain<D>> domain, int num_components)
	: domain(domain), num_components(num_components)
	{
	}
	std::shared_ptr<Vector<D>> getNewVector() const override
	{
		return ValVector<D>::GetNewVector(domain, num_components);
	}
};
} // namespace ThunderEgg
#endif
