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

#ifndef THUNDEREGG_SCHUR_VALVECTORGENERATOR_H
#define THUNDEREGG_SCHUR_VALVECTORGENERATOR_H
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Schur/InterfaceDomain.h>
#include <ThunderEgg/ValVector.h>
#include <ThunderEgg/VectorGenerator.h>
namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief Generates new ValVector objects for a given InterfaceDomain
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> class ValVectorGenerator : public VectorGenerator<D>
{
	private:
	/**
	 * @brief The dimensions of an interface
	 */
	std::array<int, D> iface_ns;
	/**
	 * @brief The InterfaceDomain
	 */
	std::shared_ptr<InterfaceDomain<D + 1>> iface_domain;

	public:
	/**
	 * @brief Construct a new ValVectorGenerator object
	 *
	 * @param iface_domain the InterfaceDomain to generate ValVector objects for
	 */
	ValVectorGenerator(std::shared_ptr<InterfaceDomain<D + 1>> iface_domain)
	: iface_domain(iface_domain)
	{
		std::array<int, D + 1> ns = iface_domain->getDomain()->getNs();
		for (int i = 1; i < D + 1; i++) {
			if (ns[0] != ns[i]) {
				throw RuntimeError(
				"Cannot form Schur compliment vector for Domain with non-square patches");
			}
		}

		iface_ns.fill(ns[0]);
	}
	std::shared_ptr<Vector<D>> getNewVector() const override
	{
		return std::make_shared<ValVector<D>>(MPI_COMM_WORLD, iface_ns, 0,
		                                      iface_domain->getNumLocalInterfaces());
	}
};
} // namespace Schur
} // namespace ThunderEgg
#endif
