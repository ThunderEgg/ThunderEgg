/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2017-2021 Scott Aiton
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

#ifndef THUNDEREGG_POISSON_MATRIXHELPER_H
#define THUNDEREGG_POISSON_MATRIXHELPER_H
/**
 * @file
 *
 * @brief MatrixHelper class
 */

#include <ThunderEgg/Domain.h>
#include <petscmat.h>

namespace ThunderEgg::Poisson
{
/**
 * @brief Create a matrix for the 3D second-order Laplacian operator
 *
 * This is equivalent to using StarPatchOperator<3> with TriLinGhostFiller
 */
class MatrixHelper
{
	private:
	/**
	 * @brief the domain
	 */
	Domain<3> domain;
	/**
	 * @brief boundary conditions
	 */
	std::bitset<6> neumann;

	public:
	/**
	 * @brief Create a MatrixHelper for a given 3D domain.
	 *
	 * @param domain the Domain
	 * @param neumann boundary conditions
	 */
	explicit MatrixHelper(const Domain<3> &domain, std::bitset<6> neumann);

	/**
	 * @brief Form the matrix for the domain
	 *
	 * @return the formed matrix
	 */
	Mat formCRSMatrix();
};
} // namespace ThunderEgg::Poisson
#endif
