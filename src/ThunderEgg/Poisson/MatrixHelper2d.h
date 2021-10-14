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

#ifndef THUNDEREGG_POISSON_MATRIXHELPER2D_H
#define THUNDEREGG_POISSON_MATRIXHELPER2D_H
/**
 * @file
 *
 * @brief MatrixHelper2D class
 */

#include <ThunderEgg/Domain.h>
#include <petscmat.h>

namespace ThunderEgg::Poisson
{
/**
 * @brief Create a matrix for the 2D second-order Laplacian operator
 */
class MatrixHelper2d
{
	private:
	Domain<2>      domain;
	std::bitset<4> neumann;

	public:
	/**
	 * @brief Create a MatrixHelper for a given domain
	 *
	 * @param domain the domain
	 * @param neumann the boundary conditions
	 */
	MatrixHelper2d(const Domain<2> &domain, std::bitset<4> neumann);

	/**
	 * @brief Form the matrix
	 *
	 * @return the formed matrix
	 */
	Mat formCRSMatrix(double lambda = 0);
};
} // namespace ThunderEgg::Poisson
#endif
