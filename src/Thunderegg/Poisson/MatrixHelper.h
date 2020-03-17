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

#ifndef THUNDEREGG_POISSON_MATRIXHELPER_H
#define THUNDEREGG_POISSON_MATRIXHELPER_H

#include <Thunderegg/Domain.h>
#include <Thunderegg/PetscVector.h>
#include <petscmat.h>

namespace Thunderegg
{
namespace Poisson
{
/**
 * @brief Create a matrix for the 3D second-order Laplacian operator
 */
class MatrixHelper
{
	private:
	/**
	 * @brief the domain
	 */
	std::shared_ptr<Domain<3>> domain;

	public:
	/**
	 * @brief Create a MatrixHelper for a given 3D domain.
	 *
	 * @param domain the Domain
	 */
	MatrixHelper(std::shared_ptr<Domain<3>> domain);

	/**
	 * @brief Form the matrix for the domain
	 *
	 * @return the formed matrix
	 */
	PW_explicit<Mat> formCRSMatrix(double lambda = 0);
};
} // namespace Poisson
} // namespace Thunderegg
#endif
